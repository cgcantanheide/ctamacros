# --- Imports ---------------------------------------------------------------- #
from ROOT import gROOT,gSystem,TFile,TGraphAsymmErrors,TH1D,TF1,TString,TH2D, TGraph, Double, TSpline3
import ctamacropy
from os.path import *
# compile upon every load:
gROOT.SetMacroPath("{0:s}".format(dirname(globals()['ctamacropy'].__file__)))
gROOT.LoadMacro("makeCTAspec_v6_pyROOT.C+")
#gSystem.AddIncludePath("{0:s}".format(dirname(globals()['cta'].__file__)))
#gSystem.AddLinkedLibs("{0:s}".format(join(dirname(globals()['cta'].__file__),"makeCTAspec_v6_pyROOT_C.so")))
gSystem.Load("makeCTAspec_v6_pyROOT_C.so")
from ROOT import makeCTAspec
from ctypes import POINTER, c_float, c_double, c_int
from array import array
from ctamacropy import spectra 
from ctamacropy import convertroot2py as cr2py
from copy import deepcopy
import numpy as np
# ---------------------------------------------------------------------------- #

class CTAObsSim(object):
    def __init__(self, irf, **kwargs):
	"""
	Init the observation simulation class

	Parameters
	----------
	irf: str
	    full path to CTA IRF root file

	kwargs
	------
	threshold:	float
	    energy threshold in TeV (default: 0.1)
	eMin:		float
	    minimum energy in TeV (default: 0.05)
	eMin:		float
	    maximum energy in TeV (default: 30.)
	"""
	kwargs.setdefault('threshold',0.1)
	#kwargs.setdefault('specBinningOut',100)
	kwargs.setdefault('eMin',0.05)
	kwargs.setdefault('eMax',30.)

	self.__dict__.update(kwargs)

	if not exists(irf):
	    raise IOError('IRF file not found at {0:s}'.format(irf))

	self._threshold		= Double(self.threshold)
	self._irf		= irf
	self._irfFile		= TFile(irf)
	self._useSpline		= False
	self._xMin		= 1e16		# dummy, very big
	self._xMax		= 1e-16		# dummy, very small
	self._spline		= TSpline3()

	return

    def __init_histograms(self):
	"""Initialize the histograms for the simulation"""
	self._specgraph		= TGraphAsymmErrors(1)
	self._ifluxgraph	= TGraphAsymmErrors(1)	# integral flux and its error
	self._excessgraph	= TGraphAsymmErrors(1)	# distribution of excess events
	self._gammaExp		= TH1D()	# distribution of excess events
	self._bkgExp		= TH1D()	# distribution of excess events
	return

    def getHistIRF(self, histname, kind = "TH1D"):
	"""
	Get a histogram from the irf, convert to python and return it

	Parameters
	----------
	histname:	str, 
			name of the histogram

	kwargs
	------
	kind:	str, 
		kind of the histogram, either TH1D or TGraph for 
		1d histograms or (asymmetric) errors

	Returns
	-------
	The converted histogram or graph as numpy array

	Notes
	-----
	For standard IRF, the standard names are:
	- EffectiveArea (TH1)
	- EffectiveAreaEtrue (TH1)
	- EffectiveArea80 (TH1)
	- EffectiveArea (TH1)
	- AngRes (TH1)
	- AngRes80 (TH1)
	- BGRatePerSqDeg (TH1)
	- BGRate (TH1)
	- ERes (TH1)
	- MigMatrix (TH2F) - not implemented for python conversion
	"""

	if kind == 'TH1D':
	    return cr2py.H1D2Py(self._irfFile.Get(histname))
	elif kind == 'TGraph':
	    return cr2py.GraphAsymmErrors2Py(self._irfFile.Get(histname))
	else:
	    raise ValueError(
		"kind keyword must either be 'TH1D' or 'TGraph' not {0:s}".format(kind))
	return


    def __setObsSpec(self, nbins = 30):
	"""
	Initialize the observed spectrum

	kwargs
	------
	nbins:	int, number of bins (default : 30)
	"""
	self._lowEdges 		= np.logspace(np.log10(self.eMin), 
					np.log10(self.eMax), nbins + 1)
	try:
	    del self._spObserved
	except AttributeError:
	    pass
	# expected spectrum to be 
	# observed with CTA in dN / dE (1/TeV/cm^2/s)
	self._spObserved		= TH1D("spObserved","",nbins, array('d', self._lowEdges))
	return

    def setSpline(self, x, y, title = 'attenuation'):
	"""
	Initialize a spline that is multiplied to the intrinsic spectrum 
	(e.g. for EBL absorption)

	Parameters
	----------
	x:	n-dim np.array with energies in TeV
	y:	n-dim np.array with values of Spline (e.g. EBL attenuation)

	kwargs
	------
	title:	str, title of spline
	"""
	self._spline	= TSpline3(title,x,y,x.shape[0])
	self._useSpline	= True
	self._xMin	= np.min(x)
	self._xMax	= np.max(x)
	return



    def setIntrSpec(self, params, eMin = 0.01, eMax = 200.):
    	"""
	Set the intrinsic source spectrum

	Parameters
	----------
	params:		dict
			dictionary with function parameters. Has to include key
			"func" which is a string with function identifier, either 'pl','lp','bpl',or 'epl',
			given in spectra.py
			The remaining kwargs depend on the chosen function.

	kwargs
	------
	eMin:	float, minimum energy in TeV (default: 0.01)
	eMax:	float, maximum energy in TeV (default: 200.)
	"""

	spec, self.func 	= spectra.set_root_spec(params)
	self._spIntr		= TF1("origSpec",spec,eMin, eMax)

	return 

    def makeCTAspec(self,**kwargs):
	"""
	Run the simulation

	kwargs
	------
	rebin:		bool 
			should rebinning be applied if the number 
			of excess events too low? (default: False)
	enres:		bool 
			should the energy resolution be applied to the data? 
			(default: False) Note: Deprecated, energy resolution always applied
	verbose:	bool 
			print output? (default: False)
	size_deg:	float
			radius of the emission in degrees (default: 0.01)
	effOnTime:	float 
			observation time in seconds (default: 20 hours)
	useExtended:	bool 
			flag between point-like and extended source 
			(default: False, i.e. point-like)
	useRandom:	bool 
			flag between BG from spectrum and fixed BG in events 
			(default: True)
	alpha:		float 
			fraction of exposure between OFF and ON regions (default: 0.2)
	minEvt:		float, 
			minimum number of events per bin (default: 7)
	minSig:		float 
			minimum significance per bin in sigma (default: 3)
	minBkg:		float 
			minimum sigmal above background in per cent (default: 0.03)
	minAeff:	float 
			minimum effective area in cm^2 (default: 1e4)
	seed:		int
			random seed for the simulation
	nbins		int
			number of bins for observed spectrum (default: 30)
	
	Returns
	-------
	`~numpy.ndarray` It has the dimension 5 x n where n is the number of energy bins. The rows are:
	    1. left edge of energy bins in TeV
	    2. center of energy bin (log scale) in TeV
	    3. right edge of energy bins in TeV
	    4. Flux in TeV^-1 cm^-2 s^-1
	    5. Uncertainty in Flux in TeV^-1 cm^-2 s^-1
	"""
	kwargs.setdefault('rebin',False)
	kwargs.setdefault('enres',False)
	kwargs.setdefault('verbose',False)
	kwargs.setdefault('size_deg',0.01)
	kwargs.setdefault('effOnTime',20. * 60. * 60.)
	kwargs.setdefault('useExtended',False)
	kwargs.setdefault('useRandom',True)
	kwargs.setdefault('minEvt',7)
	kwargs.setdefault('minSig',3.)
	kwargs.setdefault('minBkg',0.03)
	kwargs.setdefault('minAeff',1e4)
	kwargs.setdefault('seed',0)
	kwargs.setdefault('alpha',0.2)
	kwargs.setdefault('nbins',30)

	self.__init_histograms()
	self.__setObsSpec(kwargs['nbins'])

	self.threshold	= makeCTAspec(self._spObserved, 	# histogram to store 
                                                        # the output spectrum in dN / dE
			    self._specgraph,		
			    self._ifluxgraph,  # TGraphAsymmErrors integral flux graph
			    self._excessgraph, # TGraphAsymmErrors excess events graph
			    self._gammaExp,    # TH1D expected gamma rays, 
			   		      # if useRandom == True, 
					      # they will be shuffled with Poisson statistic
			    self._bkgExp,      # TH1D expected bkg rate, 
			   		      # if useRandom == True, 
					      # they will be shuffled with Poisson statistic
			    self._irf,	      # string with IRF root file
			    deepcopy(self._spline),      # spline with attenuation
			    self._xMin, 	      # min energy of spline
			    self._xMax,	      # max energy of spline
			    kwargs['rebin'],		
			    kwargs['enres'],
			    self._useSpline,
			    kwargs['verbose'],
			    self._spIntr,      # the intrinsic spectrum graph 
			                      #(will be integrated in each bin)
			    kwargs['effOnTime'],
			    kwargs['size_deg'],
			    POINTER(c_float)(c_float(self.threshold)),
			    kwargs['useExtended'],
			    kwargs['useRandom'],
			    kwargs['alpha'],
			    kwargs['minEvt'],
			    kwargs['minBkg'],
			    kwargs['minAeff'],
			    kwargs['minSig'],
			    int(kwargs['seed'])
			 )

	self.specGr	= cr2py.GraphAsymmErrors2Py(self._specgraph)
	self.spec	= cr2py.H1D2Py(self._spObserved)
	self.excess	= cr2py.GraphAsymmErrors2Py(self._excessgraph,xerr = False)
	self.gamma	= cr2py.H1D2Py(self._gammaExp)
	self.bkg	= cr2py.H1D2Py(self._bkgExp)
	return self.spec
