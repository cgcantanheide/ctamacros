# --- Imports ---------------------------------------------------------------- #
from ROOT import gROOT,gSystem,TFile,TGraphAsymmErrors,TH1D,TF1,TString,TH2D, TGraph, Double, TSpline3
import cta
from os.path import *
# compile upon every load:
gROOT.SetMacroPath("{0:s}".format(dirname(globals()['cta'].__file__)))
gROOT.LoadMacro("makeCTAspec_v6_pyROOT.C+")
#gSystem.AddIncludePath("{0:s}".format(dirname(globals()['cta'].__file__)))
#gSystem.AddLinkedLibs("{0:s}".format(join(dirname(globals()['cta'].__file__),"makeCTAspec_v6_pyROOT_C.so")))
gSystem.Load("makeCTAspec_v6_pyROOT_C.so")
from ROOT import makeCTAspec
from ctypes import POINTER, c_float, c_double, c_int
from array import array
from cta.utils.spectra import *
from cta.utils.convertRoot2Py import *
import numpy as np
# ---------------------------------------------------------------------------- #

class CTAObsSim(object):
    def __init__(self, irf, **kwargs):
	"""
	Init the observation simulation class

	Parameters
	----------
	irf:	str, full path to CTA IRF root file

	kwargs
	------
	threshold:	float, energy threshold in TeV (default: 0.1)
	eMin:		float, minimum energy in TeV (default: 0.05)
	eMin:		float, maximum energy in TeV (default: 30.)
	"""
	kwargs.setdefault('threshold',0.1)
	#kwargs.setdefault('specBinningOut',100)
	kwargs.setdefault('eMin',0.05)
	kwargs.setdefault('eMax',30.)

	self.__dict__.update(kwargs)

	if not exists(irf):
	    raise IOError('IRF file not found at {0:s}'.format(irf))

	self.threshold		= Double(self.threshold)
	self.irf		= irf
	self.irfFile		= TFile(irf)
	self.specgraph		= TGraphAsymmErrors(1)
	self.ifluxgraph		= TGraphAsymmErrors(1)	# integral flux and its error
	self.excessgraph	= TGraphAsymmErrors(1)	# distribution of excess events
	self.gammaExp		= TH1D()	# distribution of excess events
	self.bkgExp		= TH1D()	# distribution of excess events
	self.r80VsEn		= TGraph(1)	# PSF: 80% containment radius
	self.spline		= TSpline3()
	self.useSpline		= False
	self.xMin		= 1e16		# dummy, very big
	self.xMax		= 1e-16		# dummy, very small

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
	    return H1D2Py(self.irfFile.Get(histname))
	elif kind == 'TGraph':
	    return GraphAsymmErrors2Py(self.irfFile.Get(histname))
	else:
	    raise ValueError(
		"kind keyword must either be 'TH1D' or 'TGraph' not {0:s}".format(kind))
	return


    def setObsSpec(self, nbins = 30):
	"""
	Initialize the observed spectrum

	kwargs
	------
	nbins:	int, number of bins (default : 30)
	"""
	self.lowEdges 		= np.logspace(np.log10(self.eMin), 
					np.log10(self.eMax), nbins + 1)
	# expected spectrum to be 
	# observed with CTA in dN / dE (1/TeV/cm^2/s)
	self.spObserved		= TH1D("spObserved","",nbins, array('d', self.lowEdges))
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
	self.spline	= TSpline3(title,x,y,x.shape[0])
	self.useSpline	= True
	self.xMin	= np.min(x)
	self.xMax	= np.max(x)
	return



    def setIntrSpec(self, params, eMin = 0.01, eMax = 200.):
    	"""
	Set the intrinsic source spectrum

	Parameters
	----------
	params:		dict, function parameters:
			func:	str, function identifier, either 'pl','lp','bpl',or 'epl'
			The remaining kwargs depend on the chosen function.

	kwargs
	------
	eMin:	float, minimum energy in TeV (default: 0.01)
	eMax:	float, maximum energy in TeV (default: 200.)
	"""

	spec, self.func 	= set_root_spec(params)
	self.spIntr		= TF1("origSpec",spec,eMin, eMax)

	return 

    def makeCTAspec(self,**kwargs):
	"""
	Run the simulation

	kwargs
	------
	rebin:		bool, 
			should rebinning be applied if the number 
			of excess events too low? (default: False)
	enres:		bool, 
			should the energy resolution be applied to the data? 
			(default: False) Note: Deprecated, energy resolution always applied
	verbose:	bool, 
			print output? (default: False)
	size_deg:	float,
			radius of the emission in degrees (default: 0.01)
	effOnTime:	float, 
			observation time in seconds (default: 20 hours)
	useExtended:	bool, 
			flag between point-like and extended source 
			(default: False, i.e. point-like)
	useRandom:	bool, 
			flag between BG from spectrum and fixed BG in events 
			(default: True)
	alpha:		float, 
			fraction of exposure between OFF and ON regions (default: 0.2)
	minEvt:		float, 
			minimum number of events per bin (default: 7)
	minSig:		float, 
			minimum significance per bin in sigma (default: 3)
	minBkg:		float, 
			minimum sigmal above background in per cent (default: 0.03)
	minAeff:	float, 
			minimum effective area in cm^2 (default: 1e4)
	"""
	kwargs.setdefault('rebin',False)
	kwargs.setdefault('enres',False)
	kwargs.setdefault('verbose',False)
	kwargs.setdefault('size_deg',0.01)
	kwargs.setdefault('effOnTime',20. * 60. * 60.)
	kwargs.setdefault('useExtended',False)
	kwargs.setdefault('useRandom',True)
	kwargs.setdefault('alpha',0.2)
	kwargs.setdefault('minEvt',7)
	kwargs.setdefault('minSig',3.)
	kwargs.setdefault('minBkg',0.03)
	kwargs.setdefault('minAeff',1e4)

	self.threshold	= makeCTAspec(self.spObserved, 	# histogram to store 
                                                        # the output spectrum in dN / dE
			    self.specgraph,		
			    self.ifluxgraph,  # TGraphAsymmErrors integral flux graph
			    self.excessgraph, # TGraphAsymmErrors excess events graph
			    self.gammaExp,    # TH1D expected gamma rays, 
			   		      # if useRandom == True, 
					      # they will be shuffled with Poisson statistic
			    self.bkgExp,      # TH1D expected bkg rate, 
			   		      # if useRandom == True, 
					      # they will be shuffled with Poisson statistic
			    self.irf,	      # string with IRF root file
			    self.r80VsEn,     # TGraph with 80% containment radius
			    self.spline,      # spline with attenuation
			    self.xMin, 	      # min energy of spline
			    self.xMax,	      # max energy of spline
			    kwargs['rebin'],		
			    kwargs['enres'],
			    self.useSpline,
			    kwargs['verbose'],
			    self.spIntr,      # the intrinsic spectrum graph 
			                      #(will be integrated in each bin)
			    kwargs['effOnTime'],
			    kwargs['size_deg'],
			    POINTER(c_float)(c_float(self.threshold)),
	#		    self.specBinningOut,
			    kwargs['useExtended'],
			    kwargs['useRandom'],
			    kwargs['alpha'],
			    kwargs['minEvt'],
			    kwargs['minBkg'],
			    kwargs['minAeff'],
			    kwargs['minSig']
			 )

	self.specGr	= GraphAsymmErrors2Py(self.specgraph)
	self.spec	= H1D2Py(self.spObserved)
	self.excess	= GraphAsymmErrors2Py(self.excessgraph,xerr = False)
	self.gamma	= H1D2Py(self.gammaExp)
	self.bkg	= H1D2Py(self.bkgExp)
	return
