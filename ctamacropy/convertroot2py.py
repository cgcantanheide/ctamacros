from ROOT import Double
import numpy as np

def GraphAsymmErrors2Py(graph, xerr = True, yerr = True):
    """
    Extract x and y values from a ROOT::TGraphAsymmErrors object
    (works also with TGraph objects)

    Parameters
    ----------
    graph:	the ROOT::TGraphAsymmErrors object

    kwargs
    ------
    xerr:	bool, extract the asymmetrical errors on x (default: True)
    yerr:	bool, extract the asymmetrical errors on y (default: True)

    Returns
    -------
    list with 2+1+1 np.ndarrays, containing the x,y[,xerr,yerr] values
    """
    nbins	= graph.GetN()
    x, y	= Double(0), Double(0)

    result	= []
    resX, resY	= np.empty(nbins),np.empty(nbins)

    for i in range(nbins):
	graph.GetPoint(i,x,y)
	resX[i] = x
	resY[i] = y
    result.append(resX)
    result.append(resY)

    if xerr:
	resXerr = np.empty((nbins,2))
	for i in range(nbins):
	    resXerr[i,0] = graph.GetErrorXlow(i)
	    resXerr[i,1] = graph.GetErrorXhigh(i)
	result.append(resXerr)

    if yerr:
	resYerr = np.empty((nbins,2))
	for i in range(nbins):
	    resYerr[i,0] = graph.GetErrorYlow(i)
	    resYerr[i,1] = graph.GetErrorYhigh(i)
	result.append(resYerr)

    return result

def H1D2Py(h1d):
    """
    Extract bins and their contents from a  ROOT::TH1D object
    (works also with TGraph objects)

    Parameters
    ----------
    h1d:	the ROOT::TH1D object

    Returns
    -------
    (5 x nbins) x dim np.ndarray, 
    containing the low edges, bin centers, bin width, bin content, and error on bin content
    """
    nbins	= h1d.GetNbinsX()
    result	= np.empty((5,nbins))

    for i in range(nbins):
	result[0,i] = h1d.GetBinLowEdge(i)
	result[1,i] = h1d.GetBinCenter(i)
	result[2,i] = h1d.GetBinWidth(i)
	result[3,i] = h1d.GetBinContent(i)
	result[4,i] = h1d.GetBinError(i)

    return result[:,1:]
