from numpy import exp,log,power,linspace, ones, isscalar
from scipy.integrate import simps

def set_root_spec(params):
    """
    Set the root spectrum

    Parameters 
    -----------
    params:	dict, function parameters:
		func:	str, function identifier, either 'pl','lp','bpl',or 'epl'

		The remaining kwargs depend on the chosen function.

    Returns
    -------
    tuple with str that can be interpreted by root and python function pointer
    """
    if params['func'] == 'const':
	spec		= "{0[norm]}".format(params)
	func		= const
    elif params['func'] == 'lp':
	spec		= "{0[norm]}*pow(x / {0[scale]},-({0[alpha]} + {0[beta]}*log(x/{0[scale]})))".format(params)
	func		= lp
    elif params['func'] == 'pl':
	spec		= "{0[norm]}*pow(x / {0[scale]},{0[index]})".format(params)
	func		= pl
    elif params['func'] == 'epl':
	spec		= "{0[norm]}*pow(x / {0[scale]},{0[index])*exp(x * {0[cut]})".format(params)
	func 		= epl
    elif params['func'] == 'bpl':
	spec		= "{0[norm]}*pow(x / {0[scale]},{0[index1]})*pow(1. + pow(x * {0[break]},{0[f]}),({0[index2]} - ({0[index1]}))/{0[f]})".format(params)
	func 		= bpl
    else:
	raise ValueError
    return spec,func

def const(x,**p):
    """
    Return a constant value

    Paramaters
    ----------
    p: dictionary with parameters:
	    - norm: flux
    x: ndim array 

    Returns
    -------
    n dim array with const values
    p['norm']
    """
    if isscalar(x):
	return p['norm']
    return p['norm']*ones(x.shape[0])

def pl(x,**p):
    """
    Return power law for x with params p

    Paramaters
    ----------
    p: dictionary with parameters:
	    - index: power-law index
	    - scale: pivot energy
	    - norm: flux at x = p['scale']
    x: ndim array with values to evaluate power law

    Returns
    -------
    n dim array with power law values
    p['norm']*(power(x/p['scale'],p['index']))
    """
    return p['norm']*(power(x/p['scale'],p['index']))

def epl(x,**p): 
    """
    Return power law with exponential cut off for x with params p

    Paramaters
    ----------
    p: dictionary with parameters:
	    - index: power-law index
	    - scale: pivot energy
	    - norm: flux at x = p['scale']
	    - cut: 1 / cut parameter
    x: ndim array with values to evaluate power law

    Returns
    -------
    n dim array with power law values
    p['norm']*(power(xp['scale'],p['index']))*exp(x*p['cut'])
    """
    return p['norm']*(power(x/p['scale'],p['index']))*exp(x*p['cut'])

def lp(x,**p):
    """
    Return logarithmic parabola for x with params p

    Paramaters
    ----------
    p: dictionary with parameters:
	    - alpha: power-law index
	    - scale: pivot energy
	    - norm: flux at x = p['scale']
	    - beta: curvature
    x: ndim array with values to evaluate power law

    Returns
    -------
    n dim array with power law values
    p['norm']*power(x/p['scale'], -(p['alpha'] + p['beta'] * log(x/p['scale'])) )
    """
    return p['norm']*power(x/p['scale'], -(p['alpha'] + p['beta'] * log(x/p['scale'])) )

def bpl(x,**p):
    """
    Return broken power law for x with params p

    Paramaters
    ----------
    p: dictionary with parameters:
	    - index1: first power-law index
	    - index2: second power-law index
	    - scale: pivot energy
	    - break: inverse of break 
	    - f: smoothing parameter
	    - norm: flux at x = p['scale']
    x: ndim array with values to evaluate broken power law

    Returns
    -------
    n dim array with broken power law values
    p['norm']*(power(x/p['scale'],p['index1'])) * 
    (1. + (x * p['break']) ** p['f']) ** ((p['index2'] - p['index1'])/p['f'])
    """
    return p['norm']*(power(x/p['scale'],p['index1'])) * \
	    (1. + (x * p['break']) ** p['f']) ** ((p['index2'] - p['index1'])/p['f'])

def crab_hess(x):
    """
    Return flux of crab in TeV^-1 cm^-2 s^-1

    Parameters
    ----------
    x: ndim array, energies in TeV

    Returns
    -------
    ndim array with flux of Crab nebula as given in Aharonian et al. 2006
    """
    p = {'norm': 3.76e-11, 'index':-2.39, 'scale': 1., 'cut': -1./14.3}
    return epl(x,**p)
