from __future__ import print_function
import os, sys
from functools import partial
from time import sleep
import numpy as np

def _obj_wrapper(func, args, kwargs, x):
    return func(x, *args, **kwargs)

def _is_feasible_wrapper(func, x):
    return np.all(func(x)>=0)

def _cons_none_wrapper(x):
    return np.array([0])

def _cons_ieqcons_wrapper(ieqcons, args, kwargs, x):
    return np.array([y(x, *args, **kwargs) for y in ieqcons])

def _cons_f_ieqcons_wrapper(f_ieqcons, args, kwargs, x):
    return np.array(f_ieqcons(x, *args, **kwargs))
    
def pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={}, 
        swarmsize=200, omega=0.5, phip=2.5, phig=1.5, maxiter=100, 
        minstep=1e-6, minfunc=1e-7, debug=False,
        particle_output=False):
    """
    Perform a particle swarm optimization (PSO)
   
    Parameters
    ==========
    func : function
        The function to be minimized
    lb : array
        The lower bounds of the design variable(s)
    ub : array
        The upper bounds of the design variable(s)
   
    Optional
    ========
    ieqcons : list
        A list of functions of length n such that ieqcons[j](x,*args) >= 0.0 in 
        a successfully optimized problem (Default: [])
    f_ieqcons : function
        Returns a 1-D array in which each element must be greater or equal 
        to 0.0 in a successfully optimized problem. If f_ieqcons is specified, 
        ieqcons is ignored (Default: None)
    args : tuple
        Additional arguments passed to objective and constraint functions
        (Default: empty tuple)
    kwargs : dict
        Additional keyword arguments passed to objective and constraint 
        functions (Default: empty dict)
    swarmsize : int
        The number of particles in the swarm (Default: 100)
    omega : scalar
        Particle velocity scaling factor (Default: 0.5)
    phip : scalar
        Scaling factor to search away from the particle's best known position
        (Default: 0.5)
    phig : scalar
        Scaling factor to search away from the swarm's best known position
        (Default: 0.5)
    maxiter : int
        The maximum number of iterations for the swarm to search (Default: 100)
    minstep : scalar
        The minimum stepsize of swarm's best position before the search
        terminates (Default: 1e-8)
    minfunc : scalar
        The minimum change of swarm's best objective value before the search
        terminates (Default: 1e-8)
    debug : boolean
        If True, progress statements will be displayed every iteration
        (Default: False)
    particle_output : boolean
        Whether to include the best per-particle position and the objective
        values at those.
   
    Returns
    =======
    g : array
        The swarm's best known position (optimal design)
    f : scalar
        The objective value at ``g``
    p : array
        The best known position per particle
    pf: arrray
        The objective values at each position in p
   
    """
   
    assert len(lb)==len(ub), 'Lower- and upper-bounds must be the same length'
    assert hasattr(func, '__call__'), 'Invalid function handle'
    lb = np.array(lb)
    ub = np.array(ub)
    assert np.all(ub>lb), 'All upper-bound values must be greater than lower-bound values'
   
    vhigh = np.abs(ub - lb)
    vlow = -vhigh

    # Initialize objective function
    obj = partial(_obj_wrapper, func, args, kwargs)
    
    # Check for constraint function(s) #########################################
    if f_ieqcons is None:
        if not len(ieqcons):
            if debug:
                print('No constraints given.')
            cons = _cons_none_wrapper
        else:
            if debug:
                print('Converting ieqcons to a single constraint function')
            cons = partial(_cons_ieqcons_wrapper, ieqcons, args, kwargs)
    else:
        if debug:
            print('Single constraint function given in f_ieqcons')
        cons = partial(_cons_f_ieqcons_wrapper, f_ieqcons, args, kwargs)
    is_feasible = partial(_is_feasible_wrapper, cons)

    # Initialize the particle swarm ############################################
    S = swarmsize
    D = len(lb)  # the number of dimensions each particle has
    x = np.random.rand(S, D)  # particle positions
    v = np.zeros_like(x)  # particle velocities
    p = np.zeros_like(x)  # best particle positions
    fx = np.zeros(S)  # current particle function values
    fs = np.zeros(S, dtype=bool)  # feasibility of each particle
    fp = np.ones(S)*np.inf  # best particle function values
    g = []  # best swarm position
    fg = np.inf  # best swarm position starting value
    
    # Initialize the particle's position
    x = lb + x*(ub - lb)

    # Calculate objective and constraints for each particle
    # get a first guess for the global minimum by running though all particles once
    ## JRB: NOTE this assumes that for calls to func:
    ##      1) if a run for x doesn't exist, it schedules it
    ##      2) returns None if the answer isn't done calculating yet
    complete = False
    while not complete:
        complete = True
        for i in range(S):
            fx[i] = obj(x[i, :])
            fs[i] = is_feasible(x[i, :])

        if np.any(np.isnan(fx)):
            complete = False
            sleep(120)
       
    # Store particle's best position (if constraints are satisfied)
    i_update = np.logical_and((fx < fp), fs)
    p[i_update, :] = x[i_update, :].copy()
    fp[i_update] = fx[i_update]

    # Update swarm's best position
    i_min = np.argmin(fp)
    if fp[i_min] < fg:
        fg = fp[i_min]
        g = p[i_min, :].copy()
    else:
        # At the start, there may not be any feasible starting point, so just
        # give it a temporary "best" point since it's likely to change
        g = x[0, :].copy()
       
    # Initialize the particle's velocity
    v = vlow + np.random.rand(S, D)*(vhigh - vlow)
       
    # Iterate until termination criterion met ##################################
    it = 1.0
    report_it = 1.0
    while it <= maxiter:
        for i in range(S):
            # get function value
            fx[i] = obj(x[i, :])
            fs[i] = is_feasible(x[i, :])

            # if fx is NaN, it is still running
            if np.isnan(fx[i]):
                continue

            ## otherwise we update things and go to the next (fractional) iteration
            it += 1.0/S

            # Store particle's best position (if constraints are satisfied)
            if (fx[i] < fp[i]) and fs[i]:
                p[i, :] = x[i, :].copy()
                fp[i] = fx[i]

            # Compare particle's best position with global best position
            if fp[i] < fg:
                if debug:
                    print('New best for swarm at iteration {:}: {:} {:}'\
                        .format(it, p[i, :], fp[i]))
    
                p_min = p[i, :].copy()
                stepsize = np.sqrt(np.sum((g - p_min)**2))
    
                if np.abs(fg - fp[i]) <= minfunc:
                    print('Stopping search: Swarm best objective change less than {:}'\
                        .format(minfunc))
                    if particle_output:
                        return p_min, fp[i], p, fp
                    else:
                        return p_min, fp[i]
                elif stepsize <= minstep:
                    print('Stopping search: Swarm best position change less than {:}'\
                        .format(minstep))
                    if particle_output:
                        return p_min, fp[i], p, fp
                    else:
                        return p_min, fp[i]
                else:
                    g = p_min.copy()
                    fg = fp[i]

            ## set new particle position
            # generate random numbers
            rp = np.random.uniform(size=(1, D))
            rg = np.random.uniform(size=(1, D))

            # Update the particle velocity
            v[i,:] = omega*v[i,:] + phip*rp*(p[i,:] - x[i,:]) + phig*rg*(g - x[i,:])
            # Update the particles' positions
            x[i,:] = x[i,:] + v[i,:]
            # Correct for bound violations
            maskl = x[i,:] < lb
            masku = x[i,:] > ub
            x[i,:] = x[i,:]*(~np.logical_or(maskl, masku)) + lb*maskl + ub*masku


        if debug and it >= report_it:
            print('Best after iteration {:}: {:} {:}'.format(it, g, fg))
            sys.stdout.flush()
            report_it += 1.0

        sleep(120)

    print('Stopping search: maximum iterations reached --> {:}'.format(maxiter))
    sys.stdout.flush()
    
    if not is_feasible(g):
        print("However, the optimization couldn't find a feasible design. Sorry")
        sys.stdout.flush()
    if particle_output:
        return g, fg, p, fp
    else:
        return g, fg
