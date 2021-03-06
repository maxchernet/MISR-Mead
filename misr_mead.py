#
#
# ****************************************************
#
# Chernetskiy Maxim 2016
# EO-LDAS inversion of ETM+ and MISR together
#
#
#
#
# ****************************************************
#
#
# System modules
import numpy as np
import scipy.stats as ss
from collections import OrderedDict
import os
import time
import pickle
import matplotlib.pyplot as plt

# EO-LDAS modules
from gp_emulator import *
# from eoldas_ng import *
from state import *
from operators import *
import semidiscrete as rt

# My modules
from read_bounds import *
#from first_guess import *
from fapar_7 import *
import julday
import fapar_py

timeBefore1 = time.clock()
# *************************************************************************************
# Redefine the Prior class because I need 365 days prior for LAI "model"
# *************************************************************************************
class PriorMax(Prior):
       def der_cost ( self, x_dict, state_config ):
        """Calculate the cost function and its partial derivatives for the prior object.
        Assumptions of normality are clear# Constant value for all times
        Takes a parameter dictionary, and a state configuration dictionary as
        inputs.

        Parameters
        -----------
        x_dict: ordered dict
            The state as a dictionary
        state_config: oredered dict
            The configuration dictionary

        Returns
        --------
        cost: float
            The value of the cost function
        der_cost: array
            An array with the partial derivatives of the cost function
        """
        if sp.issparse ( self.inv_cov ):
            x = self.pack_from_dict ( x_dict, state_config )
            err = sp.lil_matrix ( x - self.mu )
            cost = err.dot ( self.inv_cov ).dot ( err.T )
            der_cost  = np.array( err.dot ( self.inv_cov ).todense()).squeeze()
            cost = float(np.array(cost.todense()).squeeze())
            return cost, der_cost

        # Find out about problems size
        n, n_elems = get_problem_size ( x_dict, state_config )
        # Allocate space/initialise the outputs
        der_cost = np.zeros ( n )
        cost = 0
        # The next loop calculates the cost and associated partial derivatives
        # Mainly based on the parameter type
        i = 0 # Loop variable

        for param, typo in state_config.iteritems():
            if typo == FIXED:
                # Doesn't do anything so we just skip
                pass
            elif typo == CONSTANT:
                cost = cost + 0.5*( x_dict[param] - \
                            self.mu[param])**2*self.inv_cov[param]
                der_cost[i] = ( x_dict[param] - self.mu[param]) * \
                            self.inv_cov[param]
                i += 1
            elif typo == VARIABLE:

                if self.inv_cov[param].size == 1:
                    # Single parameter for all sites/locations etc
                    # This should really be in the __init__ method!
                    sigma = self.inv_cov[param]

                    self.inv_cov[param] = sp.dia_matrix ( ( np.ones(n_elems)*sigma, 0 ), shape=(n_elems, n_elems))
                # Max Edit
                else:
                    sigma = self.inv_cov[param]
                    self.inv_cov[param] = sp.dia_matrix ( ( np.ones(n_elems)*sigma, 0 ), shape=(n_elems, n_elems))
                # End of Max Edit

                cost_m = ( x_dict[param].flatten() - self.mu[param]) * ( \
                            self.inv_cov[param] )
                cost = cost + 0.5*(cost_m*(x_dict[param].flatten() - \
                            self.mu[param])).sum()
                der_cost[i:(i+n_elems)] = cost_m

                i += n_elems

        return cost, der_cost
# *************************************************************************************
# *************************************************************************************
# Redefine time_step method without 5 grad step
# and change n_obs n_bands in init
class ObservationOperatorMax ( ObservationOperatorTimeSeriesGP ):

    def __init__ ( self, state_grid, state, observations, mask, emulators, bu, \
            band_pass=None, bw=None ):
        """
         observations is an array with n_bands, nt observations. nt has to be the
         same size as state_grid (can have dummny numbers in). mask is nt*4
         (mask, vza, sza, raa) array.


        """
        self.state = state
        self.observations = observations
        try:
            self.n_obs, self.n_bands = self.observations.shape
        except:
            raise ValueError, "Typically, obs should be (n_obs, n_bands)"
        self.mask = mask
        assert ( self.n_obs ) == mask.shape[0]
        self.state_grid = state_grid
        self.nt = self.state_grid.shape[0]
        self.emulators = emulators
        self.bu = bu
        self.band_pass = band_pass
        self.bw = bw
    # ************************************************************************************
    def time_step ( self, this_loc ):
        """Returns relevant information on the observations for a particular time step.
        """
        tag = np.round( self.mask[ this_loc, 2:])
        tag = tuple ( (tag[:2].astype(np.int)).tolist() )
        this_obs = self.observations[ this_loc, :]
        return self.emulators[tag], this_obs, [ self.band_pass, self.bw ]

        # tag = np.round(self.mask[ this_loc, 1:].astype (np.int))
        # tag = tuple ( (tag[:2].astype(np.int)).tolist() )
        # this_obs = self.observations[ this_loc, :]
       	# return self.emulators[tag], this_obs, [ self.band_pass, self.bw ]
# ********************************************************************************************
# Redifine class State to save *.pkl in other folder
# ********************************************************************************************
class StateMax(State):
    def _create_output_file ( self, output_name ):

        self.netcdf = False
        if output_name is None:
            tag = time.strftime( "%04Y%02m%02d_%02H%02M%02S_", time.localtime())
            tag += platform.node()
            self.output_name = "output/retval/eoldas_retval_%s" % tag

        elif isinstance ( output_name, basestring ):
            self.output_name = output_name + ".pkl"
        else:
            self.output_name = output_name.fname
            self.retval_file = output_name
            self.netcdf = True


        print "Saving results to %s" % self.output_name
    #
    # *************************************************************************************************
    #
    def do_uncertainty ( self, x ):
        """A method to calculate the uncertainty. Takes in a state vector.

        Parameters
        -----------
        x: array
            State vector (see ``_self._pack_to_dict``)

        Returns
        ---------
        A dictionary with the values for the posterior covariance function
        (sparse matrix), 5, 25, 75 and 95 credible intervals, and the
        main diagonal standar deviation. In each of these (apart from the
        posterior covariance sparse matrix), we get a new dictionary with
        parameter keys and the parameter estimation represented in the
        selected state grid."""


        the_hessian = sp.lil_matrix ( ( x.size, x.size ) )
        x_dict = self._unpack_to_dict ( x )
        #cost, der_cost = self.operators["Obs"].der_cost ( x_dict, \
            #self.state_config )
        #this_hessian = self.operators["Obs"].der_der_cost ( x_dict, \
                        #self.state_config, self, epsilon=1e-10 )

        #for epsilon in [ 10e-10, 1e-8, 1e-6, 1e-10, 1e-12, ]:
            # print "Hessian with epsilon=%e" % epsilon
        # epsilon is defined in order to use der_der_cost methods that
        # evaluate the Hessian numerically
        epsilon = 1e-8
        for op_name, the_op in self.operators.iteritems():
            # The try statement is here to allow der_der_cost methods to
            # take either a state dictionary or a state vector
            try:
               this_hessian = the_op.der_der_cost ( x_dict, \
                    self.state_config, self, epsilon=epsilon )
	    except OperatorDerDerTypeError:
                this_hessian = the_op.der_der_cost ( x, self.state_config, \
                    self, epsilon=epsilon )
            if self.verbose:
                print "Saving Hessian to %s_%s.pkl" % ( self.output_name, \
                    op_name )
            # Save the individual Hessian contributions to disk
            cPickle.dump ( this_hessian, open("%s_%s_hessian.pkl" \
                % ( self.output_name, op_name ), 'w'))
            # Add the current Hessian contribution to the global Hessian
            the_hessian = the_hessian + this_hessian
        # Need to change the sparse storage format for the Hessian to do
        # the LU decomposition
        a_sps = sp.csc_matrix( the_hessian )
        # LU decomposition object
        lu_obj = sp.linalg.splu( a_sps )
        # Now, invert the Hessian in order to get the main diagonal elements
        # of the inverse Hessian (e.g. the variance)
        main_diag = np.zeros_like ( x )
        for k in xrange(x.size):
            b = np.zeros_like ( x )
            b[k] = 1
            main_diag[k] = lu_obj.solve ( b )[k]

        post_cov = sp.dia_matrix(main_diag,0 ).tolil() # Sparse purely diagonal covariance matrix
        post_sigma = np.sqrt ( main_diag ).squeeze()
        # Calculate credible intervals, transform them back to real units, and
        # store in a dictionary.
        _ci_5 = self._unpack_to_dict( x - 1.96*post_sigma, do_invtransform=True )
        _ci_95 = self._unpack_to_dict( x + 1.96*post_sigma, do_invtransform=True )
        _ci_25 = self._unpack_to_dict( x - 0.67*post_sigma, do_invtransform=True )
        _ci_75 = self._unpack_to_dict( x + 0.67*post_sigma, do_invtransform=True )
        # There intervals are OK in transformed space. However, we need to ensure that
        # e.g. ci_5 <= ci_95 in real coordinates. We do this in the next loop
        ci_5 = {}
        ci_95 = {}
        ci_25 = {}
        ci_75 = {}

        for k in self.state_config.iterkeys():
            ci_5[k]  = np.min(np.array([_ci_5[k],_ci_95[k]]),axis=0)
            ci_95[k] = np.max(np.array([_ci_5[k],_ci_95[k]]),axis=0)
            ci_25[k] = np.min(np.array([_ci_25[k],_ci_75[k]]),axis=0)
            ci_75[k] = np.max(np.array([_ci_25[k],_ci_75[k]]),axis=0)


        # Now store uncertainty, and return it to the user in a dictionary
        retval = {}
        retval['post_cov'] = post_cov
        retval['real_ci5pc'] = ci_5
        retval['real_ci95pc'] = ci_95
        retval['real_ci25pc'] = ci_25
        retval['real_ci75pc'] = ci_75
        retval['post_sigma'] = post_sigma
        retval['hessian'] =  the_hessian
        return retval
#
# *******************************************************************************************
# *******************************************************************************************
#
# Create an instance of the EO-LDAS_ng state class
#
def get_state():
        # Read bounds of the state variables from an old EO-LDAS *.conf file
        bounds = read_bounds('conf/misr_02')
        # Do log. transformation for LAI, Cab, water content and dry matter
        bounds[0, :] = -2.*np.log(bounds[0, :])[::-1]
        bounds[3, :] = -100.*np.log(bounds[3, :])[::-1]
        bounds[5, :] = -(1./50.)*np.log(bounds[5, :])[::-1]
        bounds[6, :] = -(1./100.)*np.log(bounds[6, :])[::-1]
        print bounds[0, :]
        print bounds[3, :]
        print '%.5f %.5f' % (bounds[5, 0], bounds[5, 1])
        print '%.5f %.5f' % (bounds[6, 0], bounds[6, 1])
        # **************************************************************
        # Define state and its parameters
        # ***************************************************************
        # xlai,xhc,rpl,xkab,scen,xkw,xkm,xleafn,xs1,xs2
        state_config = OrderedDict()
        # Define are they fix, constant or variable
        state_config['xlai'] = VARIABLE
        state_config['xhc'] = VARIABLE
        state_config['rpl'] = VARIABLE
        state_config['xkab'] = VARIABLE
        state_config['scen'] = VARIABLE
        state_config['xkw'] = VARIABLE
        state_config['xkm'] = VARIABLE
        state_config['xleafn'] = VARIABLE
        state_config['xs1'] = VARIABLE
        state_config['xs2'] = VARIABLE
        # Set up default values
        default_par = OrderedDict()
        default_par['xlai'] = 2.
        default_par['xhc'] = 1.
        default_par['rpl'] = 0.01
        default_par['xkab'] = 40.
        default_par['scen'] = 0.0
        default_par['xkw'] = (-1/50.)*np.log ( 0.995 )
        default_par['xkm'] = (-1/100.)*np.log ( 0.995 )
        default_par['xleafn'] = 1.5
        default_par['xs1'] = 1.0
        default_par['xs2'] = 0.0
        #
        # Set min and max values of the parameters
        parameter_min = OrderedDict()
        parameter_max = OrderedDict()
        # bounds = read_bounds('/home/max/s3vt/conf/misr')
        min_vals = bounds[:, 0]
        max_vals = bounds[:, 1]
        for i, param in enumerate(state_config.keys()):
            parameter_min[param] = min_vals[i]
            parameter_max[param] = max_vals[i]
        #
        # Define parameter transformations
        transformations = {
            'xlai': lambda x: np.exp(-x/2.),\
            'xkab': lambda x: np.exp(-x/100.),\
            'xkw': lambda x: np.exp( -50.*x),\
            'xkm': lambda x: np.exp(-100.*x)}
        inv_transformations = {
            'xlai': lambda x: -2.*np.log(x),\
            'xkab': lambda x: -100.*np.log(x),\
            'xkw': lambda x: (-1/50.)*np.log(x),\
            'xkm': lambda x: (-1/100.)*np.log(x)}

        # Set state grid in time (year)
        state_grid = np.arange(1, 366)

        # Create class instance
        state = StateMax(state_config, state_grid, default_par, parameter_min, parameter_max, verbose=True)

        # Set the transformations for lai, cab, kw and km
        # I think it would be better to do in class initialization
        state.set_transformations(transformations, inv_transformations)

        return state


# ***********************************************
# Cross validation i.e. finding of optimal gamma
# ***********************************************
def cross_valid(state_file, cross_val=0, is_etm=False):
        # Load a data file
        d = np.loadtxt(state_file, skiprows=1)
        c=3
        if is_etm == False:
                rho = d[c::7, 6:]
        else: rho = d[:, 6:]
        print 'rho.shape:', rho.shape
        # Cross validation if cross_val > 0
        cross = 100/float(100-cross_val)
        # generate an array of integers of size rho
        a = np.arange(rho.shape[0])
        # randomly mix it up
        np.random.shuffle(a)
        # take % of initial size
        ind_cross = a[:rho.shape[0]/cross]
        # sort index incrementaly
        ind_cross = np.sort(ind_cross)
        return ind_cross


# ***************************************************************************
# Define observational operator for MISR
# ***************************************************************************
def get_obs_operator_misr(emu_file, state, year, ind_cross, state_file, cam='An'):
    #
    # Parameters:
    # emu_file - emulator file name which determines by SZA, VAA and RAA
    # state - an instance of state
    # year - an year
    # ind_cross - a vector indices which is used for cross validation
    # state_file -a data file with actual satellite reflectances
    # cam - a MISR camera
    #
    wl_misr = [443, 555, 670, 865]
    wl_width = [21,15,11,20]
    wl_full = np.arange(400,2501)
    # Load a data file
    d = np.loadtxt(state_file, skiprows=1)
    if cam == 'An': c = 3
    if cam == 'Af': c = 2
    if cam == 'Aa': c = 4
    if cam == 'Bf': c = 1
    if cam == 'Ba': c = 5
    if cam == 'Cf': c = 0
    if cam == 'Ca': c = 5
    doys = d[::7, 0].astype(int)
    vza = np.round(d[c::7, 2])
    sza = np.round(d[c::7, 4])
    raa = np.round(abs(d[c::7, 3] - d[c::7, 5]))
    rho = d[c::7, 6:]
    print 'rho.shape:', rho.shape
    #446nm +-21, 558nm +-15, 672nm +-11, 866nm +-20
    b_min = np.array(wl_misr) - np.array(wl_width)
    b_max = np.array(wl_misr) + np.array(wl_width)

    # if cros_vall==0 they will be the same size
    rho = rho[ind_cross, :]
    vza = vza[ind_cross]
    sza = sza[ind_cross]
    raa = raa[ind_cross]
    doys = doys[ind_cross]
    print 'rho.shape:', rho.shape

    n_bands = b_min.shape[0]
    band_pass = np.zeros((n_bands, 2101), dtype=np.bool)
    bw = np.zeros(n_bands)
    bh = np.zeros(n_bands)
    for i in xrange(n_bands):
        band_pass[i, :] = np.logical_and(wl_full >= b_min[i], wl_full <= b_max[i])
        bw[i] = b_max[i] - b_min[i]
        bh[i] = (b_max[i] + b_min[i])/2.

    sigma_min = 0.004
    sigma_max = 0.015
    sigma_obs = (sigma_max - sigma_min)*(bh-bh.min())/(bh.max()-bh.min())
    sigma_obs += sigma_min
    bu = sigma_obs
    print bu
    # Uncertainties of observations
   # bu = [0.0040, 0.0044, 0.0047, 0.0054]
    # Try to change low resolution uncertainty
    # bu = [0.012, 0.0132, 0.0141, 0.0162]
    print bu
    rho_big = np.zeros((365, n_bands))
    mask = np.zeros((365, 4))

    for i in state.state_grid:
        if i in doys:
            r = rho[doys==i, :].squeeze()
            rho_big[i-1,:] = r
            # mask[ i, :] = [ i, raa[doys==i], np.round(sza[doys==i]/5.)*5, np.round(vza[doys==i]/5.)*5 ]
            mask[ i-1, :] = [ i, raa[doys==i], np.round(sza[doys==i]), np.round(vza[doys==i]) ]

    ind = np.where(mask!=0)

    emulators = {}
    for i in range(doys.shape[0]):

        #tag = tuple([ np.round(sza[i]/5.)*5, np.round(vza[i]/5.)*5 ])
        tag = tuple([ np.round(sza[i]), np.round(vza[i]) ])

        file_emu = emu_file%(sza[i], vza[i], raa[i])
        print sza[i], vza[i], raa[i]
        if os.path.isfile(file_emu):
            print file_emu
            print tag
            emulators[tag] = MultivariateEmulator ( dump=file_emu)
        else: print 'emulator file '+file_emu+' does not exist'
    # Create an instance the ObservationOperatorTimeSeriesGP class
    obs = ObservationOperatorMax(state.state_grid, state, rho_big, mask, emulators, bu, band_pass, bw)
    return obs, doys


# ***************************************************************************
# Define observational operator for Landsat
# ***************************************************************************
def get_obs_operator_etm(emu_file, state, year, ind_cross, state_file, cost_weight=1):

    wl_misr = [483.0, 560.0, 662.0, 835.0, 1648.0, 2206.0]
    wl_width=[65, 80, 60, 150, 200, 270]
    wl_full = np.arange(400,2501)
    # Load data a file
    d = np.loadtxt(state_file, skiprows=1)
    d = d[ind_cross, :]

    jul = d[:,0].astype(int)
    dt = datetime
    # Find difference between first day of an year and array of all julian days
    day_diff = jul - julday.date2jul(dt.date(year, 1, 1))
    #Find days for this year
    ind = np.logical_and(day_diff>0, day_diff<366)
    # Find unique days in an year. Landsat data can contain two observations for one day.
    # This is not supported in the new generation...
    (arr, ind_u) = np.unique(jul, return_index=True)
    # array of not unique indices
    ind_z = np.zeros(jul.shape[0], dtype=bool)
    ind_z[:] = True
    # set all unique indices to False
    ind_z[ind_u] = False
    # set all not unique indices to False
    ind[ind_z] = False
    doys = []
    for day1 in jul[ind]:
        doys = np.append(doys, julday.jul2doy(day1))
    print 'doys=', doys

    vza = np.round(d[ind, 2])
    sza = np.round(d[ind, 4])
    vaa = np.round(d[ind, 3])
    saa = np.round(d[ind, 5])
    raa = np.round(abs(vaa - saa))
    rho = d[ind, 6:]
    #print 'vza=', vz
    print 'rho.shape (etm):', rho.shape
    # 446nm +-21, 558nm +-15, 672nm +-11, 866nm +-20
    b_min = np.array(wl_misr) - np.array(wl_width)
    b_max = np.array(wl_misr) + np.array(wl_width)

    n_bands = b_min.shape[0]
    band_pass = np.zeros(( n_bands,2101), dtype=np.bool)
    bw = np.zeros( n_bands )
    bh = np.zeros( n_bands )
    for i in xrange( n_bands ):
        band_pass[i,:] = np.logical_and ( wl_full >= b_min[i], wl_full <= b_max[i] )
        bw[i] = b_max[i] - b_min[i]
        bh[i] = ( b_max[i] + b_min[i] )/2.

    sigma_min = 0.001
    sigma_max = 0.04
    sigma_obs = (sigma_max - sigma_min)*(bh-bh.min())/(bh.max()-bh.min())
    sigma_obs += sigma_min
    bu = sigma_obs
    print bu
    bu = [0.0041, 0.0045, 0.0047, 0.0053, 0.0079, 0.0097]
    print bu
    rho_big = np.zeros(( 365, n_bands ))
    mask = np.zeros(( 365, 4))
    for i in state.state_grid:
        if i in doys:
            r = rho[doys==i, :].squeeze()
            print doys
            rho_big[i-1,:] = r
            #mask[ i, :] = [ i, raa[doys==i], np.round(sza[doys==i]/5.)*5, np.round(vza[doys==i]/5.)*5 ]
            mask[ i-1, :] = [ i, raa[doys==i], np.round(sza[doys==i]), np.round(vza[doys==i]) ]

    ind = np.where(mask!=0)

    emulators = {}
    for i in range(doys.shape[0]):

        tag = tuple([ np.round(sza[i]), np.round(vza[i]) ])

        file_emu = emu_file%(sza[i], vza[i], raa[i])
        print sza[i], vza[i], raa[i]
        if os.path.isfile(file_emu):
            print file_emu
            print tag
            emulators[tag] = MultivariateEmulator ( dump=file_emu)
        else:
                print 'emulator file '+file_emu+' does not exist'
                print sza[i], vza[i], raa[i]
                print vaa[i], saa[i]
    obs = ObservationOperatorMax(state.state_grid, state, rho_big, mask, emulators, bu, band_pass, bw, cost_weight)
    return obs, doys


# ******************************************************************************************************************
# Define prior values of state variables
# ******************************************************************************************************************
def get_prior(state):
        mu_prior = OrderedDict()
        prior_inv_cov = OrderedDict()
        prior_inv_cov['xlai'] = np.array([1.0])
        prior_inv_cov['xhc'] = np.array([1.0])
        prior_inv_cov['rpl'] = np.array([1.0])
        prior_inv_cov['xkab'] = np.array([1.0])
        prior_inv_cov['scen'] = np.array([1.0])
        prior_inv_cov['xkw'] = np.array([1.0])
        prior_inv_cov['xkm'] = np.array([1.0])
        prior_inv_cov['xleafn'] = np.array([1.0])
        prior_inv_cov['xs1'] = np.array([2.0])
        prior_inv_cov['xs2'] = np.array([2.0])
        # Set transformations for some variables
        for param in state.parameter_min.iterkeys():
            if state.transformation_dict.has_key(param):
                mu_prior[param] = state.transformation_dict[param](np.array([state.default_values[param]]))
            else:
                mu_prior[param] = np.array([state.default_values[param]])
            prior_inv_cov[param] = 1./prior_inv_cov[param]**2

        #load lai "phenological model"
        #prior_lai = np.load('/home/max/s3vt_ng/data/fit_ts_mean.npz')['arr_0']
        #mu_prior['xlai'] = state.transformation_dict['xlai'](prior_lai)

        # lai_p = np.load('/home/max/phenology/lai_model_2.npz')['arr_0']
        # lai_sd_p = np.load('/home/max/phenology/lai_model_2.npz')['arr_1']
        # mu_prior['xlai'] = lai_p
        # prior_inv_cov['xlai'] = 1/lai_sd_p**2


        # Define a "phenological model" for LAI
        lai_p = np.load('lai_model_2.npz')['arr_0']
        lai_sd_p = np.load('lai_model_2.npz')['arr_1']
        mu_prior['xlai'] = lai_p
        prior_inv_cov['xlai'] = 1/lai_sd_p**2



        # Create an instance of the prior class and return it
        return PriorMax(mu_prior, prior_inv_cov)


# *******************************************************************
# Do a first guess to define initial values (inversion)
# *******************************************************************
def get_first_guess(state, obs, year, save_dir, n_site=1, etm=True):
        # state - an instance of state
        # obs - an instance of observational operator
        # year - an year
        # save_dir - where to save stuff
        # n_site - site number
        # etm - a flag is it Landsat or MISR
        xlai1 = []
        xhc1 = []
        rpl1 = []
        xkab1 = []
        scen1 = []
        xkw1 = []
        xkm1 = []
        xleafn1 = []
        xs11 = []
        xs21 = []
        doys = []
        if obs.observations[0, :].shape[0]>4:
                wl = np.array([483, 560, 662, 835, 1648, 2206])
        else: wl = np.array([443, 555, 670, 865])
        print 'wl=', wl
        # Make a distribution
        distr_file = os.path.expanduser('~') + '/s3vt/data/US_Ne%d_misr_distr'%1
        # read train dataset from a file
        npzfile = np.load(distr_file+'.npz')
        # load train state parameters
        x_train = npzfile['arr_0']
        # train reflectance
        brf_train = npzfile['arr_1']
        train_misr=[]
        # Reduce number of bands of trained dataset to
        # a number of bands of satellite data
        for i in xrange(wl.size):
            tmp = brf_train[:,wl[i]-400]
            train_misr.append(tmp)
        train_misr = np.array(train_misr)
        n_train = brf_train.shape[0]

        for i in state.state_grid:
                if obs.mask[i-1, 0]!=0:
                        # Find sum of difference
                        sum_misr=[]
                        for j in xrange(n_train):
                                sum_misr.append(np.sum(abs(obs.observations[i-1,:] - train_misr[:,j])))
                        # Find a minimum and this is an end of first guess
                        min_i = np.argmin(sum_misr)
                        # Initial values = first guess
                        xlai1.append(-2.*np.log(x_train[min_i, 0]))
                        xhc1.append(x_train[min_i, 1])
                        rpl1.append(x_train[min_i, 2])
                        xkab1.append(-100*np.log( x_train[min_i, 3]))
                        scen1.append( x_train[min_i, 4] )
                        xkw1.append(-(1./50.)*np.log( x_train[min_i, 5]))
                        xkm1.append(-(1./100.)*np.log( x_train[min_i, 6]))
                        xleafn1.append( x_train[min_i, 7])
                        xs11.append( x_train[min_i, 8])
                        xs21.append( x_train[min_i, 9])
                        doys.append(i)
        x_dict = {}
        print state.state_grid.shape
        print doys
        print len(xlai1)
        x_dict['xlai'] = np.interp(state.state_grid, doys, xlai1)
        x_dict['xhc'] = np.interp(state.state_grid, doys, xhc1)
        x_dict['rpl'] = np.interp(state.state_grid, doys, rpl1)
        x_dict['xkab'] = np.interp(state.state_grid, doys, xkab1)
        x_dict['scen'] = np.interp(state.state_grid, doys, scen1)
        x_dict['xkw'] = np.interp(state.state_grid, doys, xkw1)
        x_dict['xkm'] = np.interp(state.state_grid, doys, xkm1)
        x_dict['xleafn'] = np.interp(state.state_grid, doys, xleafn1)
        x_dict['xs1'] = np.interp(state.state_grid, doys, xs11)
        x_dict['xs2'] = np.interp(state.state_grid, doys, xs21)
        pickle.dump( x_dict, open(save_dir+'first_guess_Ne%d_%d.pkl'%(n_site, year), "wb"))
        return x_dict

# .............................................................................................................................
# ************************************************************************************************************************
# .............................................................................................................................
# US_Ne1_%d_ang7.brf
# data_dir = '/home/max/s3vt/data_misr/'
def run_mead(f_retval, save_dir, misr_dir, etm_dir, emul_dir, n_site=1, year=2001, n_ang=1, lad=2):
    # n_ang = 1
    # Define Leaf Angle Distribution
    # lad = 2
    # A directory for results
    #save_dir = '/home/max/s3vt_ng/output_misr_etm/prior_model_lad5/'
    # Loop by test sites numbers
    #for n_site in [1]:
            # Loop by years
    #        for year in [2001]:
    print 'year: %d, site: %d' % (year, n_site)
    state_file_misr = misr_dir + 'US_Ne%d_%d_ang7.brf' % (n_site, year)
    state_file_etm = etm_dir + 'ls7_us%d.dat' % n_site
    state = get_state()
    # % of taken off data for cross validation
    # If cross_val =0, no cross validation
    cross_val = 0
    ind_cross = cross_valid(state_file_misr, cross_val=cross_val)
    ind_cross_etm = cross_valid(state_file_etm, cross_val=cross_val, is_etm=True)
    # Get observation operators for all the MISR cameras and all sensors (if any)
    print 'An'
    obs_misr_an, doys = get_obs_operator_misr(emul_dir+\
                                                'nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='An')
    print 'Af'
    obs_misr_af, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Af')
    print 'Aa'
    obs_misr_aa, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Aa')
    print 'Bf'
    obs_misr_bf, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Bf')
    print 'Ba'
    obs_misr_ba, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Ba')
    print 'Cf'
    obs_misr_cf, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Cf')
    print 'Ca'
    obs_misr_ca, doys = get_obs_operator_misr(emul_dir+\
                                              '/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
                                                state, year, ind_cross, state_file_misr, cam='Ca')


    # obs_etm, doys = get_obs_operator_etm(emul_dir+'/nad_%03d_sza_%03d_vza_%03d_raa_gp.npz',\
    #        state, year, ind_cross_etm, state_file_etm, cost_weight=1)

    prior = get_prior(state)
    pickle.dump(prior, open(save_dir+'prior_%d_%d.pkl' % (year, n_site), "wb"))

    # Loop by the regularization parameter gamma
    for gamma in [10000000]:
            # ***************************************************************************
            # Define temporal constraint
            # ***************************************************************************
            temporal = TemporalSmoother(state.state_grid, gamma=gamma, order=2, required_params=\
                        ["xlai", "xhc", "rpl", "xkab", "scen", "xkw", "xkm", "xleafn", "xs1", "xs2"])

            # ********************************************************************
            # Add all operators to state
            # ********************************************************************
	    # misr_cam = ["MISR_An", "MISR_Af", "MISR_Aa", "MISR_Bf", "MISR_Ba", "MISR_Cf", "MISR_Ca"]
            state.add_operator("Prior", prior)
            state.add_operator("Temporal", temporal)
            state.add_operator("MISR_An", obs_misr_an)
            state.add_operator("MISR_Af", obs_misr_af)
            state.add_operator("MISR_Aa", obs_misr_aa)
            state.add_operator("MISR_Bf", obs_misr_bf)
            state.add_operator("MISR_Ba", obs_misr_ba)
            state.add_operator("MISR_Cf", obs_misr_cf)
            state.add_operator("MISR_Ca", obs_misr_ca)
            # state.add_operator("ETM", obs_etm )

            # **********************************************************
            # First guess
            # **********************************************************
            # We can estimate starting value without first guess. Random.
#                        x_dict = {}
#                        x_dict['xlai'] = np.random.randn(365)*5.
#                        x_dict['xhc'] = np.random.randn(365)*5.
#                        x_dict['rpl'] = np.random.randn(365)*0.1
#                        x_dict['xkab'] = 40.*np.random.randn(365)
#                        x_dict['scen'] = np.random.randn(365)
#                        x_dict['xkw'] = 0.07*np.random.randn(365)
#                        x_dict['xkm'] = 0.02*np.random.randn(365)
#                        x_dict['xleafn'] = 0.02*np.random.randn(365)*1.5
#                        x_dict['xs1'] = 0.02*np.random.randn(365)*4.0
#                        x_dict['xs2'] = 0.02*np.random.randn(365)*3.0

            # Or with a first guess
            x_dict = get_first_guess(state, obs_misr_an, year, save_dir, n_site=n_site, etm=False)
            # state.output_name = f_retval
            # ******************************************************************
            # Do Optimization
            # ******************************************************************
            timeBefore2 = time.clock()
            retval = state.optimize(x_dict, do_unc=True)
            timeAfter2 = time.clock()
            elapsed_time2 = timeAfter2 - timeBefore2
            print 'optimization elapsed time (s): ', elapsed_time2
            # f_retval = save_dir+'misr_etm_all_Ne%d_%d_ang%d_uncreal.pkl'%(n_site, year, n_ang)
            # Save a dictionary into a pickle file.
            pickle.dump(retval, open(f_retval, "wb"))

    timeAfter1 = time.clock()
    elapsed_time1 = timeAfter1 - timeBefore1
    print 'Total elapsed time (s): ', elapsed_time1
    print 'Done!!!\n'

    return 1


#
# *****************************************************************************************************
#
def do_fapar(save_dir, f_retval, n_site=1, year=2001, lad=1):
    # Write files which will be used by Fortran for calculation of fAPAR
    dir_params = save_dir
    # dir_fapar_mean = save_dir+'output_fapar_mean/'
    # f_fapar_mean = save_dir + 'output_fapar_mean/fapar_in_Ne%d_%d.dat' % (n_site, year)
    dir_fapar_sd = save_dir+'output_fapar_sd/'
    if not os.path.exists(save_dir + 'output_fapar_nadim_mean'):
        os.makedirs(save_dir + 'output_fapar_nadim_mean')
    if not os.path.exists(save_dir + 'output_fapar_nadim_sd'):
        os.makedirs(save_dir + 'output_fapar_nadim_sd')
    prior_orig = pickle.load(open(save_dir+'prior_%d_%d.pkl' % (year, n_site), "rb"))
    try:
        (f_fapar_mean, f_fapar_sd) = write_fapar_in(f_retval, year, dir_params, dir_fapar_sd, n_train=50,\
                           n_site=n_site, lad=lad, prior=prior_orig)
    except ValueError:
            print 'Something wrong with write_fapar_in()...'
    # compile fortran code using the system call
    os.system("f2py -c -m fapar_py nadimtools.f nadimbrf.f fapar_py.f90")
    # Call fortran functions for fapar
    out_fapar_mean = save_dir + 'output_fapar_nadim_mean/fapar_out_Ne%d_%d.dat' % (n_site, year)

    fapar_py.fapar_mean(f_fapar_mean, out_fapar_mean)
    out_fapar_sd = []
    # print 'dir_fapar_sd=', dir_fapar_sd
    # print 'dir_fapar_sd[-1]=', dir_fapar_sd[-1]
    for day_sd in range(1, 366):
        out_fapar_sd = np.append(out_fapar_sd, save_dir + 'output_fapar_nadim_sd/fapar_out_Ne%d_%d_%d.dat' %\
                            (n_site, year, day_sd))
        # print 'out_fapar_sd=', out_fapar_sd
        fapar_py.fapar_mean(dir_fapar_sd + 'fapar_in_Ne%d_%d_%d.dat' % (n_site, year, day_sd), out_fapar_sd[-1])

    return f_fapar_mean, f_fapar_sd, out_fapar_mean, out_fapar_sd


#
# ***********************************************************************************************************
#
# MaX user is ucfamc3
#
#************************************************************************************************************
if __name__ == "__main__":
    save_dir = 'output/'
    misr_dir = 'data_misr/'
    etm_dir = 'data_landsat/'
    n_site = 1
    year = 2002
    lad = 2
    n_ang = 7
    lad = 2
    for n_site in [2, 3]:
        for year in range(2001, 2009):
            f_retval = save_dir + 'misr_etm_all_Ne%d_%d_ang%d.pkl' % (n_site, year, n_ang)
            emul_dir = os.path.expanduser('~') + '/DATA/semidiscrete/lad%d/'%lad
            run_mead(f_retval, save_dir, misr_dir, etm_dir, emul_dir, n_site=n_site, year=year, n_ang=n_ang, lad=lad)

            (in_fapar_mean, in_fapar_sd, out_fapar_mean, out_fapar_sd) = do_fapar(save_dir, f_retval,\
                                                                              n_site=n_site, year=year, lad=lad)