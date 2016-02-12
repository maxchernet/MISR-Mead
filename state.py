#!/usr/bin/env python

"""
The eoldas_ng state class
"""

__author__  = "J Gomez-Dans"
__version__ = "1.0 (1.12.2013)"
__email__   = "j.gomez-dans@ucl.ac.uk"

import cPickle
import platform
from collections import OrderedDict
import time


import numpy as np
import scipy.optimize
import scipy.sparse as sp

from eoldas_utils import *

FIXED = 1
CONSTANT = 2
VARIABLE = 3


class State ( object ):
    
    """A state-definition class
    
    
       In EO-LDAS, the state requires the following:
        
       1. a configuration dictionary,
       2. a state grid
       3. a dictionary with default parameter values
        
       The state grid is an array that defines the domain of the
       problem. For example, for a temporal DA problem, it will be
       a vector of however many elements are required timesteps. For
       a spatial problem, it will be a 2D array with the locations
       of the gridcells. the configuration diciontary stores whether a
       particular parameter is variable over the grid (i.e. is this
       parameter estimated for all timesteps or grid positions?),
       constant (so it is constant in time and/or space), or whether we
       just prescribe some default value."""
       
    def __init__ ( self, state_config, state_grid, default_values, \
            parameter_min, parameter_max, optimisation_options=None, \
            output_name=None, verbose=False ):
        """State constructor
        
        
        """
        #!!!Max Edit
        import os
        c_dir = os.getcwd()
        print 'current dir:', c_dir
        self.f_obj = open(c_dir + '/object_func.dat', 'w')
        #!!!!!!!!!!!!!!!!!!!!!!!

        self.state_config = state_config
        self.state_grid = state_grid
        self.n_elems =  self.state_grid.size
        self.default_values = default_values
        self.operators = {}
        self.n_params = self._state_vector_size ()
        self.parameter_min = parameter_min
        self.parameter_max = parameter_max
        self.verbose = verbose
        self.bounds = []
        for ( i, param ) in enumerate ( self.state_config.iterkeys() ):
            self.bounds.append ( [ self.parameter_min[param], \
                self.parameter_max[param] ] )
        self.invtransformation_dict = {}
        self.transformation_dict = {}
        if optimisation_options is None:
            self.optimisation_options = {"factr": 1000, \
                "m":400, "pgtol":1e-12, "maxcor":200, \
                "maxiter":1500, "disp":True }
        else:
            self.optimisation_options = optimisation_options
        if output_name is None:
            tag = time.strftime( "%04Y%02m%02d_%02H%02M%02S_", time.localtime())
            tag += platform.node()
            self.output_name = "eoldas_retval_%s.pkl" % tag
            
        else:
            self.output_name = output_name
        print "Saving results to %s" % self.output_name
        
    def set_transformations ( self, transformation_dict, \
            invtransformation_dict ):
        """We can set transformations to the data that will be
        applied automatically when required."""
        self.transformation_dict = transformation_dict
        self.invtransformation_dict = invtransformation_dict
        # Recalculate boundaries
        self.bounds = []
        for ( i, param ) in enumerate ( self.state_config.iterkeys() ):
            if transformation_dict.has_key ( param ):
                tmin = transformation_dict[param] ( self.parameter_min[param] )
                tmax = transformation_dict[param] ( self.parameter_max[param] )
            else:
                tmin = self.parameter_min[param]
                tmax = self.parameter_max[param]
            if tmin > tmax:
                self.bounds.append ([ tmax, tmin ] )
            else:
                self.bounds.append ([ tmin, tmax ] )
                
                
    def _state_vector_size ( self ):
        n_params = 0
        for param, typo in self.state_config.iteritems():
            if typo == CONSTANT:
                n_params  += 1
            elif typo == VARIABLE:
                n_params  += self.n_elems
        return n_params
        
    def pack_from_dict ( self, x_dict, do_transform=False ):
        the_vector = np.zeros ( self.n_params )
        # Now, populate said vector in the right order
        # looping over state_config *should* preserve the order
        i = 0
        for param, typo in self.state_config.iteritems():
            if typo == CONSTANT: # Constant value for all times
                if do_transform and self.transformation_dict.has_key ( param ):
                    the_vector[i] = self.transformation_dict[param] ( \
                        x_dict[param] )
                else:
                    the_vector[i] = x_dict[param]
                i = i+1        
            elif typo == VARIABLE:
                # For this particular date, the relevant parameter is at location iloc
                if do_transform and self.transformation_dict.has_key ( param ):
                    the_vector[i:(i + self.n_elems)] =  \
                        self.transformation_dict[param] ( \
                            x_dict[param].flatten() )
                else:
                    the_vector[i:(i + self.n_elems)] =  \
                        x_dict[param].flatten() 
                i += self.n_elems
        return the_vector 
    
    def _unpack_to_dict ( self, x, do_transform=False, do_invtransform=False ):
        """Unpacks an optimisation vector `x` to a working dict"""
        the_dict = OrderedDict()
        i = 0
        for param, typo in self.state_config.iteritems():
            
            if typo == FIXED: # Default value for all times                                
                if self.transformation_dict.has_key ( \
                        param ):
                    the_dict[param] = self.transformation_dict[param]( \
                         self.default_values[param] )
                else:
                #elif do_transform and self.transformation_dict.has_key ( \
                        #param ):
                    
                    the_dict[param] = self.default_values[param] 
                #else:
                    #the_dict[param] = self.default_values[param]
                
            elif typo == CONSTANT: # Constant value for all times
                if do_transform and self.transformation_dict.has_key ( \
                        param ):
                    the_dict[param] = self.transformation_dict[param]( x[i] )
                elif do_invtransform and self.invtransformation_dict.has_key ( \
                        param ):
                    p1 = self.transformation_dict[param] ( self.parameter_max[param] )
                    p2 = self.transformation_dict[param] ( self.parameter_min[param] )
                    pmax = max ( p1, p2 )
                    pmin = min ( p1, p2 )
                    xx = np.clip ( x[i], p1, p2 )
                    the_dict[param] = self.invtransformation_dict[param]( x[i] )

                else:
                    the_dict[param] = x[i]
                i += 1
                
            elif typo == VARIABLE:
                if do_transform and self.transformation_dict.has_key ( \
                        param ):
                    the_dict[param] = self.transformation_dict[param] ( \
                        x[i:(i+self.n_elems )]).reshape( \
                        self.state_grid.shape )
                elif do_invtransform and self.invtransformation_dict.has_key ( \
                        param ):
                    p1 = self.transformation_dict[param] ( self.parameter_max[param] )
                    p2 = self.transformation_dict[param] ( self.parameter_min[param] )
                    pmax = max ( p1, p2 )
                    pmin = min ( p1, p2 )
                    xx = np.clip ( x[i:(i+self.n_elems )].reshape( \
                        self.state_grid.shape), p1, p2 )
                    

                    the_dict[param] = self.invtransformation_dict[param] ( xx )

                else:
                    the_dict[param] = x[i:(i+self.n_elems )].reshape( \
                        self.state_grid.shape )
                i += self.n_elems
            
        return the_dict
    
    def add_operator ( self, op_name, op ):
         """Add operators to the state class
         
         This method will add operator classes (e.g. objects with a `der_cost` and a
         `der_der_cost` method)"""
         the_op = getattr( op, "der_cost", None)
         if not callable(the_op):
             raise AttributeError, "%s does not have a der_cost method!" % op_name     
         self.operators[ op_name ] = op
     
    def _get_bounds_list ( self ):
        the_bounds = []
        for i, ( param, typo ) in enumerate(self.state_config.iteritems()):
            if typo == CONSTANT:
                the_bounds.append ( self.bounds[i] )
            elif typo == VARIABLE:
                
                [ the_bounds.append ( self.bounds[i] ) \
                    for j in xrange ( self.n_elems )]
        return the_bounds
    
    def optimize ( self, x0=None, bounds=None, do_unc=False ):
        
        """Optimise the state starting from a first guess `x0`"""
        start_time = time.clock()
        
        if type(x0) == type ( {} ):
            # We get a starting dictionary, just use that
            x0 = self.pack_from_dict ( x0, do_transform=True )
        elif x0 is None:
            # No starting point, start from random location?
            raise NotImplementedError
        elif type( x0 ) is str:
            # Use a single operator that has a ``first_guess`` method
            x0 = self.operators[x0].first_guess( self.state_config, self.state_grid.size )
            #x0 = self._unpack_to_dict ( x0 )
            ####for param, ptype in self.state_config.iteritems():
                ####if ptype == CONSTANT:
                    ####if not x0.has_key ( param ):
                        ####try:
                            ####x0[param] = self.operators ['Prior'].mu[param]
                        ####except KeyError:
                            ####x0[param] = self.operators ['prior'].mu[param]
                        
            ####x0 = self.pack_from_dict ( x0, do_transform=False )
        #####elif type( x0 ) is list:
            ###### We get one or more operators that provide one with a first guess,
            ###### and hopefully, uncertainty. The idea is to combine these first
            ###### guess estimates
            #####guesses = []
            ###### We loop over the operators we wish to use, and we pass the do_unc
            ###### keyword so that they return a mean + uncertainty
            #####for op in x0:
                #####try:
                    #####guesses.append ( self.operators[x0].first_guess( \
                        #####self.state_config, do_unc=True ) )
                #####except NotImplementedError:
                    ###### No code!
                    #####pass
            #####if len( guesses ) == 1:
                #####x0 = guesses[0][0] # Effectively, same as usual
                #####for param, ptype in self.state_config.iteritems():
                    #####if ptype == CONSTANT:
                        #####if not x0.has_key ( param ):
                            #####try:
                                #####x0[param] = self.operators ['Prior'].mu[param]
                            #####except KeyError:
                                #####x0[param] = self.operators ['prior'].mu[param]
            ######elif len( guesses ) == 2:
                ######for param, ptype in self.state_config.iteritems():
                    ######if ptype == CONSTANT:
                        ######if not x0.has_key ( param ):
                            ######try:
                                ######x0[param] = self.operators ['Prior'].mu[param]
                            ######except KeyError:
                                ######x0[param] = self.operators ['prior'].mu[param]
                    ######elif ptype == VARIABLE:
                        ######if not x0.has_key ( param ):
                            ######try:
                                ######x0[param] = self.operators ['Prior'].mu[param]
                            ######except KeyError:
                                ######x0[param] = self.operators ['prior'].mu[param]
                        
                
                ####### x0 = (x1*s2 + x2*s1)/(s1+s2)
                ######x0 = ( guesses[0][0]*guesses[1][1] + \
                    ######guesses[1][0]*guesses[0][1] ) / \
                    ######( guesses[0][1] + guesses[1][1] )
            ######for param, ptype in self.state_config.iteritems():
                ######if ptype == CONSTANT:
                    ######if not x0.has_key ( param ):
                        ######try:
                            ######x0[param] = self.operators ['Prior'].mu[param]
                        ######except KeyError:
                            ######x0[param] = self.operators ['prior'].mu[param]
                        
            #####x0 = self.pack_from_dict ( x0, do_transform=False )
        if bounds is None:
            the_bounds = self._get_bounds_list()
            
#            r = scipy.optimize.fmin_l_bfgs_b( self.cost, x0, m=100, disp=1, \
#                 factr=1e-3, maxfun=1500, pgtol=1e-20, bounds=the_bounds )
            r = scipy.optimize.minimize ( self.cost, x0, method="L-BFGS-B", \
                jac=True, bounds=the_bounds, options=self.optimisation_options)
            end_time = time.time()
            if self.verbose:
                if r.success:
                    print "Minimisation was successful: %d \n%s" % \
                        ( r.status, r.message )
                else:
                    print "Minimisation was NOT successful: %d \n%s" % \
                        ( r.status, r.message )
                print "Number of iterations: %d" % r.nit
                print "Number of function evaluations: %d " % r.nfev
                print "Value of the function @ minimum: %e" % r.fun
                print "Total optimisation time: %.2f (sec)" % ( time.time() - start_time )
        else:
            r = scipy.optimize.minimize ( self.cost, x0, method="L-BFGS-B", \
                jac=True, bounds=the_bounds, options=self.optimisation_options )
        retval_dict = {}
        retval_dict['real_map'] = self._unpack_to_dict ( r.x, do_invtransform=True )
        retval_dict['transformed_map'] = self._unpack_to_dict ( r.x, \
            do_invtransform=False )
        if do_unc:
            retval_dict.update ( self.do_uncertainty ( r.x ) )
        if self.verbose:
            print "Saving results to %s" % self.output_name
            cPickle.dump ( retval_dict, open( self.output_name, 'wb' ) )

        #!!! Max Edit !!!
        self.f_obj.close()
        import os
        c_dir = os.getcwd()
        f_optim = open(c_dir + '/optimization_results.dat', 'w')
        f_optim.write('%d %d %.4f %.2f'%(r.nit, r.nfev, r.fun, ( time.time() - start_time )))
        f_optim.close()
        #!!!!!!!!!!!!!!!!!!!!!!!

        return retval_dict
    
    def do_uncertainty ( self, x ):
        
        the_hessian = sp.lil_matrix ( ( x.size, x.size ) )
        x_dict = self._unpack_to_dict ( x )
        #cost, der_cost = self.operators["Obs"].der_cost ( x_dict, \
            #self.state_config )
        #this_hessian = self.operators["Obs"].der_der_cost ( x_dict, \
                        #self.state_config, self, epsilon=1e-10 )
        
        #for epsilon in [ 10e-10, 1e-8, 1e-6, 1e-10, 1e-12, ]:
            # print "Hessian with epsilon=%e" % epsilon
        epsilon = 1e-5
        for op_name, the_op in self.operators.iteritems():
            try:
                this_hessian = the_op.der_der_cost ( x, self.state_config, \
                    self, epsilon=epsilon )
            except:
                this_hessian = the_op.der_der_cost ( x_dict, \
                    self.state_config, self, epsilon=epsilon )
            if self.verbose:
                print "Saving Hessian to %s_%s.pkl" % ( self.output_name, \
                    op_name )
            cPickle.dump ( this_hessian, open( "%s_%s_hessian.pkl" \
                % ( self.output_name, op_name ), 'w'))
            the_hessian = the_hessian + this_hessian
        a_sps = sp.csc_matrix( the_hessian )

        lu_obj = sp.linalg.splu( a_sps )
        
        main_diag = np.zeros_like ( x )
        for k in xrange(x.size):
            b = np.zeros_like ( x )
            b[k] = 1
            main_diag[k] = lu_obj.solve ( b )[k]
            
        post_cov = sp.dia_matrix(main_diag,0 ).tolil() # Sparse purely diagonal covariance matrix 
        post_sigma = np.sqrt ( np.abs(main_diag) ).squeeze()
        
        ci_5 = self._unpack_to_dict( x - 1.96*post_sigma, do_invtransform=True )
        ci_95 = self._unpack_to_dict( x + 1.96*post_sigma, do_invtransform=True )
        ci_25 = self._unpack_to_dict( x - 0.67*post_sigma, do_invtransform=True )
        ci_75 = self._unpack_to_dict( x + 0.67*post_sigma, do_invtransform=True )
        retval = {}
        retval['post_cov'] = post_cov
        retval['real_ci5pc'] = ci_5
        retval['real_ci95pc'] = ci_95
        retval['real_ci25pc'] = ci_25
        retval['real_ci75pc'] = ci_75
        retval['post_sigma'] = post_sigma
        
        return retval
        
    def cost ( self, x ):
         """Calculate the cost function using a flattened state vector representation"""
         
         x_dict = self._unpack_to_dict ( x )
         # Store the parameter dictionary in case we need it later for e.g.
         # crossvalidation
         self.parameter_dictionary = x_dict
         aggr_cost = 0
         aggr_der_cost = x*0.0
         self.cost_components = {}
         start_time = time.time()
         
#         import os
#         c_dir = os.getcwd()
#         f_obj = open(c_dir + '/object_func.dat', 'w')

         for op_name, the_op in self.operators.iteritems():
             
             cost, der_cost = the_op.der_cost ( x_dict, self.state_config )
             #!!!Max Edit!!!
             #print 'cost before weighting = ', cost
             cost = cost / the_op.cost_weight
             der_cost = der_cost / the_op.cost_weight
             self.f_obj.write('%.4f '%(cost))
             #!!!!!!!!!!!!!!
             aggr_cost = aggr_cost + cost
             aggr_der_cost = aggr_der_cost + der_cost
             self.cost_components[op_name] = der_cost
             if self.verbose:
                 print "\t%s %8.3e" % ( op_name, cost )
         self.f_obj.write('%.4f\n'%(aggr_cost))
         self.the_cost = aggr_cost


         
         if self.verbose:
             print "Total cost: %8.3e" % aggr_cost
             print 'Elapsed: %.2f seconds' % (time.time() - start_time)
             
             
         return aggr_cost, aggr_der_cost
