import prosail
from gp_emulator import *
import scipy.stats as ss
from get_soil import *
import pickle 
from Solar import *
import datetime
import os
import fapar_py
import matplotlib.pyplot as plt

#***************************************************************************
def write_fapar_in(f_retval, year, dir_params, dir_fapar_sd, n_train=50, n_ang=7, gamma=0, n_site=1,\
                   lad=1, prior=None):

        # print '***'
        #Find Illumination geometry.........................................

        lat = 41.165
        lon = -96.476
        zone = -6  # time zone (+ to E)
        time=12./24.
        d = datetime
        date = d.date(year, 1, 1)
        theta_i = np.zeros(365) 
        phi_i = np.zeros(365)
        soil_int_std = np.zeros(365)
        lr_int_std = np.zeros(365)
        lt_int_std = np.zeros(365)
        for i in range(365):
            date = date + d.timedelta(days=1)
            theta_i[i], phi_i[i] = CalcSun(lat, lon, date, time, zone)
        theta_i[:] = 60.

        wl_full = np.arange(400, 2501)
        wl_misr = np.array([443, 555, 670, 865])
        # array of time series of soil reflectance
        rs_ts = np.zeros((365, wl_misr.shape[0]))
        rl_ts = np.zeros((365, wl_misr.shape[0]))
        tl_ts = np.zeros((365, wl_misr.shape[0]))

        # Load dictionary from a pickle file..........................................
        retval = pickle.load(open(f_retval, "rb"))

        # Not transformed uncertainties for all parameters............................

        post_sigma = np.zeros((10, 365))
        for i in range(0, 10):
                post_sigma[i, :] = retval['post_sigma'][i*365:i*365+365]

        # Get distribution for parameters...............................................

        if not os.path.exists(dir_params + 'output_fapar_mean/'):
                os.makedirs(dir_params + 'output_fapar_mean/')
        f_fapar_mean = dir_params + 'output_fapar_mean/fapar_in_Ne%d_%d.dat' % (n_site, year)
        f_mean = open(f_fapar_mean, 'w')

        #day = 200
        f_fapar_sd = []
        for day in range(365):
                # print '***'
                # print 'day:', day
                dist = []
                n_params = 10
                p = 0
                #n_train=50
                key_name=[]
                # Get a distrubution for sampling
                for key in retval['transformed_map'].keys():
                #for k in xrange(n_params):
                    param_sd=0
                    key_name.append(key)
                    if post_sigma[p, day] != 0:
                        if np.isnan(post_sigma[p, day]):
                                print 'NaN: %d'%day
                                param_sd = 1/np.sqrt(np.mean(prior.inv_cov[key]))
                        else:
                                param_sd = post_sigma[p, day]
                        #print post_sigma[p,day]
                        dist.append(ss.norm(loc=retval['transformed_map'][key][day], scale=param_sd))
                        # print retval['transformed_map'][key][day], param_sd
                    p += 1
                # The training dataset is obtained by a LatinHypercube Design
                try:
                        x_train = lhd(dist=dist, size=n_train)
                except ValueError:
                        print 'LHD doesnt work'

                # Run prospect for trained parameters.....................................
                
                #prosail.prospect_5b(xleafn[i], xkab[i], car, scen[i], xkw[i], xkm[i])
                leaf_rt_0 = prosail.prospect_5b(retval['real_map']['xleafn'][day], retval['real_map']['xkab'][day],0,\
                        retval['real_map']['scen'][day], retval['real_map']['xkw'][day], retval['real_map']['xkm'][day])
                #copy leaf trans. and refl. at misr wave length to time series arrays
                for l in range(wl_misr.shape[0]):
                        rl_ts[day, l] = leaf_rt_0[wl_misr[l]-400,0]
                        tl_ts[day, l] = leaf_rt_0[wl_misr[l]-400,1]
                lr_arr = np.zeros((n_train, 2101))
                lt_arr = np.zeros((n_train, 2101))
                for i in range(n_train):
                    leaf_rt = prosail.prospect_5b(x_train[i,7], -100.*np.log(x_train[i,3]), 0,\
                                x_train[i,4], -1/50.*np.log(x_train[i,5]), -1/100.*np.log(x_train[i,6])) 
                    if ~np.any(np.isnan(leaf_rt)):
                        lr_arr[i,:] = leaf_rt[:,0]
                        lt_arr[i,:] = leaf_rt[:,1]
                    else:
                        if i>0:
                                lr_arr[i,:] = lr_arr[i-1,:]
                                lt_arr[i,:] = lt_arr[i-1,:]
                lr_std = np.std(lr_arr, axis=0)
                lt_std = np.std(lt_arr, axis=0)
                #print lr_std
                lr_mean = np.mean(lr_arr, axis=0)
                lt_mean = np.mean(lt_arr, axis=0)
                # print 'leaf rerlectance mean shape:', lr_mean.shape

                #Estimate integral of leaf reflectance and transmittance (mean and std)............................

                #lr.append(np.trapz(rt[0:300,0])/300.)
                lr_mean_a = np.trapz(leaf_rt_0[0:300, 0])/300.
                lr_std_up = leaf_rt_0[: ,0] + lr_std
                lr_std_down = leaf_rt_0[:, 0] - lr_std
                lr_std_up_a = np.trapz(lr_std_up[0:300])/300.
                lr_std_down_a = np.trapz(lr_std_down[0:300])/300.
                #print lr_std_up_a, lr_std_down_a
                #print lr_std_up[0:300]
                #print lr_mean[0:300]
                # print 'leaf reflectance mean integral = %.5f'%lr_mean_a
                # print lr_std_up_a, lr_std_down_a
                lr_std_a = lr_std_up_a - lr_std_down_a
                # print 'leaf reflectance std integral = %.5f'%lr_std_a

                lt_mean_a = np.trapz(leaf_rt_0[0:300,1])/300.
                lt_std_up = leaf_rt_0[:,1] + lt_std
                lt_std_down = leaf_rt_0[:,1] - lt_std
                lt_std_up_a = np.trapz(lt_std_up[0:300])/300.
                lt_std_down_a = np.trapz(lt_std_down[0:300])/300.
                #print lt_std_up_a, lt_std_down_a
                # print 'leaf transmittance mean integral = %.5f'%lt_mean_a
                lt_std_a = lt_std_up_a - lt_std_down_a
                # print 'leaf transmittance std integral = %.5f'%lt_std_a

                #Get soil albedo mean and std integral..................................................

                soil_0 = get_soil(retval['real_map']['xs1'][day], retval['real_map']['xs2'][day])
                #Soil reflectance time series at misr bands
                for l in range(wl_misr.shape[0]):
                        rs_ts[day, l] = soil_0[wl_misr[l]-400]
                soil_arr = np.zeros((n_train, 2101))
                for i in range(n_train):
                    soil_arr[i,:] = get_soil(x_train[i,8], x_train[i,9])
                soil_std = np.std(soil_arr, axis=0)
                soil_mean = np.mean(soil_arr, axis=0)
                soil_mean_a = np.trapz(soil_0[0:300])/300.
                soil_std_up = soil_0 + soil_std
                soil_std_down = soil_0 - soil_std
                soil_std_up_a = np.trapz(soil_std_up[0:300])/300.
                soil_std_down_a = np.trapz(soil_std_down[0:300])/300.
                #print soil_std_up_a, soil_std_down_a
                # print 'soil albedo mean integral = %.5f'%soil_mean_a
                soil_std_a = soil_std_up_a - soil_std_down_a
                # print 'soil albedo std integral = %.5f'%soil_std_a

                #fAPAR: THETA_I,PHI_I,NV,THETA_V,PHI_V,LAD,XRS,XHC,XLAI,RPL,XRL,XTL,BRF.............................

                dist = []
                dist.append(ss.norm(loc=soil_mean_a, scale=soil_std_a))
                # Check if posterior sigma is Not a Number. If yes get value from prior
                if np.isnan(post_sigma[1, day]):
                        param_sd = 1/np.sqrt(np.mean(prior.inv_cov['xhc']))
                else:
                        param_sd = post_sigma[1, day]
                dist.append(ss.norm(loc=retval['transformed_map']['xhc'][day], scale=param_sd))
                # print 'sd xhc=', param_sd

                if np.isnan(post_sigma[0, day]):
                        param_sd = 1/np.sqrt(np.mean(prior.inv_cov['xlai']))
                else:
                        param_sd = post_sigma[0, day]
                dist.append(ss.norm(loc=retval['transformed_map']['xlai'][day], scale=param_sd))
                # print 'sd xlai=', param_sd

                if np.isnan(post_sigma[2, day]):
                        param_sd = 1/np.sqrt(np.mean(prior.inv_cov['rpl']))
                else:
                        param_sd = post_sigma[2, day]
                dist.append(ss.norm(loc=retval['transformed_map']['rpl'][day], scale=param_sd))
                # print 'sd rpl=', param_sd
                # Sometimes lr_std_up_a and lr_std_down_a are almost the same and produce zero as difference.
                if lr_std_a == 0.:
                        lr_std_a = 0.0001
                dist.append(ss.norm(loc=lr_mean_a, scale=lr_std_a))
                # print 'sd leaf refl.=', lr_std_a
                if lt_std_a == 0.:
                        lt_std_a = 0.0001
                dist.append(ss.norm(loc=lt_mean_a, scale=lt_std_a))
                # print 'sd leaf trans.=', lt_std_a
                # The training dataset is obtaiend by a LatinHypercube Design
                try:
                        x_train = lhd(dist=dist, size=n_train)
                except ValueError:
                        print 'LHD doesnt work'
                # print 'x_train.shape:', x_train.shape

                # Write an input file for fapar calculation................................................
                #THETA_I,PHI_I,NV,THETA_V,PHI_V,LAD,XRS,XHC,XLAI,RPL,XRL,XTL,BRF

                theta_v = 80.
                phi_v = 0.
                #lad = lad
                if not os.path.exists(dir_fapar_sd):
                        os.makedirs(dir_fapar_sd)
                f_fapar_sd = np.append(f_fapar_sd, dir_fapar_sd+'fapar_in_Ne%d_%d_%d.dat' % (n_site, year, day+1))
                f_sd = open(f_fapar_sd[-1], 'w')
                for i in range(n_train):
                        f_sd.write('%.4f %.4f %d %.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %.4f\n'\
                                %(theta_i[day], phi_i[day], 1, theta_v, phi_v, lad, \
                                x_train[i,0], x_train[i,1], -2.*np.log(x_train[i,2]), \
                                x_train[i,3], x_train[i,4], x_train[i,5]))
                f_sd.close()
                f_mean.write('%.4f %.4f %d %.4f %.4f %d %.4f %.4f %.4f %.4f %.4f %.4f\n'\
                        %(theta_i[day], phi_i[day], 1, theta_v, phi_v, lad,\
                        soil_mean_a, retval['real_map']['xhc'][day], \
                        retval['real_map']['xlai'][day], retval['real_map']['rpl'][day],\
                        lr_mean_a, lt_mean_a))
                
                soil_int_std[day] = soil_std_a
                lr_int_std[day] = lr_std_a
                lt_int_std[day] = lt_std_a
        f_mean.close()
        #Write time series of the integrals
        np.savez('data/prospect_in_sd_Ne%d_%d'%(n_site, year), soil_int_std, lr_int_std, lt_int_std)
        #Write time series of soil reflectance and leaf reflectance and transmittance
        np.savez('data/soil_leaf_ts_Ne%d_%d'%(n_site, year), rs_ts, rl_ts, tl_ts)

        return f_fapar_mean, f_fapar_sd


#*************************************************************************************************************************************************

if __name__ == "__main__":

        f_fapar_mean = ''
        save_dir = os.path.expanduser('~') + '/s3vt_ng/sza_60/' # '/home/max/s3vt_ng/output_misr_etm/misr_prior_model/'
        save_dir_nad = os.path.expanduser('~') + '/s3vt_ng/sza_60/'
        dir_params = save_dir
        gamma = 10000000
        n_train = 50
        lad = 2
        # prior = pickle.load(open('/home/max/s3vt_ng/output_misr_etm/misr_prior_model/prior_2003_1.pkl', 'rb'))
        # compile fortran code using the system call
        os.system("f2py -c -m fapar_py nadimtools.f nadimbrf.f fapar_py_2.f90")
        for n_site in [3]:
                for year in [2007]:

                        f_fapar_mean = save_dir + 'output_fapar_mean/fapar_in_Ne%d_%d_mean.dat' % (n_site, year)
                        dir_fapar_sd = save_dir+'output_fapar_sd/'
                        f_retval = save_dir+'misr_etm_all_Ne%d_%d_ang7.pkl'%(n_site, year)

                        prior_orig = pickle.load(open(save_dir+'prior_%d_%d.pkl' % (year, n_site), "rb"))

                        write_fapar_in(f_retval, year, dir_params, f_fapar_mean, dir_fapar_sd, n_train=50,\
                                       n_site=n_site, lad=lad, prior=prior_orig)

                        # Call fortran functiooutput_fapar_nadim_meanns for fapar
                        fapar_py.fapar_mean(f_fapar_mean, save_dir_nad +\
                                              'output_fapar_nadim_mean/fapar_out_Ne%d_%d_mean.dat' %\
                                              (n_site, year))
                        for day_sd in range(1, 366):
                                fapar_py_2.fapar_sd(dir_fapar_sd +\
                                                    'fapar_in_Ne%d_%d_%d_ang7.dat'%(n_site, year, day_sd),\
                                                    save_dir_nad +\
                                                    'output_fapar_nadim_sd/fapar_out_Ne%d_%d_%d.dat' %\
                                                    (n_site, year, day_sd))

                                # fapar = np.loadtxt(save_dir + 'output_fapar_nadim_sd/fapar_out_Ne%d_%d_%d.dat' % (n_site, year, day_sd))
                                # plot retrieved FAPAR
                                # plt.plot(-2.*np.log(prior.mu['xlai']), label='prior')
                                # plt.title('year %d' % year)
                                # plt.plot(fapar, label='retrieved year %d' % year)
                                # plt.xlim(1, 366)
                                # plt.legend()
                                # plt.grid()
                                # plt.show()

        print 'fAPAR input is done!!!'
