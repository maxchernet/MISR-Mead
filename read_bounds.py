import numpy as np

#**************************************************
#Read bounds of the parameters from .conf file
def read_bounds(conf_file):
        bounds = np.zeros((10,2))
        f = open(conf_file)
        list_opt = f.read().split('\n')
        for i in range(0,len(list_opt)):
                if '[parameter.x.assoc_bounds]' in list_opt[i]:
                        tmp_lst =  [ss.split('=') for ss in list_opt[i+2:i+12]]
                        j=0
                        for line in tmp_lst:
                                bounds[j,:] =  np.array(line[1].split(',')).astype(float)
                                j+=1
        return bounds
#****************************************************************
