import numpy as np
import subprocess
import sys
import time
import matplotlib.pyplot as plt
sys.path.append('../../scripts/')
from observabledata import ObservableData
from fig_settings import configure_fig_settings

if __name__=="__main__":

    configure_fig_settings()

    width = 3.37
    height = 3.37
    
    gamma = str(0.04)
    k24 = str(0.5)
    Lambda = str(600.0)
    omega = str(20.0)

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\Lambda']=Lambda
    scan['\\omega']= omega
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


    # read in file name info
    FILErf_pf = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,name="rf_pf")
    FILEsplines = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,name="splines")


    rs = FILErf_pf.data[:,0]
    rf_pfs = FILErf_pf.data[:,1]

    ts = FILEsplines.data[:,0]
    cs = FILEsplines.data[:,1]
    
    fig = plt.figure()
    fig.set_size_inches(width,height)
    

    ax1 = fig.add_subplot(1,1,1)

    ax1.plot(rs,rf_pfs,'.',label='discrete-calc')
    ax1.plot(ts,cs,'-',label='cubic-spline')
    ax1.set_xlabel(r'$r$')
    ax1.set_ylabel(r'$r\times f$')
    ax1.legend(frameon=False)

    fig.savefig(FILErf_pf.observable_sname('rf-vs-r'))
