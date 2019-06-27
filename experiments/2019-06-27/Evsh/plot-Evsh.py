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

    
    fig = plt.figure()
    fig.set_size_inches(width,3*height)
    

    ax1 = fig.add_subplot(4,1,1)

    ax2 = fig.add_subplot(4,1,2)

    ax3 = fig.add_subplot(4,1,3)

    ax4 = fig.add_subplot(4,1,4)


    # load in data from timings
    time_data = np.loadtxt("data/timings.txt")
    mpts = time_data[:,0].astype(int)
    times = time_data[:,1]

    R_cs = np.empty([2],float)

    Es = np.empty([len(mpts),2],float)

    dEdR_cs = np.empty([len(mpts),2],float)

    errs = np.empty([len(mpts),2],float)


    for i,mpt in enumerate(mpts):

        scan['mpt'] = mpt

        # read in file name info
        obs = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                             name=f"EvsR_c-{mpt}")


        R_cs[0] = obs.data[0,0]
        R_cs[1] = obs.data[1,0]

        Es[i,0] = obs.data[0,1]
        Es[i,1] = obs.data[1,1]

        dEdR_cs[i,0] = obs.data[0,2]
        dEdR_cs[i,1] = obs.data[1,2]

        errs[i,0] = obs.data[0,3]
        errs[i,1] = obs.data[1,3]




    for i,R_c in enumerate(R_cs[:-1]):

        ax1.plot(mpts,Es[:,i],'.',label=rf'$R_c={R_c}$')

        ax2.plot(mpts,dEdR_cs[:,i],'.',label=rf'$R_c={R_c}$')

        ax3.plot(mpts,errs[:,i]*1e4,'.',label=rf'$R_c={R_c}$')

    ax4.plot(mpts,times,'.')


    ax1.set_ylabel(r'$E(R_c)$')

    ax1.set_xscale('log')

    ax2.set_xscale('log')
    ax2.set_ylabel(r'$\frac{\partial E}{\partial R_c}$')

    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylabel('derivative error (' + r'$\times \num{1e4}$' + ')')

    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set_xlabel('number of points')
    ax4.set_ylabel('time taken')

    fig.subplots_adjust(left=0.3)

    fig.savefig(obs.observable_sname('EvsR_c'))
