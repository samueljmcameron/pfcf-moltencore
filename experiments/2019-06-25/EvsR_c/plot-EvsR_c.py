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
    obs = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)
        
    R_c = obs.data[:,0]
    E = obs.data[:,1]
    dEdR_c = obs.data[:,2]
    err = obs.data[:,3]

    
    fig = plt.figure()
    fig.set_size_inches(width,3*height)
    

    ax1 = fig.add_subplot(3,1,1)

    ax2 = fig.add_subplot(3,1,2)

    ax3 = fig.add_subplot(3,1,3)

    ax1.plot(R_c,E,'.')
    ax1.set_ylim(-4,4)
    ax1.set_ylabel(r'$E$')
    ax2.plot(R_c,dEdR_c,'.',label='qromb')
    ax2.plot(R_c,np.gradient(E,R_c),label='py deriv')
    ax2.legend(frameon=False)
    ax2.set_ylabel(r'$\frac{\partial E}{\partial R_c}$')
    ax2.set_ylim(-1,dEdR_c[-1])
    ax3.plot(R_c,err*1e4,'.')
    ax3.set_ylabel('derivative error (' + r'$\times \num{1e4}$' + ')')
    ax3.set_xlabel(r'$R_c$')
    fig.subplots_adjust(left=0.2)

    fig.savefig(obs.observable_sname('EvsR_c'))


    plt.show()
