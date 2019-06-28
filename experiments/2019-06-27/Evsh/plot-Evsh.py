import numpy as np
import subprocess
import sys
import time
import matplotlib.pyplot as plt
sys.path.append('../../scripts/')
from observabledata import ObservableData
from fig_settings import configure_fig_settings
import seaborn as sns

if __name__=="__main__":

    configure_fig_settings()

    colors = sns.color_palette()

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
    fig.set_size_inches(width,2*height)
    

    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)

    mpts = np.array([16384,32768,65536],int)#[1024,2048,4096,8192],int)#,

    for i,mpt in enumerate(mpts):

        scan['mpt'] = str(mpt)
        print(i)

        # read in file name info
        obs = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                             name=f"Evsh-{mpt}")

        hs = obs.data[:,0]
        dEdR_cs = obs.data[:,1]
        err_trunc = obs.data[:,2]
        err_round = obs.data[:,3]



        ax1.plot(hs,dEdR_cs,'-',label=rf'$m={mpt}$',color=colors[i])
        ax2.plot(hs,err_trunc,'-',label=rf'$m={mpt}$',color=colors[i])
        ax2.plot(hs,err_round,'-',color=colors[i])


    ax1.set_ylabel(r'$\frac{\partial E}{\partial R_c}$')
    ax1.set_xscale('log')
    ax1.legend(frameon=False)


    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('derivative error (' + r'$\times \num{1e4}$' + ')')

    fig.subplots_adjust(left=0.3)

    fig.savefig(obs.observable_sname('Evsh'))
