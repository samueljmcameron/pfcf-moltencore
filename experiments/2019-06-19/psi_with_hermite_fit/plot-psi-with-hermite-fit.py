import numpy as np
import subprocess
import sys
import time
import matplotlib.pyplot as plt
sys.path.append('../../scripts/')
from psidata import PsiData
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
    psidata = PsiData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)

    hermitedata = PsiData(scan=scan,loadsuf=loadsuf,savesuf=loadsuf,
                          name = "hermite-psivsr")


    
    fig = plt.figure()
    fig.set_size_inches(width,3*height)
    

    ax1 = fig.add_subplot(3,1,1)

    ax2 = fig.add_subplot(3,1,2)

    ax3 = fig.add_subplot(3,1,3)

    ax1.plot(psidata.r(),psidata.psi(),'.',label='actual')
    ax1.plot(hermitedata.r(),hermitedata.psi(),'-',label='fit')
    ax1.set_ylabel(r'$\psi(r)$')

    ax2.plot(psidata.r(),psidata.psiprime(),'.',label='actual')
    ax2.plot(hermitedata.r(),hermitedata.psiprime(),'-',label='fit')
    ax2.set_ylabel(r'$\frac{d\psi}{dr}$')

    ax3.plot(psidata.r(),psidata.rf_fibril(),'.',label='actual')
    ax3.plot(hermitedata.r(),hermitedata.rf_fibril(),'-',label='fit')
    ax3.legend(frameon=False)
    ax3.set_ylabel(r'$r\times f_{\mathrm{fibril}}(r)$')
    ax3.set_xlabel(r'$r$')


    fig.subplots_adjust(left=0.2)

    fig.savefig(psidata.psivsr_sname())

    plt.show()
