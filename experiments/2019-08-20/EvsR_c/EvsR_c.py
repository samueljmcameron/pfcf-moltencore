import numpy as np
import subprocess
import sys
import time
sys.path.append('../../scripts/')
from singlerun import OneRun
from readparams import ReadParams
from psidata import PsiData
from scipy.integrate import simps

def integrand(rs,psis):

    return 2*rs*psis/(rs[-1]*rs[-1])

    

if __name__=="__main__":



    start_time = time.time()

    if len(sys.argv)<5:

        user_input = input("input string of a gamma,k24,omega values, "
                           "using comma as delimiter: ")
        gamma,k24,Lambda,omega = user_input.split(',')

    else:

        gamma,k24,Lambda,omega = sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]

    
    FAILED_E = 1e300

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\omega']= omega
    scan['\\Lambda']= Lambda
    
    loadsuf=savesuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    scan_dir = "scanforward"

    # first, load the minimum for delta = 0 case, so you know the upper bound for
    # the energy minimum.

    run = OneRun(scan=scan,loadsuf=loadsuf,savesuf=savesuf)


    R_cs = np.linspace(0.02,0.05,num=50,endpoint=True)

        
    for i,R_c in enumerate(R_cs):


            
        executable = "../../../bin/Retadelta3var_onerun"
        scan['R_cguess'] = str(R_c)

        

        # create a class to do calculations with current parameters in scan.
        run = OneRun(scan=scan,loadsuf=loadsuf,savesuf=savesuf,
                     scan_dir=scan_dir,executable=executable,valgrind=True)

        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file(f'observables')

        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,R_ci,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)

        run.concatenate_observables(None,externalparam=[R_c])

        # now just adjust my guess for delta
    
        if np.abs(deltai)<1e-5:
            break


        scan['deltaguess'] = str(deltai)
        scan['deltalower'] = str(deltai*0.95)
        scan['deltaupper'] = '0.817'
        scan['etaguess'] = str(etai)
        scan['etalower'] = str(etai-0.01)
        scan['etaupper'] = str(etai+0.02)
        scan['Rguess'] = str(Ri)
        scan['Rlower'] = str(0.8*Ri)
        scan['Rupper'] = str(1.1*Ri)


    print(f"Took {(time.time()-start_time)/3600} hours to complete.")
