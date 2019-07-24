import numpy as np
import subprocess
import sys
import time
sys.path.append('../../scripts/')
from singlerun import OneRun

if __name__=="__main__":

    start_time = time.time()


    FAILED_E = 1e300

    k24,Lambda,omega = sys.argv[1],sys.argv[2],str(20.0)


    scan = {}
    scan['k_{24}'] = k24

    scan['\\Lambda']= Lambda
    scan['\\omega']= omega

    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]

    savesuf=["K_{33}","k_{24}","\\Lambda","\\omega"]

    gammas = np.linspace(0.03,0.05,num=51,endpoint=True)

    i = len(gammas)-1

    while (i >= 0):

        scan['\\gamma_s'] = str(gammas[i])

        # create a class to do calculations with current parameters in scan.
        run = OneRun(scan=scan,loadsuf=loadsuf,savesuf=savesuf,
                     scan_dir="scanbackward")
        # run C executable.
        run.run_exe()

        # move file written by C executable from temporary data path to true data path
        run.mv_file('observables')

        # load the final values of E, R, eta, delta, and surface twist.
        Ei,Ri,R_ci,etai,deltai,surftwisti = run.get_all_observables('observables',str2float=True)
        

        if (Ei > 0.1*FAILED_E):

            # if the energy calculation fails, this will be true.
            print('hi')
            # remove current file with observables for the current Lambda value that are higher than
            # the delta = 0 energy.
            print(Ei)
            run.remove_file("observables")

            break

        if np.isnan(Ri) or Ri <= 0:

            # if Ri is infinite, then the calculation failed.
            # Retry it with a different initial guess.

            print("Ri is NAN, trying again with Rguess = 1.0")

            # remove the current observables file, so that a new one can be written.
            run.remove_file("observables")
            if np.abs(float(scan['Rguess'])-1.0)>1e-10:
                Ri = 1.0
            else:
                break

        else:
            # calculation ran smoothly.
            run.concatenate_observables("\\gamma_s")
            i -= 1
        
        Rguess,R_cguess,etaguess,deltaguess = str(Ri),str(R_ci),str(etai),str(deltai)


        scan['Rguess'] = Rguess
        scan['Rupper'] = str(1.2*float(Rguess))
        scan['Rlower'] = str(0.75*float(Rguess))

        scan['R_cguess'] = R_cguess
        scan['R_cupper'] = str(1.2*float(R_cguess))
        scan['R_clower'] = str(0.75*float(R_cguess))

        if not np.isnan(float(etaguess)):
            scan['etaguess'] = etaguess
            scan['etaupper'] = str(float(etaguess)+0.1)
            scan['etalower'] = str(float(etaguess)-0.02)

        if not (np.isnan(float(deltaguess))
                or abs(float(deltaguess))<1e-5):
            scan['deltaguess'] = deltaguess
            scan['deltaupper'] = '0.818'
            
            if float(deltaguess) < 0.81:
                scan['deltalower'] = str(0.95*float(deltaguess))
            else:
                scan['deltalower'] = '0.81'

    print(f"Took {(time.time()-start_time)/3600} hours to complete.")

