import numpy as np
import subprocess
import sys
import time
sys.path.append('../../scripts/')
from singlerun import SingleRun
from readparams import ReadParams

if __name__=="__main__":
    
    power2 = int(sys.argv[1])
    
    gamma = str(0.04)
    k24 = str(0.5)
    Lambda = str(600.0)
    omega = str(20.0)

    scan = {}
    scan['k_{24}'] = k24
    scan['\\gamma_s'] = gamma
    scan['\\Lambda']=Lambda
    scan['\\omega']= omega
    scan['mpt'] = str(2**power2)
    
    loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"]


    # read in file name info
    rp = ReadParams(scan=scan,loadsuf=loadsuf,savesuf=loadsuf)
        
    # create a class to do calculations with current parameters in scan.
    run = SingleRun(rp,executable="../../../bin/EvsR_c-different-grid-sizes",valgrind=False)
    # run C executable.

    start = time.time()
    
    run.run_exe()

    totaltime = time.time()-start



    with open("data/timings.txt","a") as myfile:
        myfile.write(f"{scan['mpt']}\t{totaltime:e}\n")

    # move file written by C executable from temporary data path to true data path
    run.mv_file(f"EvsR_c-{scan['mpt']}")



