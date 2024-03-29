# Last updated on December 5th, 2018. I merged the original multirun.py module with
# the singlerun.py module, as it was unnecessary to have two modules which essentially
# do the same thing.


import numpy as np
import subprocess
import sys
import os
from readparams import ReadParams

class OneRun(ReadParams):
    # Perform a single run of calculation for psi(r), R, eta, ..., etc.
    # You can choose which parameters to vary by specifying a dictionary
    # of values "scan" when creating the object.

    # attributes:
    #  datfile - the data file that you read parameters from
    #  scan - a dictionary where you can pre-specify certain parameters (vs
    #         reading them in through the datfile).
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,scan={},tmp_path="../../../tmp_data/",scan_dir="",
                 executable="../../../bin/full4var_onerun",strain=None,
                 loadsuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","\\Lambda","\\omega","\\gamma_s"],
                 datfile="data/input.dat",valgrind=False,valFLAG=None):


        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.tmp_path = tmp_path
        self.executable = executable
        self.scan_dir = scan_dir
        self.strain = strain
        self.valgrind = valgrind
        self.valFLAG = valFLAG

        return

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        if self.valgrind and self.valFLAG != None:
            cmd_args = ["valgrind",self.valFLAG,self.executable]
        elif self.valgrind:
            cmd_args = ["valgrind",self.executable]
        else:
            cmd_args = [self.executable]

        if self.strain == None:
            subprocess.run([*cmd_args,self.tmp_path,
                            *self.params.values()],check=True)
            
        else:
            subprocess.run([*cmd_args,self.tmp_path,
                            *self.params.values(),self.strain],check=True)

        return

    def get_xvals(self,fname='observables',str2float=False):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.write_suffix()
        
        fname = f"_{observables}_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            if str2float:
                E,R,eta,delta,surftwist = map(float,line.split())
            else:
                E,R,eta,delta,surftwist = line.split()
        
        return R,eta,delta

    def get_all_observables(self,fname,str2float=False):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.write_suffix()
        
        fname = f"_{fname}_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            if str2float:
                observables=map(float,line.split())
            else:
                observables=line.split()

        return observables
    
    def mv_file(self,mname,newname=None,strain=None):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.write_suffix()

        if strain != None:

            suffix = suffix + f"_{float(strain):1.4e}"

        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        if newname != None:
            newfname = f"_{newname}_{suffix}.txt"
        else:
            newfname = fname

        mvto = f"data/{newfname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

    def add_datastring(self,vrs,externalparam=None):

        if isinstance(vrs,list):
            a = [float(self.params[var]) for var in vrs]
            dstring = '\t'.join(map("{:13.6e}".format,a))
        elif vrs == None:
            a = externalparam
            if isinstance(externalparam,list):
                dstring = '\t'.join(map("{:13.6e}".format,a))
            else:
                dstring = f"{a:13.6e}"

        else:
            a = float(self.params[vrs])
            dstring = f"{a:13.6e}"

        return dstring

    def remove_file(self,fname="observables",strain=None):

        suffix = self.write_suffix()

        if strain != None:

            suffix = suffix + f"_{float(strain):1.4e}"

        fname = f"data/_{fname}_{suffix}.txt"

        if os.path.isfile(fname):
            os.remove(fname)

        return

    def concatenate_observables(self,vrs,externalparam=None):

        suffix1 = self.write_suffix(suffix_type="save")

        newfname = f"data/_observables_{self.scan_dir}_{suffix1}.txt"

        with open(newfname,"a+") as f1:
    
            suffix2 = self.write_suffix()

            fname = f"data/_observables_{suffix2}.txt"

            with open(fname) as f2:

                for line in f2:

                    f1.write(f"{self.add_datastring(vrs,externalparam=externalparam)}"
                             + f"\t{line}")
            

            if os.path.isfile(fname):

                os.remove(fname)
        
        return

    def write_observables(self,E0,R0,eta0,delta0,
                          surftwist0,vrs,modulus=None,externalparam=None):

        suffix1 = self.write_suffix(suffix_type="save")

        newfname = f"data/_observables_{self.scan_dir}_{suffix1}.txt"

        with open(newfname,"a+") as f1:

            if modulus != None:
            
                line=(f"{E0:13.6e}\t{R0:13.6e}\t{eta0:13.6e}\t"
                      +f"{delta0:13.6e}\t{surftwist0:13.6e}\t{modulus:13.6e}\n")

            else:
                line=(f"{E0:13.6e}\t{R0:13.6e}\t{eta0:13.6e}\t"
                      +f"{delta0:13.6e}\t{surftwist0:13.6e}\n")

            f1.write(f"{self.add_datastring(vrs,externalparam=externalparam)}\t{line}")
            
        suffix2 = self.write_suffix()
        fname = f"data/_observables_{suffix2}.txt"

        if os.path.isfile(fname):
            
            os.remove(fname)
        
        return

################################################################################
# EVERYTHING BELOW THIS POINT IS DEPRECATED AS OF 2019/07/19.  #
################################################################################


class SingleRun(object):
    # Perform a single run of calculation for psi(r), R, eta, ..., etc.
    # You can choose which parameters to vary by specifying a dictionary
    # of values "scan" when creating the object.

    # attributes:
    #  datfile - the data file that you read parameters from
    #  scan - a dictionary where you can pre-specify certain parameters (vs
    #         reading them in through the datfile).
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,readparams,tmp_path="../../../tmp_data/",scan_dir="",
                 params=None,executable="../../../bin/full3var_onerun",strain=None,
                 valgrind=False,valFLAG=None):

        self.readparams = readparams
        self.tmp_path = tmp_path
        self.executable = executable
        self.scan_dir = scan_dir
        self.strain = strain
        self.valgrind = valgrind
        self.valFLAG = valFLAG
        if (params == None):
            self.params = self.readparams.params
        else:
            self.params = params

        return

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        if self.valgrind and self.valFLAG != None:
            cmd_args = ["valgrind",self.valFLAG,self.executable]
        elif self.valgrind:
            cmd_args = ["valgrind",self.executable]
        else:
            cmd_args = [self.executable]

        if self.strain == None:
            subprocess.run([*cmd_args,self.tmp_path,
                            *self.params.values()],check=True)
            
        else:
            subprocess.run([*cmd_args,self.tmp_path,
                            *self.params.values(),self.strain],check=True)

        return

    def get_xvals(self,fname='observables',str2float=False):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.readparams.write_suffix()
        
        fname = f"_{observables}_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            if str2float:
                E,R,eta,delta,surftwist = map(float,line.split())
            else:
                E,R,eta,delta,surftwist = line.split()
        
        return R,eta,delta

    def get_all_observables(self,fname,str2float=False):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.readparams.write_suffix()
        
        fname = f"_{fname}_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            if str2float:
                observables=map(float,line.split())
            else:
                observables=line.split()

        return observables
    
    def mv_file(self,mname,newname=None,strain=None):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.readparams.write_suffix()

        if strain != None:

            suffix = suffix + f"_{float(strain):1.4e}"

        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        if newname != None:
            newfname = f"_{newname}_{suffix}.txt"
        else:
            newfname = fname

        mvto = f"data/{newfname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

    def add_datastring(self,vrs,externalparam=None):

        if isinstance(vrs,list):
            a = [float(self.params[var]) for var in vrs]
            dstring = '\t'.join(map("{:13.6e}".format,a))
        elif vrs == None:
            a = externalparam
            if isinstance(externalparam,list):
                dstring = '\t'.join(map("{:13.6e}".format,a))
            else:
                dstring = f"{a:13.6e}"

        else:
            a = float(self.params[vrs])
            dstring = f"{a:13.6e}"

        return dstring

    def remove_file(self,fname="observables",strain=None):

        suffix = self.readparams.write_suffix()

        if strain != None:

            suffix = suffix + f"_{float(strain):1.4e}"

        fname = f"data/_{fname}_{suffix}.txt"

        if os.path.isfile(fname):
            os.remove(fname)

        return

    def concatenate_observables(self,vrs,externalparam=None):

        suffix1 = self.readparams.write_suffix(suffix_type="save")

        newfname = f"data/_observables_{self.scan_dir}_{suffix1}.txt"

        with open(newfname,"a+") as f1:
    
            suffix2 = self.readparams.write_suffix()

            fname = f"data/_observables_{suffix2}.txt"

            with open(fname) as f2:

                for line in f2:

                    f1.write(f"{self.add_datastring(vrs,externalparam=externalparam)}"
                             + f"\t{line}")
            

            if os.path.isfile(fname):

                os.remove(fname)
        
        return

    def write_observables(self,E0,R0,eta0,delta0,
                          surftwist0,vrs,modulus=None,externalparam=None):

        suffix1 = self.readparams.write_suffix(suffix_type="save")

        newfname = f"data/_observables_{self.scan_dir}_{suffix1}.txt"

        with open(newfname,"a+") as f1:

            if modulus != None:
            
                line=(f"{E0:13.6e}\t{R0:13.6e}\t{eta0:13.6e}\t"
                      +f"{delta0:13.6e}\t{surftwist0:13.6e}\t{modulus:13.6e}\n")

            else:
                line=(f"{E0:13.6e}\t{R0:13.6e}\t{eta0:13.6e}\t"
                      +f"{delta0:13.6e}\t{surftwist0:13.6e}\n")

            f1.write(f"{self.add_datastring(vrs,externalparam=externalparam)}\t{line}")
            
        suffix2 = self.readparams.write_suffix()
        fname = f"data/_observables_{suffix2}.txt"

        if os.path.isfile(fname):
            
            os.remove(fname)
        
        return
