import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from psidata import PsiData
from observabledata import ObservableData
from readparams import ReadParams

class FibrilStrain(ReadParams):

    def __init__(self,scan_dir="",scan={},
                 loadsuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","\\Lambda",
                          "\\omega","\\gamma_s"],
                 sfile_format="pdf",obsname = "observables",
                 datfile="data/input.dat",psiname = "psivsr",
                 psiloadfilepath="data",obsloadfilepath="data",
                 psisavefilepath="results",obssavefilepath="results"):
        
        ReadParams.__init__(self,datfile=datfile,
                            scan=scan,loadsuf=loadsuf,savesuf=savesuf)

        self.obsdata = ObservableData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,
                                      datfile=datfile,loadfilepath=obsloadfilepath,
                                      savefilepath=obssavefilepath,
                                      scan_dir=scan_dir,name=obsname)

        self.psidata = PsiData(scan=scan,loadsuf=loadsuf,savesuf=savesuf,
                                datfile=datfile,loadfilepath=psiloadfilepath,
                                savefilepath=obssavefilepath,
                                scan_dir = scan_dir,name=psiname)

        self.scan_dir=scan_dir
        self.sfile_format=sfile_format
        self.i_R_c = self.find_R_c_index()
        return

    def find_R_c_index(self):

        i = np.argmin(np.abs(self.psidata.r()-self.obsdata.R_c()))


        return i

    def mesh_polar(self,num_azm=50,grid_skip=1):

        azm = np.linspace(0,2*np.pi,num=num_azm)

        rs,thetas = np.meshgrid(self.psidata.r()[::grid_skip],
                                azm)

        return rs,thetas

    def strain_1d(self,denom='d(r)'):

        preferred_dband = np.cos(self.psidata.psi()[self.i_R_c:])

        true_dband = 2*np.pi/self.obsdata.eta()


        if denom=='d(r)':
            dn = preferred_dband
        elif denom=='d':
            dn = true_dband
        else:
            raise ValueError(f'setting denominator of the strain equation to '
                             '"{denom}" is not valid, it must either be "d(r)" '
                             '(the default value) or "d".')

        dums = np.full(self.i_R_c,np.nan)

        s0s = (true_dband-preferred_dband)/dn

        return np.concatenate((dums,s0s)) 

    def strain_polar(self,r_mesh,denom='d(r)',grid_skip=1):

        return np.tile(self.strain_1d(denom=denom)[::grid_skip],
                       (r_mesh.shape[0],1))


    def strain_sname(self,descriptor="polar"):

        suffix = self.write_suffix(suffix_type="save")

        if self.scan_dir != "":
            sname = f"results/_{descriptor}_strain_{self.scan_dir}_{suffix}.{self.sfile_format}"
        else:
            sname = f"results/_{descriptor}_strain_{suffix}.{self.sfile_format}"

        return sname
