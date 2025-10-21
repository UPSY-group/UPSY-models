import os
import glob
import xarray as xr

class Run(object):
    """ Properties and functions from a UFEMISM or LADDIE run """

    def __init__(self, rundir):
        """ Gather basic info from run """

        self.dir = rundir
        self.name = os.path.basename(self.dir)
        self.model = None
        self.Nmeshes = 0

        self._detect_model()
        self._get_Nmeshes()

        print(f"Loaded {self.__repr__()}")

    def __str__(self):
        return f"Run('{self.dir}')"

    def __repr__(self):
        """ Spit out info on this run """
        modstr = f"\033[32m{self.model}\033[0m"
        namestr = f"\033[36m{self.name}\033[0m"
        Nmeshstr = f"\033[33m{self.Nmeshes}\033[0m"
        if self.Nmeshes == 1:
            return f"{modstr} run {namestr} with {Nmeshstr} mesh"
        else:
            return f"{modstr} run {namestr} with {Nmeshstr} meshes"

    def _detect_model(self):
        """ Detect whether this is a UFEMISM or LADDIE run """

        for region in ['ANT','GRL','NAM']:
            fname = os.path.join(self.dir,f'main_output_{region}_00001.nc')
            if os.path.exists(fname):
                self.model = 'UFEMISM'
                self.prefix = f'main_output_{region}'
        if self.model == None:
            fname = os.path.join(self.dir,'laddie_output_00001.nc')
            if os.path.exists(fname):
                self.model = 'LADDIE'
                self.prefix = f'laddie_output'
            else:
                raise ValueError(f"No valid output files in {self.dir}")

    def _get_Nmeshes(self):
        """ Extract the number of meshes in this run """

        fnames = sorted(glob.glob(f'{self.dir}/{self.prefix}_0*.nc'))
        self.Nmeshes = len(fnames)

        if self.Nmeshes == 0:
            raise ValueError(f"No valid meshes output files in {self.dir}")
