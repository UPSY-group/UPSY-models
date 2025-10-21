import os
import glob
import xarray as xr

class Run(object):
    """ Properties and functions from a UFEMISM or LADDIE run """

    def __init__(self, rundir):
        """ Gather basic info from run """

        if rundir[-1] == '/':
            rundir = rundir[:-1]
        self.dir = rundir
        self.name = os.path.basename(self.dir)
        self.model = None
        self.fnames = []
        self.Nmeshes = 0

        self._detect_model()
        self._get_meshes()
        self._get_variables()
        self._get_times()

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

    def _get_meshes(self):
        """ Extract the number of meshes in this run """

        self.fnames = sorted(glob.glob(f'{self.dir}/{self.prefix}_0*.nc'))
        self.Nmeshes = len(self.fnames)

        if self.Nmeshes == 0:
            raise ValueError(f"No valid meshes output files in {self.dir}")

    def _get_variables(self):
        """ Extract available variables in mesh output files """

        ds = xr.open_dataset(self.fnames[0])
        self.vars_vi = [var_name for var_name, var in ds.items() if set(var.dims) == {'time', 'vi'}]
        self.vars_ti = [var_name for var_name, var in ds.items() if set(var.dims) == {'time', 'ti'}]
        self.contours = [var_name for var_name, var in ds.items() if set(var.dims) == {'time', 'two', 'ei'}]
        ds.close()

    def _get_times(self):
        """ Get time values per mesh """

        self.times = {}
        for f,fname in enumerate(self.fnames):
            ds = xr.open_dataset(fname)
            self.times[f+1] = ds.time.values
            ds.close()
            
