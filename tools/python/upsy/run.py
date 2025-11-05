import os
import glob
import xarray as xr

from upsy.mesh import Mesh

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

    # User functions

    def get_mesh(self,mesh_number):
        """ Open mesh """

        assert mesh_number <= self.Nmeshes, 'Mesh number too high, not available in output'
        assert mesh_number >= 1, 'Mesh number too low, should be at least 1'

        mesh = Mesh(self,mesh_number)

        return mesh

    def make_movie(self,variables,framerate=10):
        """ Make a movie of variables """

        f = 1 #Frame counter

        #Loop over available meshes
        for m in range(1,self.Nmeshes+1):
            mesh = self.get_mesh(m)

            #Loop over time slices
            for t in range(0,mesh.Ntimes):
                mesh.make_plot(variables,t,f,purpose='movie')
                f += 1
        
        #Make video
        moviename = f'{self.directory}/movie/'
        for varname in variables:
                moviename += f'{varname}_'
        moviename = moviename[:-1] #Remove last underscore

        os.system(f'ffmpeg -r {framerate} -f image2 -i {self.directory}/movie/frame_%03d.png -pix_fmt yuv420p -vcodec libx264 -crf 24 {moviename}.mp4')

        #Remove frames
        os.system(f'rm {self.directory}/movie/frame*.png')

    # Helper functions

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
        self.contours = [
            var_name for var_name, var in ds.items() if (
                set(var.dims) in [{'time', 'two', 'ei'}, {'two', 'ei'}]
                and var_name[0] != 'E'
            )
        ]
        ds.close()

    def _get_times(self):
        """ Get time values per mesh """

        self.times = {}
        for f,fname in enumerate(self.fnames):
            ds = xr.open_dataset(fname)
            self.times[f+1] = ds.time.values
            ds.close()
            
