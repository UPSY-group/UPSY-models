import os
import glob
import xarray as xr
from matplotlib.collections import PolyCollection

from upsy.colormaps import *
from upsy.utils import *

class Mesh(object):
    """ Properties and functions of a single mesh """

    def __init__(self, Run, mesh_number):
        """ Gather basic info from run """

        self.Run = Run
        self.dir = self.Run.dir
        self.prefix = self.Run.prefix
        self.mesh_number = mesh_number

        self.got_voronois = False
        self.got_triangles = False

        self.open()
        self.close()

    def __repr__(self):
        return f"Mesh('{self.dir}',{self.mesh_number},'{self.prefix}')"

    def __str__(self):
        return f"Mesh number {self.mesh_number} of Run '{self.dir}'"

    def open(self):
        self.ds = xr.open_dataset(f'{self.dir}/{self.prefix}_{self.mesh_number:05d}.nc')
        self.Ntimes = len(self.ds.time)
    
    def close(self):
        self.ds.close()

    def get_voronois(self):
        """ Extract Voronoi cells as patches """
        self.voronois = []

        print(f"Computing {len(self.ds.vi)} new voronoi polygons ...")

        nVVor = self.ds.nVVor.values
        VVor = self.ds.VVor.values
        Vor = self.ds.Vor.values

        for vi in range(len(nVVor)):
            v_indices = VVor[:nVVor[vi], vi] -1
            x = Vor[0, v_indices]
            y = Vor[1, v_indices]
            self.voronois.append(np.column_stack((x,y)).tolist())

        self.got_voronois = True

        return f"Finished computing voronoi polygons"

    def get_triangles(self):
        """ Extract triangles as patches """
        self.triangles = []

        print(f"Computing {len(self.ds.ti)} new triangle polygons ...")

        Tri = self.ds.Tri.values
        V = self.ds.V.values

        for ti in range(0,len(self.ds.ti)):
            tri_indices = Tri[:, ti] -1
            x = V[0, tri_indices]
            y = V[1, tri_indices]
            self.triangles.append(np.column_stack((x,y)).tolist())
        
        self.got_triangles = True

        return f"Finished computing triangle polygons"

class Timeframe(object):
    """ Single timeframe of a mesh """

    def __init__(self, Mesh, t):
        """ Gather basic info from run """

        self.Mesh = Mesh
        self.t = t

        self.ds = self.Mesh.ds.isel(time=t)
        self.time = self.ds.time.values

        self.got_gl = False
        self.got_mask = False

    def __repr__(self):
        return f"Timeframe({repr(self.Mesh)},{self.t})"

    def __str__(self):
        return f"Timeframe {self.t} of Mesh number {self.mesh_number} of Run '{self.dir}'"

    def get_gl(self):
        """ Extract grounding line """

        # Read variable

        try:
            var = self.ds['grounding_line'].values
        except KeyError:
            print(f"ERROR: 'grounding_line' not in output files")
            return

        self.gl = var
        self.got_gl = True

        return

    def get_mask(self):
        """ Get a reduced mask separating grounded ice, ice shelf and ocean """

        mask = self.ds['mask']
        # Define as ocean ( = 2):

        # Define as sheet ( = 3):
        mask = xr.where(mask==1,3,mask) 
        mask = xr.where(mask==5,3,mask)
        mask = xr.where(mask==7,3,mask)
        mask = xr.where(mask==9,3,mask)
        mask = xr.where(mask==10,3,mask)
        # Define as shelf ( = 4):
        mask = xr.where(mask==6,4,mask)
        mask = xr.where(mask==8,4,mask)

        self.mask = mask.values
        self.got_mask = True

        return

    def get_pcoll(self,varname):
        """ Get patch collection """

        #Get data
        if varname == 'Uabs_lad':
            var1 = self.get_data('U_lad')
            var2 = self.get_data('V_lad')
            var = (var1**2+var2**2)**.5
        elif varname[:3] == 'BMB':
            var = self.get_data('BMB')
        else:
            var = self.get_data(varname)

        #Get colormap info
        scalarmap = get_cmap(varname)

        # Fill array
        if varname[:3] == 'BMB':
            #Reverse values for BMB
            pcols = scalarmap.to_rgba(-var.values)
        else:
            pcols = scalarmap.to_rgba(var.values)

        #Check type (voronoi / triangle)
        if 'vi' in var.dims:
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            pcoll = PolyCollection(self.Mesh.voronois, fc=pcols)
        elif 'ti' in var.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            pcoll = PolyCollection(self.Mesh.triangles, fc=pcols)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        return pcoll

    def get_data(self,varname):
        """ Get data array of variable """

        # Read variable
        try: 
            var = self.ds[varname]
        except:
            print(f'ERROR: could not read variable {varname}, make sure dataset is open and variable exists')
            return

        return var
