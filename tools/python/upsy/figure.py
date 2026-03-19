import os
import glob
import argparse

import xarray as xr
import matplotlib as mpl
from matplotlib.collections import PolyCollection
import matplotlib.gridspec as gridspec

from upsy.run import Run
from upsy.mesh import Mesh, Timeframe
from upsy.colormaps import *
from upsy.utils import *

class Figure(object):
    """ Figure of subplots of Ufe/Lad output """

    def __init__(self,
        figsize: tuple = (7,5),
        directory: str = 'figures',
        orientation: str = 'horizontal',
        cbar_ratio: int | float = 30,
        dpi: int = 1200,
        shrink: list = [0,0,0,0],
    ):

        """ Some input parameters with default settings """

        self.figsize = figsize
        self.directory = directory
        self.orientation = orientation
        self.cbar_ratio = cbar_ratio
        self.dpi = dpi
        self.shrink = shrink

        self.fields = {}

        assert self.orientation in ['horizontal','vertical'], f'Invalid orientation: {self.orientation}'

    def __repr__(self):
        return f"Figure()"

    def __str__(self):
        return f"Figure"

    def add_field(self, tf, varname, mask=None):
        """Add a field (say BMB of the last time frame) to the figure"""

        field = Field(tf, varname, mask=mask)

        print(f'Added {varname} with hash {hash(field)}')

        self.fields[hash(field)] = field
        return 

    def add_diff(self, Timeframe1, varname1, Timeframe2, varname2, name='diff', mask=None, vmax=None, colmap='cmo.diff'):
        """ Add a difference field to the figure, 
        for example the difference of 1 variable between 2 timeframes
        or the difference of two variables on the same timeframe
        """

        field1 = Field(Timeframe1, varname1)
        field2 = Field(Timeframe2, varname2)
        field = DiffField(field1, field2, name=name, mask=mask, vmax=vmax, colmap=colmap)

        print(f'Added difference field {name} with hash {hash(field)}')

        self.fields[hash(field)] = field
        return 

    def delete_field(self,fieldhash):
        # Delete a field from the figure

        del self.fields[fieldhash]

        return

    def print_fields(self):
        """Print all fields added to the figure thus far"""

        print(self.fields.keys())

        return

    def make(self,figname, add_gl=True, add_time=False):
        """Make the actual figure with all the added fields 
        It aligns panels either horizontally or vertically,
        along with a colorbar per panel.
        If requested by add_gl, the grouding line is added to each panel
        This routine will spit out a .png file in the requested directory
        """
        

        fig = plt.figure(figsize=self.figsize,constrained_layout=True)

        # Design the layout of the panels
        if self.orientation=='horizontal':
            # Align all panels horizontally, add a row below for the colorbar
            nrows = 2
            ncols = len(self.fields)
            spec = gridspec.GridSpec(nrows,ncols,figure=fig,height_ratios=[self.cbar_ratio,1],hspace=0.01,wspace=0.01)
            for k,key in enumerate(self.fields.keys()):
                self.fields[key].ax = fig.add_subplot(spec[0,k])
                self.fields[key].cax = fig.add_subplot(spec[1,k])
        else:
            # Align all panels vertically, add a row to the right for the colorbar
            nrows = len(self.fields)
            ncols = 2
            spec = gridspec.GridSpec(nrows,ncols,figure=fig,width_ratios=[self.cbar_ratio,1],hspace=0.01,wspace=0.01)
            for k,key in enumerate(self.fields.keys()):
                self.fields[key].ax = fig.add_subplot(spec[k,0])
                self.fields[key].cax = fig.add_subplot(spec[k,1])

        # Plot fields on subpanels
        for k,key in enumerate(self.fields.keys()):
            field = self.fields[key]
            if field.mask is not None:
                # Add mask of, say, grounded ice/ocean
                field.get_mcoll()
                field.ax.add_collection(field.mcoll)
            # Add field itself
            field.get_pcoll()
            field.ax.add_collection(field.pcoll)
            # Add colorbar
            cbar = plt.colorbar(field.scalarmap,cax=field.cax,orientation=self.orientation)
            cbar.set_label(field.varname)

            if add_gl:
                # Add the grounding line with some default settings
                if not field.Timeframe.got_gl:
                    field.Timeframe.get_gl()
                field.ax.plot(field.Timeframe.gl[0,:],field.Timeframe.gl[1,:],c='k',lw=.5)

            # Shrink each subpanel as requested
            field.ax.set_xlim([field.xmin+self.shrink[0], field.xmax-self.shrink[1]])
            field.ax.set_ylim([field.ymin+self.shrink[2], field.ymax-self.shrink[3]])
            field.ax.set_aspect(1)
            field.ax.set_xticks([])
            field.ax.set_yticks([])

        if add_time:
            # Add the timestamp, mostly useful for movies
            fig.suptitle(f'Year {field.Timeframe.time:.0f}')

        # Save figure in desired directory
        fullfigname = os.path.join(self.directory,f'{figname}.png')
        plt.savefig(fullfigname, bbox_inches = 'tight', pad_inches = 0,dpi=self.dpi)
        print(f'Created {figname}')
        plt.close()

        return

def make_2dplot(
    rundir: str, #Directory of the run where output is stored
    variables: list | str, #Variables to plot
    timeindices: list | int, #Time indices (integers) to plot
    figsize: str | tuple = (7,5), #Specify the figure size: (width,height)
    orientation: str = 'horizontal', #Alignment of panels
    cbar_ratio: int | float = 30, #Relative width/height of the colorbar
    directory: str = 'figures', #Directory to plot figure
    shrink: str | list = [0,0,0,0], #Shrink the domain size of each panel
    dpi: str | int = 450, #DPI of the output figure
):

    """ Function to make the actual plot,
    allowing user input """

    # Change input types to required ones
    try:
        figsize = tuple(figsize)
        shrink = list(shrink)
        dpi = int(dpi)
        variables = list(variables)
        timeindices = list(timeindices)
        cbar_ratio = float(cbar_ratio)
    except ValueError:
        raise argparse.ArgumentTypeError("Invalid input type")

    # Check whether input is valid: either 1 variable, 1 timeindex, or the same number of both
    if (len(variables) > 1 and len(timeindices) > 1):
        assert len(variables) == len(timeindices), f'Number of variables ({len(variables)}) must be equal to number of timeindices ({len(timeindices)})'
    elif len(variables) > 1:
        timeindices = timeindices*len(variables)
    elif len(timeindices) > 1:
        variables = variables*len(timeindices)

    # Make directory
    os.makedirs(directory,exist_ok=True)

    # Prepare figure
    fig = Figure(
        figsize=figsize,
        directory=directory,
        orientation=orientation,
        cbar_ratio=cbar_ratio,
        dpi=dpi,
        shrink=shrink,
    )

    # Open run
    run = Run(rundir)

    # Define beginning of figure name
    figname = f'{os.path.basename(rundir)}_'

    # Add fields for desired variables and timeindices
    for tidx,varname in zip(timeindices,variables):

        #TODO find mesh number. For now, only works on the first one
        midx = 1
        mesh = Mesh(run,midx)
        tf = Timeframe(mesh, tidx)
        fig.add_field(tf=tf, varname=varname, mask=None)
        figname += f'{varname}_{tidx}_'

    # Make the actual figure, leave out the last underscore
    fig.make(figname[:-1])


def main():
    """ This is the command-line interface to plot a 2D figure
    Example usage:

    uspy-plot-2dfigure rundir (default: show BMB on last timeframe)
    uspy-plot-2dfigure rundir -v R BMB_v3 (show two variables: resolution and BMB)
    uspy-plot-2dfigure rundir -t 0 -1 (show two timeframes of BMB)

    """

    # Parse command-line input
    parser = argparse.ArgumentParser(
    description='Make figure of subplots'
    )

    parser.add_argument(
        'rundir',
        help='Run directory where output is stored')

    parser.add_argument(
        '-D',
        '--directory',
        dest='directory',
        type = str,
        default='',
        help='Directory to save figure'
    )

    parser.add_argument(
        '-fs',
        '--figsize',
        dest='figsize',
        type = tuple,
        default=(7,5),
        help='Figure size (width, height)'
    )

    parser.add_argument(
        '-o',
        '--orientation',
        dest='orientation',
        type = str,
        default= 'horizontal',
        help="Orientation of panels: 'horizontal' or 'vertical'"
    )

    parser.add_argument(
        '-c',
        '--cbar_ratio',
        dest='cbar_ratio',
        type = float,
        default= 30,
        help="Ratio of panels height/width versus colorbar width"
    )

    parser.add_argument(
        '-v',
        '--variables',
        dest='variables',
        type = str,
        default = ['BMB_v3'],
        nargs = '+',
        help='Variable names'
    )

    parser.add_argument(
        '-t',
        '--timeindices',
        dest='timeindices',
        type = int,
        default = [-1],
        nargs = '+',
        help='Time indices'
    )

    parser.add_argument(
        '-s',
        '--shrink',
        dest='shrink',
        nargs=4,
        type=float,
        default=[1.0,1.0,1.0,1.0],
        help='Shrink domain [xmin, xmax, ymin, ymax] in km'
    )

    parser.add_argument(
        '-d',
        '--dpi',
        dest='dpi',
        type=int,
        default=1200,
        help='DPI of figure'
    )

    args = parser.parse_args()

    # Set default location to rundir/figures
    if args.directory == '':
        args.directory = os.path.join(args.rundir,'figures')

    # Make the plot
    make_2dplot(
        rundir = args.rundir,
        variables = args.variables,
        timeindices = args.timeindices,
        figsize = args.figsize,
        orientation = args.orientation,
        cbar_ratio = args.cbar_ratio,
        directory = args.directory,
        shrink = args.shrink,
        dpi = args.dpi,
    )


class Field(object):
    """ Field to include in subplot """

    def __init__(self, Timeframe, varname, mask=None):
        """ Input:
        Timeframe: class(Timeframe)
        varname: str
        mask: desired mask, default: None
        """

        # Extract info from input
        self.Timeframe = Timeframe
        self.Mesh = self.Timeframe.Mesh
        self.t = self.Timeframe.t
        self.varname = varname
        self.name = f"{varname}_{self.t}"
        self.mask = mask

        self.xmin = self.Timeframe.ds.xmin
        self.xmax = self.Timeframe.ds.xmax
        self.ymin = self.Timeframe.ds.ymin
        self.ymax = self.Timeframe.ds.ymax

        # Define the colormap based on the variable name
        # These are pre-defined but can be overwritten
        self.get_cmap()

        # Read the field data
        self.get_data()

    def __repr__(self):
        return f"Field({repr(self.Timeframe)}, {self.name})"

    def __str__(self):
        return f"Field {self.name} at time {self.Timeframe.t}"

    def get_data(self):
        """ Extract data values """

        # Special cases
        if self.varname == 'Uabs_lad':
            # Absolute speed, based on U and V velocity
            uvar = self.Timeframe.ds['U_lad']
            vvar = self.Timeframe.ds['V_lad']
            self.data = (uvar**2+vvar**2)**.5
        elif self.varname[:3] == 'BMB':
            # BMB or melt, with different versions
            if 'BMB' in self.Timeframe.ds:
                self.data = -self.Timeframe.ds['BMB']
            elif 'melt' in self.Timeframe.ds:
                self.data = self.Timeframe.ds['melt']*3600*24*365.25
            else:
                print(f"ERROR: no valid BMB or melt variable in Timeframe")
                return
        elif self.varname == 'C_friction':
            # Friction coefficient
            var = self.Timeframe.ds['till_friction_angle']
            self.data = np.tan(np.deg2rad(var)) 
        elif self.varname == 'duabs_surf':
            # Difference in uabs_surf compared to reference
            uref = self.Timeframe.ds['bed_roughness_nudge_H_u_target_velocity']
            uvar = self.Timeframe.ds['uabs_surf']
            self.data = xr.where(uref>0,uvar-uref,0)
        else:
            # Regular case: read variable if available in output
            try:
                self.data = self.Timeframe.ds[self.varname]
            except KeyError:
                print(f"ERROR: {self.varname} not in Timeframe")
                return

        if self.mask is not None:
            # Read mask, and mask out data by setting to NaN
            if 'vi' in self.data.dims:

                # Read mask values
                if not self.Timeframe.got_mask:
                    self.Timeframe.get_mask()

                # Mask out data in desired regions
                if self.mask == 'shelf':
                    self.data = xr.where([x in [4] for x in self.Timeframe.mask], self.data, np.nan)
                elif self.mask == 'sheet':
                    self.data = xr.where([x in [3] for x in self.Timeframe.mask], self.data, np.nan)
                elif self.mask == 'ice':
                    self.data = xr.where([x in [3,4] for x in self.Timeframe.mask], self.data, np.nan)
                else:
                    raise ValueError('Invalid option for mask')

            elif 'ti' in self.data.dims:
                # For data on triangles, no clean option yet
                # Just mask out region where data values are zero

                self.data = xr.where(self.data==0, np.nan, self.data)
            else:
                print(f'ERROR: variable {varname} is not on vertices or triangles')

            message = f"Extracted {self.varname} with mask {self.mask}"
        else:
            message = f"Extracted {self.varname}"

        return message

    def get_cmap(self):
        """ Get colormap info """
        scalarmap = get_cmap(self.varname)

        self.scalarmap = scalarmap

        return

    def get_pcoll(self):
        """ Get patch collection """

        pcols = self.scalarmap.to_rgba(self.data.values)

        #Check type (voronoi / triangle)
        if 'vi' in self.data.dims:
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            self.pcoll = PolyCollection(self.Mesh.voronois, fc=pcols, ec=pcols, lw=.1)
        elif 'ti' in self.data.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            self.pcoll = PolyCollection(self.Mesh.triangles, fc=pcols, ec=pcols, lw=.1)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        return

    def get_mcoll(self):
        """ Get mask collection """

        mcmap = 'ocean'
        mnorm = mpl.colors.Normalize(vmin=1.7, vmax=3.1, clip=True)
        scalarmap = mpl.cm.ScalarMappable(norm=mnorm,cmap=mcmap)

        if not self.Mesh.got_voronois:
            self.Mesh.get_voronois()

        mcols = scalarmap.to_rgba(self.Timeframe.mask)
        self.mcoll = PolyCollection(self.Mesh.voronois,fc=mcols, ec=mcols, lw=.1)

        return

class DiffField(object):
    """ Difference field between two Fields 

        Diff == Field1 - Field2
    """

    def __init__(self, Field1, Field2, name, mask, vmax, colmap):

        # Read two fields
        self.Field1 = Field1
        self.Field2 = Field2
        # Read the timeframe of the first field, for example to plot the grounding line
        self.Timeframe = self.Field1.Timeframe

        self.name = name
        self.varname = 'diff'
        self.mask = mask

        self.xmin = self.Field1.xmin
        self.xmax = self.Field1.xmax
        self.ymin = self.Field1.ymin
        self.ymax = self.Field1.ymax

        # Define the colormap and range (assumed to be symmetric: -vmax to vmax)
        self.vmax = vmax
        self.colmap = colmap

        # Check whether the two fields are compatible (on the same mesh, both voronoi or triangle)
        self.check_compatibility()

        self.get_cmap()

        self.get_data()

    def __repr__(self):
        return f"DiffField({repr(self.Field1)}, {repr(self.Field2)}, {self.name})"

    def __str__(self):
        return f"DiffField between {self.Field1.name} at time {self.Field1.Timeframe.t} and {self.Field1.name} at time {self.Field1.Timeframe.t}"

    def check_compatibility(self):
        """ Check whether difference field can be computed """

        # Check whether the domain is equal
        assert self.xmin == self.Field2.xmin, 'ERROR: xmin is not equal between input Fields'
        assert self.xmax == self.Field2.xmax, 'ERROR: xmax is not equal between input Fields'
        assert self.ymin == self.Field2.ymin, 'ERROR: ymin is not equal between input Fields'
        assert self.ymax == self.Field2.ymax, 'ERROR: ymax is not equal between input Fields'

        # Check whether the dimensions (vi or ti) are equal
        if 'vi' in self.Field1.data.dims:
            assert 'vi' in self.Field2.data.dims, 'ERROR: cannot take difference between voronois and triangles'
            assert (self.Field1.Timeframe.ds['V'] == self.Field2.Timeframe.ds['V']).all(), 'ERROR: V is not equal between input Fields'
        elif 'ti' in self.Field1.data.dims:
            assert 'ti' in self.Field2.data.dims, 'ERROR: cannot take difference between voronois and triangles'
            assert (self.Field1.Timeframe.ds['Tri'] == self.Field2.Timeframe.ds['Tri']).all(), 'ERROR: ti is not equal between input Fields'

        print('Input Fields are compatible')

        # Inherit the mesh
        self.Mesh = self.Field1.Timeframe.Mesh

        return

    def get_cmap(self):
        """ Define cmap and norm """

        self.cmap = plt.get_cmap(self.colmap)

        if self.vmax is None:
            self.dmax = max(abs(self.Field1.data-self.Field2.data).values)
            if self.dmax == 0:
                print('WARNING: both fields are equal')
                vmax = 1
            else:
                vmax = self.dmax
        else:
            vmax = self.vmax

        self.norm = mpl.colors.Normalize(vmin=-vmax,vmax=vmax,clip=True)
        return

    def get_data(self):
        """ Get data and add mask if necessary """

        self.data = self.Field1.data - self.Field2.data

        if self.mask is not None:
            if 'vi' in self.data.dims:
                if not self.Timeframe.got_mask:
                    self.Timeframe.get_mask()

                if self.mask == 'shelf':
                    self.data = xr.where([ x in [4] for x in self.Timeframe.mask], self.data, np.nan)
                elif self.mask == 'sheet':
                    self.data = xr.where([ x in [3] for x in self.Timeframe.mask], self.data, np.nan)
                elif self.mask == 'ice':
                    self.data = xr.where([x in [3,4] for x in self.Timeframe.mask], self.data, np.nan)
                else:
                    raise ValueError('Invalid option for mask')
            elif 'ti' in self.data.dims:
                # Just mask out zero values, no clean option yet
                self.data = xr.where(self.data==0, np.nan, self.data)
            else:
                print(f'ERROR: something went wrong with vertices and triangles')

            message = f"Extracted {self.name} with mask {self.mask}"
        else:
            message = f"Extracted {self.name}"

        return message

    def get_pcoll(self):
        """ Get patch collection """

        # Fill array
        pcols = self.scalarmap.to_rgba(self.data.values)

        #Check type (voronoi / triangle)
        if 'vi' in self.data.dims:
            if not self.Mesh.got_voronois:
                self.Mesh.get_voronois()
            self.pcoll = PolyCollection(self.Mesh.voronois, fc=pcols)
        elif 'ti' in self.data.dims:
            if not self.Mesh.got_triangles:
                self.Mesh.get_triangles()
            self.pcoll = PatchCollection(self.Mesh.triangles, fc=pcols)
        else:
            print(f'ERROR: variable {varname} is not on vertices or triangles')

        return

    def get_mcoll(self):
        """ Get mask collection """

        mcmap = 'ocean'
        mnorm = mpl.colors.Normalize(vmin=1.7, vmax=3.1, clip=True)

        scalarmap = mpl.cm.ScalarMappable(norm=mnorm,cmap=mcmap)
        mcols = scalarmap.to_rgba(self.Timeframe.mask)

        if not self.Mesh.got_voronois:
            self.Mesh.get_voronois()

        self.mcoll = PolyCollection(self.Mesh.voronois, fc=mcols)

        return
