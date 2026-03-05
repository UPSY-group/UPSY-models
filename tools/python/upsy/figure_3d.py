import os
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from tqdm import tqdm

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from upsy.run import Run
from upsy.mesh import Mesh, Timeframe

def make_3dplot(
    rundir: str, # Directory of the run where output is stored
    azilight: str | int | float = 0, # Azimuth of the light source for hillshading
    elelight: str | int | float = 45, # Elevation of the light source for hillshading
    aziview: str | int | float = 190, # Azimuth of the viewpoint
    eleview: str | int | float = 15, # Elevation of the viewpoint
    vmax: str | int | float = 200, # Maximum of the colorbar
    shrink: str | list = [1,1,1,1], # Number of gridcells to cut from the domain to avoid edge effects
    dpi: str | int = 1200, # DPI of the output figure
    vaspect: str | float = 0.1, # Vertical aspect ratio to change the exaggeration of the vertical dimension
    linewidth: str | int = 0.0, # Linewidth around the cells
):

    """
    Make a fancy 3D plot of your output.
    This requires the output variables Hs, Hs_b, dHs_dx, and dHs_dy
    Also, it's now hardcoded only for BMB (UFEMISM run) or melt (LADDIE), 
    which should be available in the output
    """

    # Change input types to required ones
    try:
        azilight = float(azilight)
        elelight = float(elelight)
        aziview = float(aziview)
        eleview = float(eleview)
        vmax = float(vmax)
        shrink = list(shrink)
        dpi = int(dpi)
        vaspect = float(vaspect)
        linewidth = float(linewidth)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")

    scalarmap = _get_cmap(vmax=vmax)

    #Extract timeframe
    run = Run(rundir)
    mesh = Mesh(run,run.Nmeshes)
    tf = Timeframe(mesh,-1)

    dirname = os.path.join(run.dir,'figures')
    os.makedirs(dirname,exist_ok=True)

    print('Getting verts')
    #Get surface faces: the sloped voronoi cells
    verts_oce = _verts_vi(tf, 0*tf.ds.Hs_b)
    verts_Hs = _verts_vi(tf, tf.ds.Hs_b)

    print('Getting curts')
    #Get curtains: the vertical faces around the voronoi cells
    curts_cf_fl, curts_vi_cf_fl = _curts_vi(tf, tf.ds.Hs_b, 0*tf.ds.Hs_b, [4,6], [2,8])
    curts_cf_gr, curts_vi_cf_gr = _curts_vi(tf, tf.ds.Hs_b, 0*tf.ds.Hs_b, [3,5,7,9,10], [2,8])

    print('Getting hillshades')
    #Get hillshades
    hn_oce = hillshade(0*tf.ds.dHs_dx,0*tf.ds.dHs_dy,azimuth=azilight,altitude=elelight)
    hn_Hs = hillshade(tf.ds.dHs_dx,tf.ds.dHs_dy,azimuth=azilight,altitude=elelight)

    print('Adding colors to verts')
    #Add colors to verts (= surface faces)

    # Allocate colors
    cols_oce = [(0,0,0,0)]*len(verts_oce)
    cols_ice = [(0,0,0,0)]*len(verts_Hs)
    cols_bmb = [(0,0,0,0)]*len(verts_Hs)
    cols_bed = [(0,0,0,0)]*len(verts_Hs)

    # Determine mask
    mask = tf.ds.mask.values

    # Read the melt rates
    try:
        melt_vals = tf.ds.melt.values * 3600*24*365.25
    except:
        melt_vals = -tf.ds.BMB.values

    # This seems to be necessary...
    hn_Hs = np.array(hn_Hs)
    hn_oce = np.array(hn_oce)
    
    # Change matplotlib colors to rgba, and use predefined colors for masked surfaces
    rgba_bg_bmb = scalarmap.to_rgba(melt_vals)
    rgba_bg_oce = mpl.colors.to_rgba('darkslategray')
    rgba_bg_ice = mpl.colors.to_rgba('lightblue')
    rgba_bg_bed = mpl.colors.to_rgba('saddlebrown')

    # Determine actual colors based on colormap (or predefined color) and hillshade
    for vi in tqdm(range(len(mask))):
        #BMB
        if mask[vi] in [4,6]:
            cols_bmb[vi] = rgba(rgba_bg_bmb[vi],hn_Hs[vi])
    
        #Ocean
        elif mask[vi] in [2,8]:
            cols_oce[vi] = rgba(rgba_bg_oce,hn_oce[vi])
    
        # Ice
        elif mask[vi] in [3,5,7,9,10]:
            cols_ice[vi] = rgba(rgba_bg_ice,hn_Hs[vi],alpha=.4)
        # Bed
        elif mask[vi] in [1]:
            cols_bed[vi] = rgba(rgba_bg_bed,hn_Hs[vi],alpha=.4)

    print('Adding colors to curts')

    #Add colors to curts, same procedure as verts

    #Curtains (vertical faces) at the calving front of floating cells 
    cols_curts_cf_fl = [(0,0,0,0)]*len(curts_cf_fl)
    for i,curt in tqdm(enumerate(curts_cf_fl)):
        vi = curts_vi_cf_fl[i]
        dy = curt[2][1] - curt[1][1]
        dx = curt[2][0] - curt[1][0]
        hn = curtshade(dy,dx,azilight,elelight)
        cols_curts_cf_fl[i] = rgba(rgba_bg_bmb[vi],hn)

    # Curtains (vertical faces) at the calving front of grounded cells (= cliff)
    cols_curts_cf_gr = [(0,0,0,0)]*len(curts_cf_gr)
    for i,curt in tqdm(enumerate(curts_cf_gr)):
        vi = curts_vi_cf_gr[i]
        dy = curt[2][1] - curt[1][1]
        dx = curt[2][0] - curt[1][0]
        hn = curtshade(dy,dx,azilight,elelight)
        cols_curts_cf_gr[i] = rgba(rgba_bg_ice,hn,alpha=.4)

    print('Making figure')

    #Make figure

    #Prepare figure
    fig = plt.figure(figsize=(7,3))#, constrained_layout=True)
    ax = fig.add_subplot(projection='3d')
    
    # Add faces (verts, curts) and associated colors together
    allv = verts_Hs*3 + verts_oce + curts_cf_fl + curts_cf_gr
    allc = cols_ice + cols_bmb + cols_bed + cols_oce + cols_curts_cf_fl + cols_curts_cf_gr
    
    # Create a polygon of the thing and add to the axis
    poly = Poly3DCollection(allv,fc=allc,lw=linewidth,edgecolor='k',axlim_clip=True)
    ax.add_collection3d(poly)
    
    # Set the boundary domain, aspect, and dimension limits
    ax.set_aspect('equalxy')
    ax.set_box_aspect((1,1,vaspect))
    ax.set_xlim([tf.ds.xmin+shrink[0]*1e3, tf.ds.xmax-shrink[1]*1e3])
    ax.set_ylim([tf.ds.ymin+shrink[2]*1e3, tf.ds.ymax-shrink[3]*1e3])
    ax.set_zlim(zmin=min(tf.ds.Hs.values)-10,zmax=max(tf.ds.Hs.values)+10)
    
    # Set the viewpoint
    ax.view_init(elev=eleview, azim=aziview)
    ax.set_axis_off()

    # This is a nasty hardcoded adjustment to make sure the figure is fully visible
    # without too large white spaces around it. There doesn't seem to be a clean
    # way to do this, so sticking to it for now.
    fig.subplots_adjust(left=-.8,right=1.8,top=1.8,bottom=-.8)

    filename = os.path.join(dirname,'3Dplot.png')

    # Add metadata of choices to the figure for reproducibility and optimisation, 
    # it's not so easy to find though...
    metadata = {
        'azilight': str(azilight),
        'elelight': str(elelight),
        'aziview': str(aziview),
        'eleview': str(eleview),
        'vmax': str(vmax),
        'shrink': str(shrink),
        'dpi': str(dpi),
        'vaspect': str(vaspect),
        'linewidth': str(linewidth)
    }

    # Just print the metadata to screen so adjustments (viewpoint, vertical aspect, etc) can be optimised
    print(metadata)

    # Save the thing
    plt.savefig(filename,dpi=dpi, metadata=metadata)#, bbox_inches='tight')

    print(f'Finished {filename}')

def main():
    """
    Command-line interface to plot the 3D figure

    Example usage:

    upsy-plot-3dfigure rundir
    """

    # Parse command-line input
    parser = argparse.ArgumentParser(
    description='Make 3D plot'
    )

    parser.add_argument(
        'rundir',
        help='Run directory where output is stored')

    parser.add_argument(
        '-al',
        '--azilight',
        dest='azilight',
        default=0,
        help='Azimuth of light source for hillshade, 0 = south'
    )

    parser.add_argument(
        '-el',
        '--elelight',
        dest='elelight',
        default=45,
        help='Elevation of light source for hillshade, between 0 and 90'
    )

    parser.add_argument(
        '-av',
        '--aziview',
        dest='aziview',
        default=190,
        help='Azimuth of viewpoint, 0 = east'
    )

    parser.add_argument(
        '-ev',
        '--eleview',
        dest='eleview',
        default=15,
        help='Elevation of view point, between 0 and 90'
    )

    parser.add_argument(
        '-vm',
        '--vmax',
        dest='vmax',
        default=200,
        help='Maximum value of BMB colormap'
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

    parser.add_argument(
        '-va',
        '--vaspect',
        dest='vaspect',
        type=float,
        default=0.1,
        help='Vertical aspect ratio /1000. Value 0.1 equals a 100x vertical exaggeration'
    )

    parser.add_argument(
        '-lw',
        '--linewidth',
        dest='linewidth',
        type=float,
        default=0.0,
        help='Linewidth of mesh'
    )

    args = parser.parse_args()

    # Make the plot
    make_3dplot(
        rundir = args.rundir,
        azilight = args.azilight,
        elelight = args.elelight,
        aziview = args.aziview,
        eleview = args.eleview,
        vmax = args.vmax,
        shrink = args.shrink,
        dpi = args.dpi,
        vaspect = args.vaspect,
        linewidth = args.linewidth,
    )

# Helper functions

def _verts_vi(tf, H_b):
    """
    Determine the corner points (horizontal and vertical of the verts (surface faces).
    These are essentially the voronoi cells with their vertical dimension.
    Input:
    tf: Timeframe class
    H_b: Height variable on b-grid, typically Hs_b
    """

    verts = []

    # Extract necessary data from the dataset at once
    nVVor = tf.ds.nVVor.values
    VVor = tf.ds.VVor.values
    Vor = tf.ds.Vor.values
    vori2ti = tf.ds.vori2ti.values

    for vi in tqdm(range(0,len(nVVor))):
        # Get number of Voronoi vertices for current index
        num_vvor = nVVor[vi]
        
        # Get indices of Voronoi vertices
        vvor_indices = VVor[:num_vvor, vi] - 1

        # Get 3D locations
        x = Vor[0, vvor_indices]
        y = Vor[1, vvor_indices]
        
        # Get triangle indices and corresponding H
        ti = vori2ti[vvor_indices] -1
        z = H_b[ti]
        
        # Combine x, y, z into a vertex array
        v = np.column_stack((x, y, z)).tolist()
        verts.append(v)
    
    return verts

def _curts_vi(tf, H_b0, H_b1, maskvals, neighbs):
    """
    Determine the cornerpoints of the curtains (vertically oriented rectangles)
    Input:
    tf: Timeframe class
    H_b0: Height variable on b-grid of the bottom two cornerpoints
    H_b1: Height variable on b-grid of the top two cornerpoints
    maskvals: List of mask values, 
        only form a curtain for cells containing a mask value equal to one of these values
    neighbs: List of mask values for the opposite cell,
        only form a curtain if the opposite neighbour has one of these mask values

    Example for the calving front of floating ice shelves: 
        H_b0 would typically be Hs_b (ice surface)
        H_b1 would typically be 0 (sea level)
        maskvals should include 4 (floating ice) and often 6 (cf_fl (?) )
        neighbs should include 2 (ocean) and often 8 (cf_oc (?) )
    """

    curts = []
    curts_vi = []

    # Extract necessary data from the dataset
    EV = tf.ds.EV.values
    ETri = tf.ds.ETri.values
    Tricc = tf.ds.Tricc.values
    mask = tf.ds.mask.values

    for ei in tqdm(range(0,len(tf.ds.ei))):
        # Determine vertex indices at both sides of the face
        v0 = EV[0,ei] - 1
        v1 = EV[1,ei] - 1

        # Determine mask values of the two cells bordering the curtain
        mask_v0 = mask[v0]
        mask_v1 = mask[v1]

        # Check whether a curtain should be determined
        # If so, extract the triangle indices (borders of the curtain)
        if mask_v0 in maskvals and mask_v1 in neighbs:
            vi = v0
            t0 = ETri[1, ei] -1
            t1 = ETri[0, ei] -1
        elif mask_v1 in maskvals and mask_v0 in neighbs:
            vi = v1
            t0 = ETri[0, ei] -1
            t1 = ETri[1, ei] -1
        else:
            continue

        # Translate triangle indices to values on the horizontal plane
        x0, x1 = Tricc[0, t0], Tricc[0, t1]
        y0, y1 = Tricc[1, t0], Tricc[1, t1]

        # Add the vertical dimension to the four corner points
        z00, z01 = H_b0[t0], H_b1[t0]
        z10, z11 = H_b0[t1], H_b1[t1]

        # Combine horizontal and vertical dimensions
        c = np.array([[x0, x0, x1, x1], [y0, y0, y1, y1], [z00, z01, z11, z10]]).T.tolist()

        # Append the curtain and the vertex index to determine the color for this curtain
        curts.append(c)
        curts_vi.append(vi)

    return curts, curts_vi

def hillshade(dHdx,dHdy,azimuth,altitude,z_fac=1000):
    """
    Determine the hillshade factor to darken cells
    based on the angle of the surface (dHdx, dHdy)
    and the azimuth and elevation of the light source
    plus some vertical exaggeration factor
    """

    # Got this whole algorithm from the website of ArcGIS Pro:
    # https://pro.arcgis.com/en/pro-app/latest/tool-reference/3d-analyst/how-hillshade-works.htm
    # It seems to work well

    zenith_deg = 90-altitude
    zenith_rad = zenith_deg * np.pi / 180.0
    azimuth_math = azimuth+90
    if azimuth_math >= 360:
        azimuth_math += -360
    azimuth_rad = azimuth_math * np.pi / 180.0
    dzdx = dHdx
    dzdy = dHdy
    slope_angle = np.pi / 2.0 - np.arctan(z_fac * np.sqrt(dzdx**2+dzdy**2))
    aspect_angle = np.arctan2(-dzdx,dzdy)
    aspect_angle[aspect_angle<0] += 2*np.pi

    hillshade = (np.sin(zenith_rad) * np.sin(slope_angle) +
        np.cos(zenith_rad) * np.cos(slope_angle) *
        np.cos(azimuth_rad - np.pi / 2.0 - aspect_angle))

    hn = (np.maximum(0,hillshade)+1)/2

    return hn

def curtshade(dy,dx,azimuth,altitude,z_fac=1000):
    """
    Derived a similar algorithm for the shading of the curtains
    """

    zenith_deg = 90-altitude
    zenith_rad = zenith_deg * np.pi / 180.0
    azimuth_math = azimuth+90
    if azimuth_math >= 360:
        azimuth_math += -360
    azimuth_rad = azimuth_math * np.pi / 180.0
    aspect_angle = np.arctan2(dx,-dy)
    slope_angle = 0

    hillshade = (np.sin(zenith_rad) * np.sin(slope_angle) +
        np.cos(zenith_rad) * np.cos(slope_angle) *
        np.cos(azimuth_rad - aspect_angle))
    
    hn = (np.maximum(0,hillshade)+1)/2

    return hn

def _get_cmap(vmax=200):
    """ BMB colormap, hardcoded for now """

    Ncols = 68
    vmax = vmax 
    vmin = -10 
    linthresh = .3
    linscale = .25 
    fracpos = (np.log10(vmax/linthresh)+linscale)/(np.log10(vmax/linthresh)+np.log10(-(vmin/linthresh))+2*linscale)
    nneg = np.int_((1-fracpos)*Ncols) + 1 
    colors1 = plt.get_cmap('cmo.diff')(np.linspace(.15,.5,nneg))
    colors2 = plt.get_cmap('afmhot_r')(np.linspace(0,.9, Ncols-nneg))
    colors = np.vstack((colors1, colors2))
    
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors, Ncols)
    norm = mpl.colors.SymLogNorm(linthresh, vmin=vmin, vmax=vmax, linscale=linscale)
    scalarmap = mpl.cm.ScalarMappable(norm=norm,cmap=cmap)

    return scalarmap

def rgba(rgba_bg,hn,alpha=.2):
    """
    Combine RGB values with the normalised hillshading using
    an alpha factor to determine how 'dark' the shadow should be 
    """
    rgb = [hn * alpha + rgba_bg[i] * (1-alpha) for i in [0,1,2]]
    rgba = tuple(np.append(rgb,1))
    return rgba
