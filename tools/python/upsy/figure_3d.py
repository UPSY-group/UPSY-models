import os
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from tqdm import tqdm

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from upsy.main_run import *
from upsy.main_mesh import *

def make_3dplot(
    rundir: str, 
    azilight: str | int | float = 0, 
    elelight: str | int | float = 45, 
    aziview: str | int | float = 190, 
    eleview: str | int | float = 15, 
    vmax: str | int | float = 200
):

    try:
        azilight = float(azilight)
        elelight = float(elelight)
        aziview = float(aziview)
        eleview = float(eleview)
        vmax = float(vmax)
    except ValueError:
        raise argparse.ArgumentTypeError("Must be a floating point number")

    scalarmap = _get_cmap(vmax=vmax)

    #Extract timeframe
    run = Run(rundir)
    run.Nmeshes = 1
    mesh = run.get_mesh(1,file='laddie_output')
    tf = Timeframe(mesh,-1)

    dirname = os.path.join(run.directory,'figures')
    os.makedirs(dirname,exist_ok=True)

    print('Getting verts')
    #Get surface faces
    verts_oce = _verts_vi(tf, 0*tf.ds.Hs_b)
    verts_Hs = _verts_vi(tf, tf.ds.Hs_b)

    print('Getting curts')
    #Get curtains
    curts_cf_fl, curts_vi_cf_fl = _curts_vi(tf, tf.ds.Hs_b, 0*tf.ds.Hs_b, [4], [2,8])
    curts_cf_gr, curts_vi_cf_gr = _curts_vi(tf, tf.ds.Hs_b, 0*tf.ds.Hs_b, [3,5,6,7,9,10], [2,8])

    print('Getting hillshades')
    #Get hillshades
    hn_oce = hillshade(0*tf.ds.dHs_dx,0*tf.ds.dHs_dy,azimuth=azilight,altitude=elelight)
    hn_Hs = hillshade(tf.ds.dHs_dx,tf.ds.dHs_dy,azimuth=azilight,altitude=elelight)

    print('Adding colors to verts')
    #Add colors to verts
    cols_oce = [(0,0,0,0)]*len(verts_oce)
    cols_ice = [(0,0,0,0)]*len(verts_Hs)
    cols_bmb = [(0,0,0,0)]*len(verts_Hs)
    cols_bed = [(0,0,0,0)]*len(verts_Hs)

    for vi in tqdm(range(0,len(tf.ds.vi))):
        #BMB
        if (tf.ds.mask[vi] in [4]):
            rgba_bg = scalarmap.to_rgba(tf.ds.melt[vi]*3600*24*365.25)
            cols_bmb[vi] = rgba(rgba_bg,hn_Hs[vi])
    
        #Ocean
        elif (tf.ds.mask[vi] in [2,8]):
            rgba_bg = mpl.colors.to_rgba('darkslategray')
            cols_oce[vi] = rgba(rgba_bg,hn_oce[vi])
    
        # Ice
        elif (tf.ds.mask[vi] in [3,5,6,7,9,10]):
            rgba_bg = mpl.colors.to_rgba('lightblue')
            cols_ice[vi] = rgba(rgba_bg,hn_Hs[vi],alpha=.4)
        # Bed
        elif (tf.ds.mask[vi] in [1]):
            rgba_bg = mpl.colors.to_rgba('saddlebrown')
            cols_bed[vi] = rgba(rgba_bg,hn_Hs[vi],alpha=.4)

    print('Adding colors to curts')
    #Add colors to curts
    cols_curts_cf_fl = [(0,0,0,0)]*len(curts_cf_fl)
    for i,curt in tqdm(enumerate(curts_cf_fl)):
        vi = curts_vi_cf_fl[i]
        rgba_bg = scalarmap.to_rgba(tf.ds.melt[vi]*3600*24*365.25)
        dy = curt[2][1] - curt[1][1]
        dx = curt[2][0] - curt[1][0]
        hn = curtshade(dy,dx,azilight,elelight)
        cols_curts_cf_fl[i] = rgba(rgba_bg,hn)

    cols_curts_cf_gr = [(0,0,0,0)]*len(curts_cf_gr)
    for i,curt in tqdm(enumerate(curts_cf_gr)):
        vi = curts_vi_cf_gr[i]
        rgba_bg = mpl.colors.to_rgba('lightblue')
        dy = curt[2][1] - curt[1][1]
        dx = curt[2][0] - curt[1][0]
        hn = curtshade(dy,dx,azilight,elelight)
        cols_curts_cf_gr[i] = rgba(rgba_bg,hn,alpha=.4)

    print('Making figure')
    #Make figure
    fig = plt.figure(figsize=(7,3))#, constrained_layout=True)
    ax = fig.add_subplot(projection='3d')
    
    allv = verts_Hs*3 + verts_oce + curts_cf_fl + curts_cf_gr
    allc = cols_ice + cols_bmb + cols_bed + cols_oce + cols_curts_cf_fl + cols_curts_cf_gr
    
    poly = Poly3DCollection(allv,fc=allc,axlim_clip=True)
    ax.add_collection3d(poly)
    
    ax.set_aspect('equalxy')
    ax.set_box_aspect((1,1,.1))
    ax.set_ylim([tf.ds.ymin+2e3,tf.ds.ymax-2e3])
    ax.set_xlim([tf.ds.xmin+2e3,tf.ds.xmax-2e3])
    ax.set_zlim(zmin=-10,zmax=2000)
    
    ax.view_init(elev=eleview, azim=aziview)
    ax.set_axis_off()
    
    fig.subplots_adjust(left=-.8,right=1.8,top=1.8,bottom=-.8)

    filename = os.path.join(dirname,'3Dplot.png')
    plt.savefig(filename,dpi=1200)#, bbox_inches='tight')

    print(f'Finished {filename}')

def main():
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

    args = parser.parse_args()

    make_3dplot(
        rundir = args.rundir,
        azilight = args.azilight,
        elelight = args.elelight,
        aziview = args.aziview,
        eleview = args.eleview,
        vmax = args.vmax
    )

def _verts_vi(tf, H_b):
    verts = []
    for vi in tqdm(range(0,len(tf.ds.vi))):
        #Get number of voronoi vertices
        nVVor = tf.ds.nVVor[vi].values
        #Get indices of voronoi vertices
        VVor = tf.ds.VVor[:nVVor,vi].values
        #Get 3D locations
        x = tf.ds.Vor[0,VVor-1].values
        y = tf.ds.Vor[1,VVor-1].values
        # Get triangle indices
        ti = tf.ds.vori2ti[VVor-1].values
        #Add Hs
        z = H_b[ti-1].values
        v = np.array([x,y,z]).T.tolist()
        verts.append(v)
    
    return verts

def _curts_vi(tf, H_b0, H_b1, maskvals, neighbs):

    curts = []
    curts_vi = []
    for ei in tqdm(range(0,len(tf.ds.ei))):
        v0 = tf.ds.EV[0,ei].values-1
        v1 = tf.ds.EV[1,ei].values-1

        
        if (tf.ds.mask[v0].values in maskvals and tf.ds.mask[v1].values in neighbs):
            vi = v0
            t0 = tf.ds.ETri[1,ei].values-1
            t1 = tf.ds.ETri[0,ei].values-1
        elif (tf.ds.mask[v1].values in maskvals and tf.ds.mask[v0].values in neighbs):
            vi = v1
            t0 = tf.ds.ETri[0,ei].values-1
            t1 = tf.ds.ETri[1,ei].values-1
        else:
            continue

        x0 = tf.ds.Tricc[0,t0].values
        x1 = tf.ds.Tricc[0,t1].values
        y0 = tf.ds.Tricc[1,t0].values
        y1 = tf.ds.Tricc[1,t1].values

        z00 = H_b0[t0].values
        z01 = H_b1[t0].values
        z10 = H_b0[t1].values
        z11 = H_b1[t1].values

        xx = [x0, x0, x1, x1]
        yy = [y0, y0, y1, y1]
        zz = [z00, z01, z11, z10]

        c = np.array([xx,yy,zz]).T.tolist()
        curts.append(c)
        curts_vi.append(vi)

    return curts, curts_vi

def hillshade(dHdx,dHdy,azimuth,altitude,z_fac=1000):
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
    """ BMB colormap """

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
    rgb = [hn * alpha + rgba_bg[i] * (1-alpha) for i in [0,1,2]]
    rgba = tuple(np.append(rgb,1))
    return rgba
