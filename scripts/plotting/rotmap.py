from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import os,imageio
import sys
import argparse

import fitsio
from astropy.table import join,Table

parser = argparse.ArgumentParser()
parser.add_argument("--types", help="tracer types to be selected",default='all')
parser.add_argument("--comp", help="comp tipe",default='COMP_TILE')
parser.add_argument("--nframe", help="number of images for gif",default=20)
args = parser.parse_args()


nframe = int(args.nframe)

outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/plots/'
#dat = fitsio.read('/Users/ashleyross/Dropbox/DESI/main/'+sys.argv[1]+'zdone_full_noveto.dat.fits',colums=['RA','DEC',sys.argv[2]])
#dat = fitsio.read('/Users/ashleyross/Dropbox/DESI/main/'+sys.argv[1]+'targets_withcomp.fits')

def gif_maker(gif_name,png_dir,gif_indx,num_gifs,dpi=90):
    #copied from
    # make png path if it doesn't exist already
    if not os.path.exists(png_dir):
        os.makedirs(png_dir)

    # save each .png for GIF
    # lower dpi gives a smaller, grainier GIF; higher dpi gives larger, clearer GIF
    plt.savefig(png_dir+'frame_'+str(gif_indx)+'_.png',dpi=dpi)
    plt.close('all') # comment this out if you're only updating the x,y data

    if gif_indx==num_gifs-1:
        # sort the .png files based on index used above
        images,image_file_names = [],[]
        for file_name in os.listdir(png_dir):
            if file_name.endswith('.png'):
                image_file_names.append(file_name)       
        sorted_files = sorted(image_file_names, key=lambda y: int(y.split('_')[1]))

        # define some GIF parameters
    
        frame_length = 0.5 # seconds between frames
        end_pause = 4 # seconds to stay on last frame
        # loop through files, join them to image array, and write to GIF called 'wind_turbine_dist.gif'
        for ii in range(0,len(sorted_files)):       
            file_path = os.path.join(png_dir, sorted_files[ii])
            if ii==len(sorted_files)-1:
                for jj in range(0,int(end_pause/frame_length)):
                    images.append(imageio.imread(file_path))
            else:
                images.append(imageio.imread(file_path))
        # the duration is the time spent on each image (1/duration is frame rate)
        imageio.mimsave(gif_name, images,'GIF',duration=frame_length)



if args.types == 'all':
    tpl = ['LRG','ELG','ELG_LOP','BGS_ANY','BGS_BRIGHT','QSO']
else:
    tpl = args.types.split(',')
for tp in tpl:
    tars = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tp+'targetsDR9v1.1.1.fits',columns=['TARGETID','RA','DEC'])
    dat = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full_noveto.dat.fits',columns=['TARGETID',args.comp])
    dat = join(tars,dat,keys=['TARGETID'],join_type='left')
    mrows = dat['COMP_TILE'].mask
    dat['OBS'] = np.zeros(len(dat))
    dat['OBS'][~mrows] = 1
    dat['COMP_TILE'][mrows] = 0



    sel = dat['OBS'] == 1
    dato = dat[sel]
    datno = dat[~sel]
    if tp == 'ELG':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.1
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]
    if tp == 'BGS_ANY':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.2
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]
    if tp == 'LRG':
         rns = np.random.random_sample(len(datno))
         sel = rns < 0.5
         #sel |= dat['COMP_TILE'] > 0
         datno = datno[sel]



    print(len(dato),len(datno))

    qt = args.comp
    vm = np.min(dat[qt])
    vx = np.max(dat[qt])
    #if qt == 'COMP_TILE' and sys.argv[1] == 'LRG':
    #    vm = 0.3
    #    vx = 1.

    # plt.ion() allows python to update its figures in real-time
    plt.ion()
    #fig = plt.figure(figsize=(9,6))

    # set perspective angle
    lat_viewing_angle = 30#50
    lon_viewing_angles = np.linspace(0,360,nframe)#[0,30,60,90,120,150,180,210,240,270]#-73

    # latitude/longitude line vectors
    lat_line_range = [-90,90]
    lat_lines = 8
    lat_line_count = (lat_line_range[1]-lat_line_range[0])/lat_lines

    merid_range = [-180,180]
    merid_lines = 8
    merid_count = (merid_range[1]-merid_range[0])/merid_lines


    # define color maps for water and land
    #ocean_map = (plt.get_cmap('ocean'))(210)
    #cmap = plt.get_cmap('gist_earth')



    # for making the gif animation
    gif_indx = 0

    os.system('rm ./png_dir/*_.png')


    for lonv in lon_viewing_angles:
        plt.cla()
        # call the basemap and use orthographic projection at viewing angle

        m = Basemap(projection='ortho', 
                  lat_0=lat_viewing_angle, lon_0=lonv)    

        xo,yo = m(dato['RA'],dato['DEC'])
        xno,yno = m(datno['RA'],datno['DEC'])

        #m.plot(x,y,'k,')
        plt.scatter(xno,yno,c=datno[args.comp],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1)
        plt.scatter(xo,yo,c=dato[args.comp],edgecolor='none',s=.1,vmax=vx,vmin=vm,zorder=1000)
        #plt.hexbin(x,y,dat['WEIGHT'],edgecolor='none',gridsize=360)
        plt.colorbar()
        # coastlines, map boundary, fill continents/water, fill ocean, draw countries
        #m.drawcoastlines()
        #m.drawmapboundary(fill_color=ocean_map)
        #m.fillcontinents(color=cmap(200),lake_color=ocean_map)
        #m.drawcountries()


        ct = 1001
        lats = np.arange(lat_line_range[0],lat_line_range[1],lat_line_count)
        m.drawparallels(lats,zorder=ct)
        ct += 1
        for lat in lats:
            plt.annotate(str(lat),xy=m(20,lat),xycoords='data',color='k',zorder=ct)
            ct += 1

        mers = np.arange(merid_range[0],merid_range[1],merid_count)
        m.drawmeridians(mers,zorder=ct)
        ct += 1
        for mer in mers:
            if mer < 0:
                mer = 360+mer
            plt.annotate(str(mer),xy=m(mer,-10),xycoords='data',color='k',zorder=ct)
            ct += 1

        # save figure at 150 dpi and show it
        #plt.savefig('orthographic_map_example_python.png',dpi=150,transparent=True)
        plt.title(tp+' '+args.comp)
        #plt.show()
        #plt.pause(0.01)
        # iterate to create the GIF animation
        gif_maker(outdir+'basemap_rotating_DESI'+tp+args.comp+'.gif','./png_dir/',gif_indx,len(lon_viewing_angles)-1,dpi=90)
        gif_indx+=1 
        print(gif_indx)