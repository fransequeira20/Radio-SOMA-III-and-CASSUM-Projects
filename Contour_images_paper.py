# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:27:16 2022

@author: Francisco Sequeira
"""

import numpy as np
import pickle
from astropy import units as u
from astropy.io import fits 
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle, add_beam
from astropy.visualization import LogStretch, ImageNormalize
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib as mpl
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from Plot_label_mark import plot_label_mark
import aplpy

# To settle the direction of the ticks (direction=in) in the inset region
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

#Paths for the sources data (table_coord) and dictionary with the contour levels
table_coord = 'Data/Tables/SOMA_coordinates.txt'
t_coord = np.genfromtxt(table_coord, names=True, dtype=None, encoding=None)
contour_levels_dic ='Data/contour_levels_dic.p'
contour_levels_dic = pickle.load(open(contour_levels_dic, 'rb'))

#Paths for the FITS images
FITS_path = 'Data/Images_FITS'

#Color for the SOMA aperture
myblue ='#0000DD'

for i, source in enumerate(t_coord['source']): 
    if source == 'IRAS_22198':
    
        #Ploting the IR image
        IR_fits_file = FITS_path+'/'+source+'_SOFIA_37.1um_cal_ast.fits'
        
        hdu_IR = fits.open(IR_fits_file)[0]
        wcs_IR = WCS(hdu_IR.header)
            
        fig, ax = plt.subplots()
        
        ax = plt.subplot(projection=wcs_IR) 
        lon = ax.coords[0]
        lat = ax.coords[1]
        
        #For better imaging, LogNorm and vmin=0.05
        #Liu et al. 2020 settings for the IR images
        #plt.imshow(hdu_IR.data, cmap='gnuplot', origin='lower',
        #           interpolation='bicubic', vmin=0, vmax=1) 
        #Rosero et al 2019 setting for the images
        vmin, vmax = 0.05, 3
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())
        plt.imshow(hdu_IR.data, cmap='gray', origin='lower', norm=norm, 
                   interpolation='nearest')
        cbar=plt.colorbar()
        cbar.set_label('Flux (Jy)', family='sans-serif')
        if vmax == 3:    
            cbar.set_ticks([0.5, 1.0, 2.0, 3.0])
            cbar.set_ticklabels([0.5, 1.0, 2.0, 3.0])
        elif vmax == 12:
            cbar.set_ticks([1.0, 2.0, 4.0, 8.0, 12.0])
            cbar.set_ticklabels([1.0, 2.0, 4.0, 8.0, 12.0])
        elif vmax == 16:
            cbar.set_ticks([1.0, 2.0, 4.0, 8.0, 12.0, 16.0])
            cbar.set_ticklabels([1.0, 2.0, 4.0, 8.0, 12.0, 16.0])
        
        #Ploting the VLA data and contour maps
        VLA_C_fits_file = FITS_path+'/'+source+'_C.image.pbcor.fits'
        VLA_K_fits_file = FITS_path+'/'+source+'_K.image.pbcor.fits'
    
        hdu_VLA_C = fits.open(VLA_C_fits_file)[0]
        data_VLA_C = hdu_VLA_C.data.squeeze()
        header_C = hdu_VLA_C.header
        wcs_VLA_C = WCS(header_C).celestial
        
        hdu_VLA_K = fits.open(VLA_K_fits_file)[0]
        data_VLA_K = hdu_VLA_K.data.squeeze()
        header_K = hdu_VLA_K.header
        wcs_VLA_K = WCS(header_K).celestial
        
        #Setting the values for the beam size
        B_min_C, B_maj_C, B_PA_C = hdu_VLA_C.header['BMIN'], hdu_VLA_C.header['BMAJ'], hdu_VLA_C.header['BPA']
        
        B_min_K, B_maj_K, B_PA_K = hdu_VLA_K.header['BMIN'], hdu_VLA_K.header['BMAJ'], hdu_VLA_K.header['BPA']
        
        #Setting the contour levels for the C and K band data
        ax.contour(data_VLA_C, transform=ax.get_transform(wcs_VLA_C), 
                   levels=[element *t_coord['C_rms'][i]*1e-06 for element in contour_levels_dic[source]['Cont_Levels_SOMA']['C']], 
                   colors='red', linewidths=1, alpha=0.6, returnlevels=True)
    
        #SOMA regions for the sources
        c = SkyCoord(ra=t_coord['RA'][i], dec=t_coord['Dec'][i], unit=(u.hourangle, u.deg))
        c_fk5 = c.transform_to('fk5')
        ra=c_fk5.ra.degree
        dec=c_fk5.dec.degree
        print(source)
        print(c_fk5)
        
        #Spherical Circles with astropy using RA, Dec and the SOMA aperture for each region.
        #All the values are in degrees (°)
        SOMA_region = SphericalCircle((ra * u.deg, dec * u.deg), 
                                      (t_coord['aperture_SOMA_deg'][i]/3600) * u.degree, 
                                      edgecolor=myblue, facecolor='none', transform=ax.get_transform('fk5'))
        ax.add_patch(SOMA_region)
      
        #Box for the zoom images with astropy using RA, Dec for each region.
        #All the values are in degrees (°)
        if 'Zoom_parameters' in contour_levels_dic[source].keys():
            x_box=contour_levels_dic[source]['Zoom_parameters']['x']/3600
            y_box=contour_levels_dic[source]['Zoom_parameters']['y']/3600
            ra_box_dic=contour_levels_dic[source]['Zoom_parameters']['RA']
            dec_box_dic=contour_levels_dic[source]['Zoom_parameters']['Dec']
            
            c = SkyCoord(ra=ra_box_dic, dec=dec_box_dic, unit=(u.hourangle, u.deg))
            c_fk5 = c.transform_to('fk5')
            ra_box=c_fk5.ra.degree
            dec_box=c_fk5.dec.degree
            
            if source == 'IRAS_22198' or source =='IRAS_21391_BIMA2' or source =='IRAS_21391_BIMA3':
                zoom_region = Rectangle(((ra_box-x_box), (dec_box-y_box/2)), 
                                        width = 2*x_box, height = y_box, edgecolor=myblue, 
                                        facecolor='none',linestyle='--', transform=ax.get_transform('fk5'))
            elif source == 'L1206_A':
                zoom_region = Rectangle(((ra_box-x_box/2), (dec_box-y_box/1.8)), 
                                        width = 1.8*x_box, height = y_box*0.75, edgecolor=myblue, 
                                        facecolor='none',linestyle='--', transform=ax.get_transform('fk5'))
            else:
                zoom_region = Rectangle(((ra_box-x_box/2), (dec_box-y_box/2)), 
                                        width = x_box, height = y_box, edgecolor=myblue, 
                                        facecolor='none',linestyle='--', transform=ax.get_transform('fk5'))
            ax.add_patch(zoom_region)
        
        #Axis limits for the images
        x, y = wcs_IR.all_world2pix(ra, dec, 0)
        ax.set_xlim(x-t_coord['aperture_SOMA_pix'][i], x+t_coord['aperture_SOMA_pix'][i])
        ax.set_ylim(y-t_coord['aperture_SOMA_pix'][i], y+t_coord['aperture_SOMA_pix'][i])
        
        #Set the levels and rms for the contour maps
        if 'Zoom_parameters' in contour_levels_dic[source].keys():
            C_levels = [str(x) for x in contour_levels_dic[source]['Cont_Levels_Zoom']['C']] 
            C_rms = str(round(contour_levels_dic[source]['C_rms'],1))
            K_levels = [str(x) for x in contour_levels_dic[source]['Cont_Levels_Zoom']['K']] 
            K_rms = str(round(contour_levels_dic[source]['K_rms'],1))
    
            anot_C = 'Levels C-Band: [' + ', '.join(C_levels) + '] x ' + C_rms + r' $\mu$' + 'Jy beam' + r'$^{-1}$'
            anot_K = 'Levels K-Band: [' + ', '.join(K_levels) + '] x ' + K_rms + r' $\mu$' + 'Jy beam' + r'$^{-1}$'
        
            ax.annotate(anot_C + '\n' + anot_K, 
                    xy=(x-t_coord['aperture_SOMA_pix'][i]*1.2, y-t_coord['aperture_SOMA_pix'][i]*1.6), 
                    annotation_clip=False, fontsize=8, family='sans-serif')
        else:
            C_levels = [str(x) for x in contour_levels_dic[source]['Cont_Levels_SOMA']['C']] 
            C_rms = str(round(contour_levels_dic[source]['C_rms'],1))
            
            anot_C = 'Levels C-Band: [' + ', '.join(C_levels) + '] x ' + C_rms + r' $\mu$' + 'Jy beam' + r'$^{-1}$'
            
            ax.annotate(anot_C,
                    xy=(x-t_coord['aperture_SOMA_pix'][i]*1.05, y-t_coord['aperture_SOMA_pix'][i]*1.45), 
                    annotation_clip=False, fontsize=8, family='sans-serif')
            
        #Set the scale bar for the main image
        scale_arc = 2*t_coord['aperture_SOMA_deg'][i]/3.    # also I could choose it based on the radius like r/5
        scale = scale_arc/3600.     #  to deg
        label = (scale_arc*t_coord['Dist_SOMA'][i]*1000)    # AU
        print('Scale bar distance: ', label)
        fontprops = fm.FontProperties(size='10', family='sans-serif')
        scalebar = AnchoredSizeBar(ax.transData, 4, str(int(round(label,-1)))+' au', 'upper left', pad=1, 
                                   color='white', frameon=False, size_vertical=0.2, fontproperties=fontprops)
        ax.add_artist(scalebar)
        
        #Labels for the marks of the detections, see Contour_levels.py for details
        plot_label_mark(source_name=source, wcs_IR=wcs_IR, ax=ax)
        
        # showing the outflow direction
        if 'Outflow_parameters' in contour_levels_dic[source].keys():
            if source == 'IRAS_22172_MIR2':
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xb'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yb'],
                          contour_levels_dic[source]['Outflow_parameters']['x'],
                          contour_levels_dic[source]['Outflow_parameters']['y'], width=0.2, head_width=0.7, color='blue')  
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xr'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yr'],
                          -contour_levels_dic[source]['Outflow_parameters']['x'],
                          -contour_levels_dic[source]['Outflow_parameters']['y'], width=0.2, head_width=0.7, color='red') 
            elif source == 'IRAS_22198':
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xb'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yb'],
                          contour_levels_dic[source]['Outflow_parameters']['x'],
                          contour_levels_dic[source]['Outflow_parameters']['y'], width=0.5, head_width=1.5, color='red')  
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xr'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yr'],
                          -contour_levels_dic[source]['Outflow_parameters']['x'],
                          -contour_levels_dic[source]['Outflow_parameters']['y'], width=0.5, head_width=1.5, color='blue') 
            else:
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xb'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yb'],
                          contour_levels_dic[source]['Outflow_parameters']['x'],
                          contour_levels_dic[source]['Outflow_parameters']['y'], width=0.5, head_width=1.5, color='blue')  
                plt.arrow(x+contour_levels_dic[source]['Outflow_parameters']['xr'],
                          y+contour_levels_dic[source]['Outflow_parameters']['yr'],
                          -contour_levels_dic[source]['Outflow_parameters']['x'],
                          -contour_levels_dic[source]['Outflow_parameters']['y'], width=0.5, head_width=1.5, color='red') 
        if 'Outflow_parameters_2' in contour_levels_dic[source].keys():
            plt.arrow(x+contour_levels_dic[source]['Outflow_parameters_2']['xb'],
                          y+contour_levels_dic[source]['Outflow_parameters_2']['yb'],
                          contour_levels_dic[source]['Outflow_parameters_2']['x'],
                          contour_levels_dic[source]['Outflow_parameters_2']['y'], width=0.5, head_width=1.5, color='red')  
            plt.arrow(x+contour_levels_dic[source]['Outflow_parameters_2']['xr'],
                          y+contour_levels_dic[source]['Outflow_parameters_2']['yr'],
                          -contour_levels_dic[source]['Outflow_parameters_2']['x'],
                          -contour_levels_dic[source]['Outflow_parameters_2']['y'], width=0.5, head_width=1.5, color='blue')
        
        #Zoom image (inset region), made using Viviana's script with APLPY
        if 'Zoom_parameters' in contour_levels_dic[source].keys():
            if source == 'L1206_A' or source == 'IRAS_21391_BIMA3':
                subplot=[0.57,0.15,0.16,0.27]
            elif source == 'IRAS_22198':
                subplot=[0.26,0.15,0.13,0.22]
            else:
                subplot=[0.26,0.15,0.16,0.27]
            
            f2 = aplpy.FITSFigure(wcs_VLA_K, figure=fig,
                           subplot=subplot)
            
            f2.recenter(ra_box, dec_box, width=x_box, height=y_box)
            
            f2.show_contour(data_VLA_C, colors='red', linewidths=1, 
                            levels=[element *t_coord['C_rms'][i]*1e-06 for element in contour_levels_dic[source]['Cont_Levels_Zoom']['C']], 
                            returnlevels=True, alpha=0.6)      
            
            f2.show_contour(data_VLA_K, colors='tab:cyan', linewidths=1, 
                            levels=[element *t_coord['K_rms'][i]*1e-06 for element in contour_levels_dic[source]['Cont_Levels_Zoom']['K']], 
                            returnlevels=True, alpha=0.6)        
            
            #Set the scale bar for the zoom image (inset region)
            scale_arc_ins = contour_levels_dic[source]['Zoom_parameters']['x']/3.    # also I could choose it based on the radius like r/5
            scale = scale_arc_ins/3600.     #  to deg
            f2.add_scalebar(scale)
            label = (scale_arc_ins*t_coord['Dist_SOMA'][i]*1000)    # AU
            f2.scalebar.set_corner('top left')
            f2.scalebar.set_linestyle('solid')
            f2.scalebar.set_linewidth(2)
            f2.scalebar.set_label(str(int(round(label,-1)))+' au')
            f2.scalebar.set_font(size='6', weight='medium', \
                                  stretch='normal', family='sans-serif', \
                                  style='normal', variant='normal')
            f2.scalebar.set_color('black')
            
            #Set the beam in the inset region
            if source == 'IRAS_21391_BIMA2' or source == 'IRAS_22198':
                f2.add_beam(major= B_maj_C, minor= B_min_C, angle=B_PA_C, fill=True, color='red', alpha=0.6, corner='bottom right')
                f2.add_beam(major= B_maj_K, minor= B_min_K, angle=B_PA_K, fill=True, color='cyan', alpha=0.6, corner='bottom right')
            else: 
                f2.add_beam(major= B_maj_C, minor= B_min_C, angle=B_PA_C, fill=True, color='red', alpha=0.6)
                f2.add_beam(major= B_maj_K, minor= B_min_K, angle=B_PA_K, fill=True, color='cyan', alpha=0.6)
            
            #Set axis and ploting of the zoom image (inset region)
            f2.ticks.set_color('black')
            f2.tick_labels.hide()
            f2.axis_labels.hide()
            f2.tick_labels.set_font(size=10, weight='medium', stretch='normal', family='sans-serif', 
                                    style='normal', variant='normal')
        else:
            #Set the beam in the plot for sources without detection (no inset region)
            add_beam(ax, major= B_maj_C, minor= B_min_C, angle=B_PA_C, fill=True, color='red', alpha=0.6)
            add_beam(ax, major= B_maj_K, minor= B_min_K, angle=B_PA_K, fill=True, color='cyan', alpha=0.6)
         
        #Set title, axis and ploting
        if source == 'IRAS_22198':    
            ax.set_title(source.replace('_',' ')+'+6336')
        elif source == 'IRAS_22391_BIMA2':
            ax.set_title('IRAS 22391+5802 BIMA2')
        elif source == 'IRAS_22391_BIMA3':
            ax.set_title('IRAS 22391+5802 BIMA3')
        elif source == 'IRAS_22391_MIR48':
            ax.set_title('IRAS 22391+5802 MIR48')
        elif source == 'IRAS_22172_MIR1':
            ax.set_title('IRAS 22172+5549 MIR1')
        elif source == 'IRAS_22172_MIR2': 
            ax.set_title('IRAS 22172+5549 MIR2')
        elif source == 'IRAS_22172_MIR3':   
            ax.set_title('IRAS 22172+5549 MIR3')   
        else:
            ax.set_title(source.replace('_',' '))   
            
        lon.set_axislabel(r'Right Ascension (J2000)', family='sans-serif')
        lon.set_ticklabel(exclude_overlapping=True)
        lon.display_minor_ticks(True)
        lon.tick_params(direction='in', length=6, color='white')
        lat.set_axislabel(r'Declination (J2000)', family='sans-serif')
        lat.set_ticklabel(exclude_overlapping=True)
        lat.display_minor_ticks(True)
        lat.tick_params(direction='in', length=6, color='white')
        #plt.savefig(source+'_VLA_contours.png', bbox_inches='tight')
        #plt.savefig(source+'_VLA_contours.pdf')
        plt.show()
        plt.close()
    
    else:
        continue
