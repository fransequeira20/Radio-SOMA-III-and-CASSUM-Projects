# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 15:56:49 2022

@author: Francisco Sequeira
"""

from astropy import units as u
from astropy.coordinates import SkyCoord

def plot_label_mark(source_name, wcs_IR, ax):
    if source_name == 'IRAS_22198':
        c = SkyCoord(ra='22:21:27.20', dec='+63:51:40.15', unit=(u.hourangle, u.deg))
        c_fk5 = c.transform_to('fk5')
        ra=c_fk5.ra.degree
        dec=c_fk5.dec.degree
        x, y = wcs_IR.all_world2pix(ra, dec, 0)
        ax.annotate('MM2', xy=(x, y), fontsize=8, color='blue')
        
        c = SkyCoord(ra='22:21:27.68', dec='+63:51:26.96', unit=(u.hourangle, u.deg))
        c_fk5 = c.transform_to('fk5')
        ra=c_fk5.ra.degree
        dec=c_fk5.dec.degree
        x, y = wcs_IR.all_world2pix(ra, dec, 0)
        ax.annotate('A', xy=(x, y), fontsize=8, color='blue')        
        
        c = SkyCoord(ra='22:21:27.48', dec='+63:51:14.92', unit=(u.hourangle, u.deg))
        c_fk5 = c.transform_to('fk5')
        ra=c_fk5.ra.degree
        dec=c_fk5.dec.degree
        x, y = wcs_IR.all_world2pix(ra, dec, 0)
        ax.annotate('B', xy=(x, y), fontsize=8, color='blue')

        c = SkyCoord(ra='22:21:24.88', dec='+63:51:50.98', unit=(u.hourangle, u.deg))
        c_fk5 = c.transform_to('fk5')
        ra=c_fk5.ra.degree
        dec=c_fk5.dec.degree
        x, y = wcs_IR.all_world2pix(ra, dec, 0)
        ax.annotate('C', xy=(x, y), fontsize=8, color='blue')        
        