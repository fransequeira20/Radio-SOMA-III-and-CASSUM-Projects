# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 11:19:43 2022

@author: tiqui
"""

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
import pylab as pl
import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np
from astropy.io import ascii
import math
from scipy.optimize import curve_fit

###
# Uploading the data from webpages and csv files
###

### For the whole data UCHII and HII
#data_CORNISH = 'https://cornish.leeds.ac.uk/cgi-bin/public/catalogue_query.py?title=hii_all&type=CSV&artifacts=0'
### For HII regions only (sample of 37 sources)
data_CORNISH_HII = 'https://cornish.leeds.ac.uk/cgi-bin/public/catalogue_query.py?title=hii_only&type=CSV&artifacts=0'
### For UCHII regions only (sample of 239 sources)
data_CORNISH = 'https://cornish.leeds.ac.uk/cgi-bin/public/catalogue_query.py?title=uchii_only&type=CSV&artifacts=0'

df_CORNISH = pd.read_csv(data_CORNISH, usecols=[0,3,4])
df_CORNISH.rename(columns = {'#Name':'Source_Name', 'RA_deg':'RA(J2000)_C','Dec_deg':'Dec(J2000)_C'}, inplace = True)

df_CORNISH_HII = pd.read_csv(data_CORNISH_HII , usecols=[0,3,4])
df_CORNISH_HII.rename(columns = {'#Name':'Source_Name', 'RA_deg':'RA(J2000)','Dec_deg':'Dec(J2000)'}, inplace = True)

path_CORNISH_HII_html = "Data/Tables/CORNISH_All_HII.html"
data_CORNISH_HII_html = pd.read_html(path_CORNISH_HII_html, header=0)
df_CORNISH_HII_html = data_CORNISH_HII_html[0]

path_CORNISH_III = "Data/Tables/CORNISH_All_UCHII.html"
data_CORNISH_III = pd.read_html(path_CORNISH_III, header=0)
df_CORNISH_III = data_CORNISH_III[0]

### CSV for the SOMA Survey
data_SOMA = 'Data/Tables/SOMA_final_table - MAiN_TABLE.csv'
df_SOMA = pd.read_csv(data_SOMA, usecols=[0,1,2])
df_SOMA = df_SOMA.dropna()
df_SOMA.rename(columns = {'DEC(J2000)':'Dec(J2000)'}, inplace = True)

###
# Changing units for the coordinates on both dataframes
###

for ind in df_CORNISH.index:
    ra = df_CORNISH['RA(J2000)_C'][ind]
    dec = df_CORNISH['Dec(J2000)_C'][ind]
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    ra = c.ra.to_string(u.hour)
    dec = c.dec.to_string(u.degree, alwayssign=True)
    df_CORNISH['RA(J2000)_C'][ind] = ra
    df_CORNISH['Dec(J2000)_C'][ind] = dec
    
for ind in df_CORNISH_HII.index:
    ra = df_CORNISH_HII['RA(J2000)'][ind]
    dec = df_CORNISH_HII['Dec(J2000)'][ind]
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    ra = c.ra.to_string(u.hour)
    dec = c.dec.to_string(u.degree, alwayssign=True)
    df_CORNISH_HII['RA(J2000)'][ind] = ra
    df_CORNISH_HII['Dec(J2000)'][ind] = dec    

###
# Working on the dataframes, joins and drop nan values and other values
###

df_CORNISH_SOMA = pd.merge(df_CORNISH_III, df_SOMA, on='RA(J2000)', how='inner', suffixes=('', '_y'))
df_CORNISH_SOMA.drop(df_CORNISH_SOMA.filter(regex='_y$').columns, axis=1, inplace=True)

df_CORNISH_alone = pd.merge(df_CORNISH_III, df_SOMA, on='RA(J2000)', how='left', indicator=True)
df_CORNISH_alone = df_CORNISH_alone[(df_CORNISH_alone['_merge'] != 'both')].reset_index(drop=True)
df_CORNISH_alone.drop(df_CORNISH_alone.filter(regex='_y$').columns, axis=1, inplace=True)
df_CORNISH_alone.drop(['SOMA SOURCE','_merge'], axis=1, inplace=True)
df_CORNISH_alone.dropna(subset=['Dist_h(kpc)'], inplace=True)

df = df_CORNISH_SOMA = pd.merge(df_CORNISH_alone, df_CORNISH, on='Source_Name', how='inner')

# Optimal aperture calculated using Zoie algorithm without the correction for crowdede regions
optimal_aperture_values = [93.0, 21.0, 93.0, 15.0, 24.0, 34.0, 28.0, 27.0, 33.0, 13.0, 65.0, 14.0, 38.0, 21.0, 13.0, 13.0, 74.0, 15.0, 14.0, 84.0, 0.0, 79.0, 17.0, 47.0, 43.0, 50.0, 31.0, 37.0, 92.0, 13.0, 24.0, 15.0, 14.0, 27.0, 44.0, 18.0, 94.0, 77.0, 14.0, 13.0, 77.0, 92.0, 18.0, 92.0, 91.0, 23.0, 38.0, 14.0, 50.0, 19.0, 15.0, 14.0, 16.0, 41.0, 14.0, 41.0, 34.0, 15.0, 13.0, 15.0, 78.0, 83.0, 65.0, 88.0, 85.0, 88.0, 89.0, 89.0, 95.0, 91.0, 91.0, 70.0, 14.0, 85.0, 94.0, 13.0, 96.0, 77.0, 17.0, 13.0, 17.0, 14.0, 15.0, 81.0, 31.0, 24.0, 15.0, 28.0, 94.0, 17.0, 41.0, 21.0, 20.0, 12.0, 37.0, 80.0, 41.0, 95.0, 18.0, 89.0, 24.0, 18.0, 44.0, 13.0, 93.0, 17.0, 18.0, 20.0, 17.0, 21.0, 79.0, 23.0, 21.0, 35.0, 25.0, 28.0, 18.0, 14.0, 21.0, 25.0, 21.0, 21.0, 92.0, 21.0, 71.0, 21.0, 82.0, 83.0, 86.0, 86.0, 86.0, 18.0, 17.0, 91.0, 13.0, 15.0, 14.0, 13.0, 70.0, 15.0, 53.0, 16.0, 16.0, 21.0, 13.0, 28.0, 78.0, 13.0, 14.0, 14.0, 14.0, 15.0, 18.0, 20.0, 45.0, 14.0, 94.0, 92.0, 89.0, 94.0, 84.0, 84.0, 84.0, 84.0, 85.0, 5.0, 85.0, 85.0, 86.0, 91.0, 91.0, 87.0, 91.0, 62.0, 20.0, 13.0, 15.0, 80.0, 14.0, 13.0, 85.0, 20.0, 15.0, 52.0, 83.0, 82.0, 87.0, 38.0, 88.0, 90.0, 18.0, 97.0, 88.0, 88.0, 68.0, 20.0, 13.0, 18.0, 31.0, 38.0, 87.0, 78.0, 14.0]
optimal_aperture_values_HII = [63.0, 35.0, 53.0, 89.0, 48.0, 89.0, 28.0, 30.0, 14.0, 28.0, 5.0, 82.0, 92.0, 25.0, 85.0, 42.0, 48.0, 33.0, 27.0, 54.0, 25.0, 27.0, 5.0, 33.0, 93.0, 5.0, 40.0, 86.0, 53.0, 5.0, 85.0, 28.0, 25.0, 93.0, 88.0, 60.0, 77.0]
df['Optimal_apert'] = optimal_aperture_values
df_CORNISH_HII['Optimal_apert'] = optimal_aperture_values_HII

# Remove the bad sources from the sample
list_bad_sources = ['G019.6087-00.2351', 'G019.6090-00.2313', 'G024.4721+00.4877', 'G032.7966+00.1909', 'G032.7982+00.1937', 'G038.6465-00.2260', 'G034.0901+00.4365', 'G011.0328+00.0274', 'G012.4317-01.1112', 'G018.1460-00.2839', 'G024.8497+00.0881', 'G030.6881-00.0718', 'G030.7532-00.0511', 'G034.2544+00.1460', 'G043.1460+00.0139', 'G043.1489+00.0130', 'G043.1520+00.0115', 'G043.1651-00.0283', 'G043.1652+00.0129', 'G043.1657+00.0116', 'G043.1665+00.0106', 'G043.1674+00.0128', 'G043.1677+00.0196', 'G043.1684+00.0124', 'G043.1699+00.0115', 'G043.1701+00.0078', 'G043.1706-00.0003', 'G043.1716+00.0001', 'G043.1720+00.0080', 'G043.1763+00.0248', 'G048.6099+00.0270', 'G049.4640-00.3511', 'G049.4891-00.3763', 'G049.4905-00.3688', 'G010.6297-00.3380', 'G014.5988+00.0198', 'G018.8250-00.4675', 'G023.8985+00.0647', 'G025.7157+00.0487', 'G026.0083+00.1369', 'G028.4518+00.0027', 'G043.1684+00.0087', 'G048.9296-00.2793']
list_bad_sources_index = []

for i in df.index:
    source = df['Source_Name'][i]
    if source in list_bad_sources:
        list_bad_sources_index.append(i)
        
df.drop(list_bad_sources_index, inplace=True)
df = df.reset_index(drop=True)

# Drop HII regions that don't have distance data, and bad sources (G045.4790+00.1294)
df_CORNISH_HII_html.dropna(subset=['Dist_h(kpc)'], inplace=True)
df_CORNISH_HII_html.drop([29], inplace=True)
df_CORNISH_HII_html = df_CORNISH_HII_html.reset_index(drop=True)

###
# Anglada plot
###

# Paths for the three different average methods
path_all = 'Data/Fitting Results'
table_name = path_all+'/table_aver_4p_bad_sources.txt'
table_UCHII_verf = ascii.read(table_name)

path_SOMA = 'Data/Fitting Results/average_goodmodels_SOMA_sources_linear_flu+bkg.txt'
SOMA_table = ascii.read(path_SOMA)

path_all = 'Data/Fitting Results'
table_name = path_all+'/table_aver_4p.txt'
table_HII = ascii.read(table_name)

path_Anglada = 'Data/Anglada Models/'

# Add results of the fitting (Lbol) to the dataframe
Source_Name = table_UCHII_verf['Source_Name']
df_Lbol = pd.DataFrame(Source_Name, columns=['Source_Name'])
df_Lbol['lbol'] = table_UCHII_verf['lbol']

df_UCHII = pd.merge(df, df_Lbol, on='Source_Name', how='inner', suffixes=('', '_y'))
df_UCHII.drop(df_UCHII.filter(regex='_y$').columns, axis=1, inplace=True)

# Remove bad Sources HII regions (G045.4790+00.1294)
table_HII.remove_row(19)

# Uploading data from models and SOMA data points
data_dict_Lit = {}
Lit_sources = []
infile='Data/Anglada Models/Reference_table.txt'

for line in open(infile, 'r'):
    l1 = line.split()
    if l1==[]: continue
    skipChars = ['#']
    if line[0] in skipChars: continue
    if not l1[0] in data_dict_Lit:       
            data_dict_Lit[l1[0]] = {}
            Lit_sources.append(l1[0])

    data_dict_Lit[l1[0]] = {'Dist':   l1[1],
                            'Flux':  l1[2],
                            'Radio_lum': l1[3],
                            'Lum_bol': l1[4],
                            'Pdot': l1[5],
                            'Refs':l1[6]}

Unresolved_Kurtz_94 = ['G10.841-2.592', 'G28.200-0.049', 'G48.61+0.02', 'G76.383-0.621',\
                      'G138.295+1.555', 'G189.030+0.784', 'G189.876+0.516',\
                      'G206.543-16.347']
scale= 1.36 
scale_k = 0.95
Anglada_Lbol = []
Anglada_Radio = []
Kur_Lbol = []
Kur_Radio = []

for source in Lit_sources:
    if  data_dict_Lit[source]['Refs'] == 'Anglada_95'   and data_dict_Lit[source]['Lum_bol'] != 'na':
        Anglada_Lbol.append(float(data_dict_Lit[source]['Lum_bol']))
        Anglada_Radio.append(float(data_dict_Lit[source]['Radio_lum'])/scale)

### Values for the Data points from the SOMA Radio I, II and III
# SOMA 1
lbol_S1    = [1.4e4, 3.9e4, 5.6e4, 8.7e4, 4.2e5, 6.9e4, 4.9e4, 2.1e5]
lbol_max_S1 = [2.1e4, 5.6e4, 8.6e4, 19e4, 7.0e5, 14e4, 7.3e4, 4.7e5]
lbol_min_S1 = [0.9e4, 2.7e4, 3.7e4, 4.0e4, 2.5e5, 3.3e4, 3.3e4, 0.9e5]
Radio_flux_S1 = [0.14, 0.77, 1.15, 0.74, 91.0, 0.06, np.nan, 0.42]
Dist_S1 = [2.0, 2.0, 1.68, 2.2, 8.4, 1.64, 0.7, 2.65]

# SOMA 2
lbol_S2 =     [5.5e5, 5.5e5, 3.0e5, 2.6e5, 2.6e5, 1.2e5, 1.4e5]
lbol_max_S2 = [9.7e5, 9.4e5, 3.0e5, 6.6e5, 4.6e5, 3.4e5, 5.1e5]
lbol_min_S2 = [3.1e5, 3.3e5, 2.9e5, 1.0e5, 1.4e5, 0.4e5, 0.4e5]
Radio_flux_S2 = [1702.69, np.nan, 89.70, 8.53, np.nan, 0.04, 0.78]
Dist_S2 = [7.4, 5.5, 10.2, 1.7, 4.1, 5.55, 2.1]

# SOMA 3
lbol_S3 =     [1.5e4, 4.2e3, 1.9e3, 6.6e2, 4.7e3, 6.72e2, 6.3e2, 5.3e3, 5.5e2, 4.8e2, 4.4e2, 6.7e2]
lbol_min_S3 = [7.2e4, 11e3, 1.9e3, 14e2, 9.3e3, 6.7e2, 21e2, 34e3, 14e2, 12e2, 9.6e2, 6.7e2]
lbol_max_S3 = [0.3e4, 1.6e3, 1.9e3, 3.2e2, 2.3e3, 6.7e2, 1.9e2, 0.8e3, 2.1e2, 2.0e2, 2.0e2, 6.7e2]
Radio_flux_S3 = [260e-3, 233e-3, 634e-2, 153e-4, 313e-4, 264e-4, 258e-4, 276e-4, 272e-4, 160e-3, 161e-3, 178e-4]  
Dist_S3 = [1.8, 0.764, 0.39, 0.73, 0.776, 0.776, 2.40, 2.40, 2.40, 0.75, 0.75, 0.75]

# Anglada plot
pl.figure(9, figsize=(9,9))
Fig = 'Anglada_plot'
fig = pl.gcf()
ax = fig.add_axes([.15,.15, .8, .75])

# To make the Anglada and Kurtz points
pl.loglog(Anglada_Lbol, Anglada_Radio, '^y', markersize=8, markeredgewidth=1.5, alpha=0.5, 
          label = 'Jets low-Mass YSO: Anglada et al. 1995')           
#pl.loglog(Kur_Lbol, Kur_Radio, 'xk', markersize=8, markeredgewidth=1.5, alpha=0.5,
#          label = 'UC/HC HII: Kurtz et al. 1994')     

# Values from Monge
T_e = 1e4    #K
nu = 6.0       #GHz at 5 cm

# Lyman Continuum from YSO, making the teal (blue) line
infile_M603 = path_Anglada+'stel_Mc60.Sigma0.3.dat'
M_star_M603, r_star_M603, L_star_M603, T_star_M603, rad_lum_star_M603, Q_star_M603 = np.loadtxt(infile_M603, unpack=True ,usecols=[0, 1, 2, 3, 4, 5])

infile_M601 = path_Anglada+'stel_Mc60.Sigma1.dat'
M_star_M601, r_star_M601, L_star_M601, T_star_M601, rad_lum_star_M601, Q_star_M601 = np.loadtxt(infile_M601, unpack=True ,usecols=[0, 1, 2, 3, 4, 5])

infile_M602 = path_Anglada+'stel_Mc60.Sigma3.dat'
M_star_M602, r_star_M602, L_star_M602, T_star_M602, rad_lum_star_M602, Q_star_M602 = np.loadtxt(infile_M602, unpack=True ,usecols=[0, 1, 2, 3, 4, 5])

pl.loglog(L_star_M601, np.multiply(rad_lum_star_M601,1.07789612350129e-44), '-', color= 'teal', linewidth=3, label='_nolegend_', alpha=0.8)
pl.loglog(L_star_M602, np.multiply(rad_lum_star_M602,1.07789612350129e-44), '-', color= 'teal', linewidth=3, label='_nolegend_', alpha=0.8)
pl.loglog(L_star_M603, np.multiply(rad_lum_star_M603,1.07789612350129e-44), '-', color= 'teal', linewidth=3, label='_nolegend_', alpha=0.8)

# Lyman Continuum from ZAMS, making the black line
cluster_file = path_Anglada+'lbin_cluster_Lbol-Nlym.dat'
logL_bol_cl, logNe_05_cl, logNe_95_cl = np.loadtxt(cluster_file, unpack=True)

rad_lum_cesaroni_data = 2.08e-46 * 10**(np.array(logNe_95_cl))* nu**-0.1 * T_e**0.45 
pl.loglog(10**pl.array(logL_bol_cl), rad_lum_cesaroni_data, 'k-', linewidth=3)

pl.annotate(r'$S_{\nu}d^{2} =\,8 \times 10^{-3}(L_{bol})^{0.6}$', xy=(5e-1,5e-0), annotation_clip=False, fontsize='20')
pl.annotate(r'Lyman Continuum from ZAMS', xy=(4e-1,1e3), annotation_clip=False, fontsize='18')
pl.annotate(r'Lyman Continuum from YSO', xy=(1.5e2,2e-4), annotation_clip=False, fontsize='18', color='teal')

# CORNISH UCHII and HII data points
pl.scatter(x=df_UCHII['lbol'], y=df_UCHII['Flux (mJy)']*df_UCHII['Dist_h(kpc)']**2, marker = 'o', label='CORNISH HII', 
           color='blue', zorder=10, s=40)

pl.scatter(x=table_HII['lbol'], y=df_CORNISH_HII_html['Flux (mJy)']*df_CORNISH_HII_html['Dist_h(kpc)']**2, marker = 'o', 
           color='blue', zorder=10, s=40)

# To plot the SOMA data points with the errorbar (See Radio SOMA I, II and III)
sct1 = pl.scatter(x=lbol_S1, y=np.multiply(Radio_flux_S1,np.multiply(Dist_S1,Dist_S1)), marker = 'o', 
                  label='SOMA I and II', color='#d62728', zorder=5, s=80)
#plt.errorbar(sct1, xerr=[[np.subtract(lbol_S1,lbol_min_S1)], [np.subtract(lbol_max_S1,lbol_S1)]], color = 'grey', 
#             elinewidth=2, capsize=7, capthick=3, alpha=0.5)
sct2 = pl.scatter(x=lbol_S2, y=np.multiply(Radio_flux_S2,np.multiply(Dist_S2,Dist_S2)), marker = 'o', 
                  color='#d62728', zorder=5, s=80)
#plt.errorbar(sct2, xerr=[[np.subtract(lbol_S2,lbol_min_S2)], [np.subtract(lbol_max_S2,lbol_S2)]], color = 'grey', 
#             elinewidth=2, capsize=7, capthick=3, alpha=0.5)
sct3 = pl.scatter(x=lbol_S3, y=np.multiply(Radio_flux_S3,np.multiply(Dist_S3,Dist_S3)), marker = 'o', 
                  label='SOMA III', color='#8c92ac', zorder=5, s=80)
#plt.errorbar(sct3, xerr=[[np.subtract(lbol_S3,lbol_min_S3)], [np.subtract(lbol_max_S3,lbol_S3)]], color = 'grey', 
#             elinewidth=2, capsize=7, capthick=3, alpha=0.5)

# To plot the data points from Tanaka 2016
infile_A = path_Anglada+'A.dat'
Mst, Lst, Sdd_1GHz, Sdd_8GHz, Sdd_10GHz = np.loadtxt(infile_A, unpack=True ,usecols=[0, 3, 4, 5, 6])
pl.loglog(Lst, Sdd_8GHz, 's', color='gold', markersize=8, markeredgewidth=1.5, alpha=0.5, label = 'Model TTZ16') 

infile_Al = path_Anglada+'Al.dat'
Mst, Lst, Sdd_1GHz, Sdd_8GHz, Sdd_10GHz = np.loadtxt(infile_Al, unpack=True ,usecols=[0, 3, 4, 5, 6])
pl.loglog(Lst, Sdd_8GHz, 's', color='gold', markersize=8, markeredgewidth=1.5, alpha=0.5, label = '_nolegend_') 

infile_Ah = path_Anglada+'Ah.dat'
Mst, Lst, Sdd_1GHz, Sdd_8GHz, Sdd_10GHz = np.loadtxt(infile_Ah, unpack=True ,usecols=[0, 3, 4, 5, 6])
pl.loglog(Lst, Sdd_8GHz, 's', color='gold', markersize=8, markeredgewidth=1.5, alpha=0.5, label = '_nolegend_') 

infile_B = path_Anglada+'B.dat'
Mst, Lst, Sdd_1GHz, Sdd_8GHz, Sdd_10GHz = np.loadtxt(infile_B, unpack=True ,usecols=[0, 3, 4, 5, 6])
pl.loglog(Lst, Sdd_8GHz, 's', color='gold', markersize=8, markeredgewidth=1.5, alpha=0.5, label = '_nolegend_') 

infile_C = path_Anglada+'C.dat'
Mst, Lst, Sdd_1GHz, Sdd_8GHz, Sdd_10GHz = np.loadtxt(infile_C, unpack=True ,usecols=[0, 3, 4, 5, 6])
pl.loglog(Lst, Sdd_8GHz, 's', color='gold', markersize=8, markeredgewidth=1.5, alpha=0.5, label = '_nolegend_') 

pl.legend(fontsize=10, loc=0, ncol=1)

pl.xlabel('$L_\mathrm{bol}\,(L_\odot)$', fontsize='20')
pl.ylabel(r'$S_{\nu}D^2$ ($mJy$ $kpc^2$)', fontsize='20')
pl.gca().set_ylim(1e-4,1.0e6)
pl.gca().set_xlim(1e-1,1.9e6)

pl.gca().tick_params(axis='both', reset=True, length=10, width=1, which='major', direction='in', labelsize=18)
pl.gca().tick_params(axis='both', reset=True, length=4, width=1, which='minor', direction='in', labelsize=18)

## Theretical fit from the low mass sources (- line)
plot_lum_v2= pl.array(pl.gcf().gca().get_xlim())
rad_lum_an_v2= 8*10**(-3)*plot_lum_v2**(0.6)
pl.loglog(plot_lum_v2, rad_lum_an_v2,'--k')

pl.show()
#fig.savefig('Anglada Plot.png')#,dpi=300,bbox_inches='tight')
