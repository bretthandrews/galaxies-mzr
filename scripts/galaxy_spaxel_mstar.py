# @Author: Brett Andrews <andrews>
# @Date:   2018-03-28 17:03:95
# @Last modified by:   andrews
# @Last modified time: 2018-03-29 16:03:56

"""
Write the spaxel stellar mass of a galaxy to csv.

Spaxel stellar masses come from the Firefly stellar population catalog
of Goddard et al. (2017).

Firefly VAC description:
http://www.sdss.org/dr14/manga/manga-data/manga-firefly-value-added-catalog/

Datamodel:
https://data.sdss.org/datamodel/files/MANGA_FIREFLY/FIREFLY_VER/manga_firefly-STELLARPOP.html
"""

import argparse
import os
from os.path import join

from astropy.io import fits
from marvin import config
from marvin.tools.cube import Cube
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(description='Write resolved Mstar from Firefly stellar '
                                             'population catalog for a single galaxy.')
parser.add_argument('plateifu', type=str)
args = parser.parse_args()
# args = parser.parse_args('8485-1901')

path_repo = os.path.split(os.path.abspath('.'))[0]
path_data = join(path_repo, 'data')
path_firefly = join(path_data, 'manga_firefly-v2_1_2-STELLARPOP.fits')

fin = fits.open(path_firefly)
plateifus = fin['GALAXY_INFO'].data['PLATEIFU']

# Ngal x Nybin x Nxbin x
#     (binid,
#      xmin [arcsec],
#      xmax [arcsec],
#      ymin [arcsec],
#      ymax [arcsec],
#      image size [units of spaxel number])
spaxel_binid = fin['SPAXEL_BINID'].data

# Ngal x Nbin x (Mstar, Mstar_err)
mstar_all = fin['STELLAR_MASS_VORONOI'].data

# Select galaxy and binids
ind1 = np.where(plateifus == args.plateifu)[0][0]
ind_binid = spaxel_binid[ind1, :, :, 0].astype(int)

# Create 2D stellar mass array
mstar = np.ones(ind_binid.shape) * np.nan
for row, inds in enumerate(ind_binid):
    ind_nans = np.where(inds == -99)
    mstar[row] = mstar_all[ind1, inds, 0]
    mstar[row][ind_nans] = np.nan

# trim mstar to match size of DAP maps and write to csv
config.forceDbOff()
cube = Cube(args.plateifu)
len_x = int(cube.header['NAXIS1'])

df = pd.DataFrame(mstar[:len_x, :len_x])
fout = join(path_data, 'manga-{}_mstar.csv'.format(args.plateifu))
df.to_csv(fout, index=False)

print('\nWrote:', fout)

fin.close()
