# ------------------------------------------------------------------
#+
#+ IMPLEMENTATION
#+
#+    Authors:  Zuben P. Brown & Shayan Huda Chowdhury
#+    Title:    Particle refiner for tilt data
#+    Version:  0.1 
#+    GitHub:   https://github.com/ZubenGithub/TiltRefiner
#+
#+ DESCRIPTION
#+
#+   Script for removing particles after per-particle defocus estimation that
#+    lie some threshold away from a plane defined by all particles.
# ------------------------------------------------------------------
# HISTORY
#    2019/02/26 -- 
# ------------------------------------------------------------------
#%
#% Usage: /path/to/python TiltRefiner.py StarFile.star
#%
#% Options:
#%                            (default)
#%  --star_file     None    Star file with micrograph data (default: None)
#%  --MicrographX   5760    Micrograph image X size (pixels). Default is for K3 image
#%  --MicrographY   4092    Micrograph image Y size (pixels).
#%  --threshold     1000    Threshold in Angstroms for filtering out "bad" particles. Should correspond to thickness of ice.
#%  --showploti     False   Plot micrographs that contain any particles outside threshold. Useful for finding appropriate threshold
#%  --test          False   Shows the fist 10 micrographs with particles outside the threshold. Should be useful for determining the threshold range.
#%
#% Comments
#%
#%      Finds a plane for all particles and excludes any that fall outside
#%      some threshold (that should be approx. ice thickness).
#%      Requires pyem package from: https://github.com/asarnow/pyem/wiki/Install-pyem
# ------------------------------------------------------------------

import argparse

import pandas as pd 

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

from pyem import star
import re

# arguments
# Add option to show scatterplot
# Add option to only show scatterplot when removing a particle

def ReadStarFile(star_file):
    """
    Reads a star file and returns a pandas dataframe.
        NB: I need to test the different programs that it can read.
    Args:
        Star file
    Returns:
        Pandas Dataframe
    """
    dataframe = star.parse_star(star_file)
    columns = dataframe.columns
    new_columns = []
    # Change column names to general column names (without numbers)
    # Just parses the header and finds the names
    for column in columns:
        no_number = column.split()[0]
        new_column_name = re.sub('rln','',no_number) #NB: removed _
        new_columns.append(new_column_name)

    # Renames the columns with the neater and nicer names
    dataframe.columns = new_columns
    return(dataframe)

def WriteStarFile(star_file, df):
    """
    Writes a star file from given pandas dataframe.
    
    Args:
        Star file destination, pandas dataframe
    """
    return star.write_star(star_file, df)

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument(
                    '--star_file',
                    help='Star file with micrograph data',
                    required=True)
parser.add_argument(
                    '--MicrographX',
                    help='Micrograph image X size (pixels). Default is for K3 image.',
                    default='5760')
parser.add_argument(
                    '--MicrographY',
                    help='Micrograph image Y size (pixels)',
                    default='4092')
parser.add_argument(
                    '--threshold',
                    help='''Threshold in Angstroms for filtering out "bad"
                     particles. Should correspond to thickness of ice.''',
                    default=(1000),
                    type=int)
parser.add_argument(
                    '--showplot',
                    help='''Plot micrographs that contain any particles outside
                     threshold. Useful for finding appropriate threshold and 
                    visually checking which particles are being removed.''',
                    action="store_true",
                    default=False)
parser.add_argument(
                    '--test',
                    help='''Shows the fist 10 micrographs with particles outside
                    the threshold. Should be useful for determining the threshold range.''',
                    action="store_true",
                    required=False,
                    default=False)

args = parser.parse_args()

# Star file read into a pandas dataframe
df = ReadStarFile(args.star_file)

"""
For some operating systems the pyem starfile reader includes an underscore in the name.
e.g., MicrographName is stored as _MicrographName
This means that df.sort_values doesn't work
"""
if df.columns[0].startswith('_') == True:
    df.columns=df.columns.str.replace('_','')
else:
    pass

# dataframe to store 'bad' particles
bad_df = pd.DataFrame()

# Specific values of the df are selected and converted to np array for convenience
arr = df.sort_values(by='MicrographName').loc[:, ['MicrographName', 'CoordinateX', 'CoordinateY', 'DefocusV']].values

# dictionary of micrographs
micrographs = {}
for i in arr:
    if i[0] not in micrographs:
        micrographs[i[0]] = []
    micrographs[i[0]].append(i[1:4])

# This is a counter for testing the threshold value
test_counter = 0

for i in micrographs:
    data = np.array(micrographs[i], dtype=np.float)

    x = data[:,0]
    # regular grid covering the domain of the data
    X,Y = np.meshgrid(np.arange(min(x)/1000, max(x)/1000, 0.5)*1000, np.arange(min(x)/1000, max(x)/1000, 0.5)*1000)
    XX = X.flatten()
    YY = Y.flatten()

    # best-fit linear plane
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients

    # evaluate it on grid
    Z = C[0]*X + C[1]*Y + C[2]
    print ('%f x + %f y + %f = z' % (C[0], C[1], C[2]))

    def shortest_distance(x1, y1, z1, a, b, c, d):  
        d = abs((a * x1 + b * y1 + c * z1 + d))  
        e = (math.sqrt(a * a + b * b + c * c)) 
        return d/e

    # classification
    good = []
    bad = []
    for j in data:
        if shortest_distance(j[0], j[1], j[2], C[0], C[1], -1, C[2]) > float(args.threshold):
            bad.append(j)
            df_loc = df.loc[(df['CoordinateX'] == j[0]) & (df['CoordinateY'] == j[1]) & (df['DefocusV'] == j[2])]
            df = (df.drop(df_loc.index))
            bad_df = bad_df.append(df_loc)
        else:
            good.append(j)

    # For visual interpretation of points, uncomment
    if args.test == True:
      args.showplot = True
    else:
      pass
    if args.showplot == True:    
        # plot points and fitted surface
        # Currently only showing removed particles
        good = np.array(good)
        bad = np.array(bad)
    
        if len(bad) == 0:  # This prevents an error if bad has length 0
            pass
        else:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
            try:
              ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r', s=50)
            except IndexError:
              pass
            try:
              ax.scatter(good[:,0], good[:,1], good[:,2], c='b', s=50)
            except IndexError:
              pass
            plt.xlabel('MicrographX')
            plt.ylabel('MicrographY')
            ax.set_zlabel('Defocus')
            ax.axis('equal')
            ax.axis('tight')
            ax.set_title(i) # Shows the micrograph title
            print(i) # Since the micrograph title is often long also print it in terminal
            plt.show()
            if args.test == True:
                test_counter=test_counter + 1
                if test_counter == 11:
                    exit()
                else:
                    pass
    else:
      pass   

"""
 Relion headers have the format of:
 _rlnHEADERNAME
 This adds the 'rln' string so that they are saved correctly
"""
df.columns = ['rln' + str(col) for col in df.columns]
bad_df.columns = ['rln' + str(col) for col in bad_df.columns]

WriteStarFile(args.star_file + '_threshold{0}_bad.star'.format(args.threshold), bad_df)
WriteStarFile(args.star_file + '_threshold{0}_good.star'.format(args.threshold), df)
