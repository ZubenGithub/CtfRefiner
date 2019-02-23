from ReadStarFile import ReadStarFile, WriteStarFile
import argparse

import pandas as pd 

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('--star_file', help='Star file with micrograph data', default='run_data_job203_10000lines.star')
parser.add_argument('--MicrographX', help='Micrograph image X size (pixels)', default='5760')
parser.add_argument('--MicrographY', help='Micrograph image Y size (pixels)', default='4092')
parser.add_argument('--pixel_to_angstrom', help='Pixel to angstrom conversion factor', default='0.95')
parser.add_argument('--threshold', help='Threshold for filtering out "bad" particles', default='400')

p = parser.parse_args()

# Star file read into a pandas dataframe
df = ReadStarFile(p.star_file)

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
    for i in data:
        if shortest_distance(i[0], i[1], i[2], C[0], C[1], -1, C[2]) > float(p.threshold):
            bad.append(i)
            df_loc = df.loc[(df['CoordinateX'] == i[0]) & (df['CoordinateY'] == i[1]) & (df['DefocusV'] == i[2])]
            df = (df.drop(df_loc.index))
            bad_df = bad_df.append(df_loc)
        else:
            good.append(i)

#    print(bad_df)
#    print(df)

    # For visual interpretation of points, uncomment
    '''
    # plot points and fitted surface
    good = np.array(good)
    bad = np.array(bad)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)

    ax.scatter(good[:,0], good[:,1], good[:,2], c='b', s=50)
    ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r', s=50)


    plt.xlabel('MicrographX')
    plt.ylabel('MicrographY')
    ax.set_zlabel('Defocus')
    ax.axis('equal')
    ax.axis('tight')
    plt.show()
    '''

WriteStarFile(p.star_file + '_bad.star', bad_df)
WriteStarFile(p.star_file + '_good.star', df)
