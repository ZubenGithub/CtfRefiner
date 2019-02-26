"""
  StarRefiner.py

  After the per particle defocus has been estimated the general plane
  that all particles lie upon. Then particles that are some threshold
  outside that plane can be removed.

"""

from ReadStarFile import ReadStarFile, WriteStarFile
import argparse

import pandas as pd 

import numpy as np
import scipy.linalg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math

# arguments
# Add option to show scatterplot
# Add option to only show scatterplot when removing a particle

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
parser.add_argument('--star_file', help='Star file with micrograph data', required=True)
parser.add_argument('--MicrographX', help='Micrograph image X size (pixels). Default is for K3 image.', default='5760')
parser.add_argument('--MicrographY', help='Micrograph image Y size (pixels)', default='4092')
parser.add_argument('--threshold',
                  help='Threshold in Angstroms for filtering out "bad" particles. Should correspond to thickness of ice.',
                  default=(1000),
                  type=int)
parser.add_argument('--showplot',
                  help='Plot micrographs that contain any particles outside threshold. Useful for finding appropriate threshold',
                  action="store_true",
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

#    print(bad_df)
#    print(df)

    # For visual interpretation of points, uncomment
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
            ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r', s=50)
            ax.scatter(good[:,0], good[:,1], good[:,2], c='b', s=50)
            plt.xlabel('MicrographX')
            plt.ylabel('MicrographY')
            ax.set_zlabel('Defocus')
            ax.axis('equal')
            ax.axis('tight')
            ax.set_title(i) # Shows the micrograph title
            print(i) # Since the micrograph title is often long also print it in terminal
            plt.show()
    else:
      pass   

print(df.columns)
df.columns = ['rln' + str(col) for col in df.columns]
bad_df.columns = ['rln' + str(col) for col in bad_df.columns]
print(df.columns)

WriteStarFile(args.star_file + '_bad.star', bad_df)
WriteStarFile(args.star_file + '_good.star', df)
