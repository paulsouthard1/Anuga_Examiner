import numpy as np
import numpy.ma as ma
import os
import argparse

parser = argparse.ArgumentParser(description='Produce ASCII raster of roughness values a la Casas et al, 2010')
parser.add_argument('ZRef_File', help='ASCII Raster of D values to calculate Epsilon')
parser.add_argument('Depth_File', help='ASCII Raster of h values to caclulate Epsilon')
parser.add_argument('N_File', help='Name for ASCII Raster output of roughness values')
parser.add_argument('Background_N',help='Background Roughness value of channel, where there is no vegetation')
args = parser.parse_args()


# Define function to load Z ref values in
def LoadZRef(ZRef_File):
    # Load array of Zref
    zref = np.genfromtxt(ZRef_File,delimiter = ' ',skip_header = 6)
    # Remove all negative values and NoData
    zref[zref < 0] = 0
    return zref

# Define function to load depth values in
def LoadDepth(Depth_File):
    # Load array of depth
    depth = np.genfromtxt(Depth_File,delimiter = ' ',skip_header = 6)
    depth[depth < 0] = 0
    return depth

# Define function to calculate array of Epsilon values
def Epsilon(zref,depth):
    newzref = ma.masked_where(zref == 0,zref)
    epsilon = depth/zref
    epsilon = ma.masked_where(epsilon <= 0.2,epsilon)
    epsilon = ma.masked_where(epsilon >= 7,epsilon)
    return epsilon

# Define fucntion to calculate array of function of alpha and epsilon
def Func(alpha,epsilon):
    function = ma.masked_where(epsilon == 0, epsilon)
    function = 1+(alpha*(1/function)*np.log((np.cosh((1/alpha)-((1/alpha)*function)))/(np.cosh(1/alpha))))
    return function

# Define function to calculate array of roughness values
def NCalc(Cu,g,depth,function):
    n = (depth**(1/6)/((g**(1/2))*Cu*function))
    return n

# Define function to read in raster ASCII header and create output raster ASCII
def ReadWriteRaster(inraster,outraster,array):
    from numpy import savetxt
    with open(inraster,"r") as file:
        header = {}
        for i in range(6):
            header[i]=file.readline()
    file.close()
    with open(outraster,'w') as file:
        for i in range(6):
            file.write(header[i])
    file.close()
    savetxt('vals.txt',array)
    datafile = open('vals.txt','r')
    data = datafile.read()
    datafile.close()
    with open(outraster,'a') as file:
        file.write(data)
    file.close()

def Main(ZRef_File,Depth_File,N_File,Background_N):
    # Define constants
    Cu = 4.5
    a = 1
    g = 9.8 # m/s^2
    # Load in values from rasters as arrays
    zref = LoadZRef(ZRef_File)
    depth = LoadDepth(Depth_File)
    depth = depth[:-1,:-1]
    # Produced masked array of roughness
    eps = Epsilon(zref,depth)
    func = Func(a,eps)
    n = NCalc(Cu,g,depth,func)
    # Fill masked values with default of 0.04
    Background_N = float(Background_N)
    n = ma.filled(n,fill_value=Background_N)
    ReadWriteRaster(ZRef_File,N_File,n)

Main(args.ZRef_File,args.Depth_File,args.N_File,args.Background_N)


