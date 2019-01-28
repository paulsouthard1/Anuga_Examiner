import anuga
import os
import anuga.file.sww as sww
import numpy as np
import argparse as ap
import pickle

parser = ap.ArgumentParser(description='Easily create ASCII Rasters from Anuga SWW Model Outputs')
parser.add_argument('Input_SWW',help='Full path of SWW file to create rasters from')
parser.add_argument('Output_DEM_Path',help='Full path of where to store output DEM')
args = parser.parse_args()

# Define function to create directory to store ASCII Rasters
def CreateDir(Input_SWW,Output_DEM_Path):
    filename = os.path.basename(Input_SWW)
    filename = os.path.splitext(filename)[0]
    path_final = Output_DEM_Path + "\\" + filename + "\\"
    if not os.path.exists(path_final):
    	os.mkdir(path_final)
    return path_final, filename

# Define function to extract ASCII Rasters from SWW file
def CreateRasts(Input_SWW,path_final,filename):
    quantities = ['stage','depth','momentum']
    for i in range(len(quantities)):
        print('Creating raster for quantity ' + quantities[i])
        outname = path_final + '\\' + str(filename) + "_" + str(quantities[i] + ".asc")
        anuga.sww2dem(Input_SWW,outname,quantity=quantities[i],reduction = max, cellsize = 1)



def Main(Input_SWW,Output_DEM_Path):
    # Create directory to store output rasters
    path_final, filename = CreateDir(Input_SWW,Output_DEM_Path)
    # Extract ASCII Rasters from SWW
    CreateRasts(Input_SWW,path_final,filename) 

Main(args.Input_SWW,args.Output_DEM_Path)
