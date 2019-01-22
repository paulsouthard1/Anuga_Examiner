import anuga
import os
import anuga.file.sww as sww
import numpy as np
import argparse as ap
import pickle

parser = ap.ArgumentParser(description='Easily create ASCII Rasters from Anuga SWW Model Outputs')
parser.add_argument(Input_SWW,help='Full path of SWW file to create rasters from')
parser.add_argument(Output_DEM_Path,help='Full path of where to store output DEM')
parser.add_argument(Full_Series,help='"Yes" if you want to extract the full series and "No" if you only want the end result')
args = parser.parse_args()

# Define function to create directory to store ASCII Rasters
def CreateDir(Input_SWW,Output_DEM_path):
    filename = os.path.basename(Input_SWW)
    filename = os.path.splitext(filename)[0]
    path_final = Output_DEM_Path + "\\" + filename + "\\"
    os.mkdir(path_final)
    return path_final, filename

# # Define function to get Time info from SWW
# def GetTime(Input_SWW):
#   sww_file = sww.Read_sww(Input_SWW)
#   timeinfo = sww_file.time
#   return timeinfo

# Define function to extract ASCII Rasters from SWW file
def CreateRasts(Input_SWW,Output_DEM_Path,Full_Series):
    quantities = ['elevation','stage','depth','momentum']
    sww_file = sww.Read_sww(Input_SWW)
    timeinfo = sww_file.time
    if Full_Series == 'Yes':
        for i in range(len(quantities)):
            for j in range(len(timeinfo)):
                outname = Output_DEM_Path + '\\' + str(filename) + "_" + str(quantities[i] + "_" + str(int(timeinfo[j])) + ".asc")
                anuga.sww2dem(Input_SWW,outname,quantity=quantities[i],reduction = j, cellsize = 5)
    else:
        for i in range(len(quantities)):
            outname = Output_DEM_Path + '\\' + str(filename) + "_" + str(quantities[i] + "_" + str(int(timeinfo[j])) + ".asc")
            anuga.sww2dem(Input_SWW,outname,quantity=quantities[i],reduction = (len(timeinfo)-1), cellsize = 5)



Main(Input_SWW,Output_DEM_Path,Full_Series):
    # Create directory to store output rasters
    path_final, filename = CreateDir(Input_SWW,Output_DEM_path)
    # Extract ASCII Rasters from SWW
    CreateRasts(Input_SWW,Output_DEM_Path,Full_Series) 

Main(args.Input_SWW,args.Output_DEM_Path,args.Full_Series)
