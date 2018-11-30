

import numpy as np
import os
import time
import sys
import anuga
import pickle
from scipy.interpolate import NearestNDInterpolator
from anuga.file.netcdf import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a,netcdf_float
from anuga.operators.vegetation_operator import Vegetation_operator


# Define function to read in veg rasters
def VegCreate(name_in,quantity_name):
    easting_min = None
    easting_max = None
    northing_min = None
    northing_max = None

    anuga.asc2dem(name_in + '.asc', use_cache=False, verbose=True)

    infile = NetCDFFile(name_in + '.dem', netcdf_mode_r)

    ncols = int(infile.ncols)
    nrows = int(infile.nrows)
    xllcorner = float(infile.xllcorner)  # Easting of lower left corner
    yllcorner = float(infile.yllcorner)  # Northing of lower left corner
    cellsize = float(infile.cellsize)
    NODATA_value = float(infile.NODATA_value)

    dem_elevation = infile.variables['elevation']

    # Assign default values
    if easting_min is None: easting_min = xllcorner
    if easting_max is None: easting_max = xllcorner + ncols*cellsize
    if northing_min is None: northing_min = yllcorner
    if northing_max is None: northing_max = yllcorner + nrows*cellsize

    y = np.arange(nrows,dtype=np.float)
    y = yllcorner + (nrows-1)*cellsize - y*cellsize

    x = np.arange(ncols,dtype=np.float)
    x = xllcorner + x*cellsize

    xx,yy = np.meshgrid(x,y)
    xx = xx.flatten()
    yy = yy.flatten()

    flag = np.logical_and(np.logical_and((xx <= easting_max),(xx >= easting_min)),
                           np.logical_and((yy <= northing_max),(yy >= northing_min)))

    dem = dem_elevation[:].flatten()

    id = np.where(flag)[0]
    xx = xx[id]
    yy = yy[id]
    dem = dem[id]

    data_flag = dem != NODATA_value
    data_id = np.where(data_flag)

    points =  np.zeros((len(data_id[0]),3))
    points[:,0] =   xx[data_id] - easting_min
    points[:,1] = yy[data_id] - northing_min
    points[:,2] = dem[data_id]

    x_ = np.arange(0, np.ceil(points[:,0].max()+1))
    y_ = np.arange(0, np.ceil(points[:,1].max()+1))

    grid_x, grid_y = np.meshgrid(x_,y_)
    grid_z = np.zeros_like(grid_x)

    dem=dem_elevation[:]
    dem=dem.data
    grid_z=dem

    points = np.vstack((grid_x.flatten() + easting_min,
                         grid_y.flatten() + northing_min,
                         grid_z.flatten())).T

    infile.close()

    interp = NearestNDInterpolator(points[:,0:2], points[:,2])

    coord = domain.get_centroid_coordinates(absolute=True)
    z_ = interp( coord )

    print 'Veg Load: saving interpolated file: ', name_in + '_interp.npy'
    np.save(name_in + '_interp.npy', z_)

# Set relevant parameters here
topo_in = 'Topo'
bound_in = 'extent'
res = 10
inlet_elev = 1422.3
inlet_stage = 1
t_to_print = 10
duration = 1000
inlet_btag = 1
ref_btag = 0
outlet_btag = 3

identifier = topo_in+'_'+bound_in+'_'+str(res)+'m_S'+str(inlet_stage)

#Create list of parameters
Parameter_names = ["topo_in:","bound_in:","res:","inlet_elev:","inlet_stage:","t_to_print:","duration:","inlet_btag:","outlet_btag:"]
Parameter_vals = [topo_in,bound_in,res,inlet_elev,inlet_stage,t_to_print,duration,inlet_btag,outlet_btag]
Parameters = ["","","","","","","","",""]
for i in np.arange(0,9):
    Parameters[i] = str(Parameter_names[i])+" "+str(Parameter_vals[i])

#Save Parameter List as Metadata
with open(identifier+'_par.txt', 'w') as f:
    for item in Parameters:
        f.write(str(item) + "\n")

#Save Parameter List as Pkl
fname = identifier+'_par.pkl'
pickle_out = open(fname,"wb")
pickle.dump(Parameters,pickle_out)

# Import Channel topography DEM
# Create DEM from asc data
anuga.asc2dem(topo_in + '.asc', use_cache=False, verbose=True)

# Create DEM from asc data
anuga.dem2pts(topo_in + '.dem', use_cache=False, verbose=True)

# Define boundaries for mesh
# Read in coordinate file
bounding_polygon = anuga.read_polygon(bound_in + '.csv')
# Create mesh with bounding polygon
domain = anuga.create_domain_from_regions(bounding_polygon,
                                    boundary_tags={'inlet': [inlet_btag],'exterior': [ref_btag],'outlet': [outlet_btag]},
                                    maximum_triangle_area=res,
                                    mesh_filename=identifier+'.msh',
                                    use_cache=False,
                                    verbose=True)

# Name domain and decide where to save
domain.set_name(identifier)
domain.set_datadir('.')


# Set quantities for domain
# Set elevation from topography file
domain.set_quantity('elevation', filename=topo_in+'.pts')
# Set Manning Roughness of bed
domain.set_quantity('friction',0.030)
# Set initial stage (if dry bed, set to elevation)
domain.set_quantity('stage',expression='elevation')

#Initiate Vegetation operator
op1 = Vegetation_operator(domain, use_diffusivity=False)

#Read in vegetation rasters
VegCreate('Veg_D','veg_diameter')
VegCreate('Veg_S','veg_spacing')

# Set Quantities
op1.set_veg_quantity('Veg_S', quantity_name='veg_spacing',convert_file=False, save_file= False, load_interp=True)
op1.set_veg_quantity('Veg_D', quantity_name='veg_diameter', convert_file= False, save_file= False, load_interp=True)



# Define and set boundaries
# Define transmissive boundary for downstream outlet
Bt = anuga.Transmissive_boundary(domain)    # Continue all values on boundary
#Define reflective boundary for channel edges
Br = anuga.Reflective_boundary(domain)
# Define Dirichlet boundary for upstream inlet flow
Bi = anuga.Dirichlet_boundary([inlet_elev + inlet_stage,0,0])
domain.set_boundary({'inlet': Bi,'exterior': Br,'outlet': Bt})



for t in domain.evolve(yieldstep=t_to_print, finaltime=duration):
    print domain.timestepping_statistics()

