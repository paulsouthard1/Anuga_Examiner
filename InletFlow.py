

import numpy
import os
import time
import sys
import anuga
import pickle

# Set relevant parameters here
topo_in = '11_fbe_usw'
bound_in = 'channel_upper'
res = 1
inlet_elev = 1422.3
inlet_stage = 1
t_to_print = 30
duration = 1500
inlet_btag = 0
outlet_btag = 8

identifier = topo_in+'_'+bound_in+'_'+str(res)+'m_S'+str(inlet_stage)

#Create list of parameters
Parameter_names = ["topo_in:","bound_in:","res:","inlet_elev:","inlet_stage:","t_to_print:","duration:","inlet_btag:","outlet_btag:"]
Parameter_vals = [topo_in,bound_in,res,inlet_elev,inlet_stage,t_to_print,duration,inlet_btag,outlet_btag]
Parameters = ["","","","","","","","",""]
for i in numpy.arange(0,9):
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
                                    boundary_tags={'inlet': [inlet_btag],'exterior': [1],'outlet': [outlet_btag]},
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

