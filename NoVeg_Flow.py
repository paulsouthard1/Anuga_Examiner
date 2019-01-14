#Import relevant modules
import numpy
import argparse
import os
import sys
import anuga
from anuga.shallow_water.boundaries import Inflow_boundary


parser = argparse.ArgumentParser(description='Run a simulation of flow over topography with uniform roughness')
parser.add_argument('topo_in', help='Name of topo data file w/o extension')
parser.add_argument('bound_in', help='Text file containing boundary coordinates w/o extension')
parser.add_argument('datadir', help='Directory to store Model Outputs')
parser.add_argument('res', help='Maximum Triangle Area for Mesh (resolution)')
parser.add_argument('t_to_print', help='Length of time before results are exported to model output file (s)')
parser.add_argument('duration', help='Length of simulation (s)')
parser.add_argument('inlet_btag', help='Boundary tag for inlet boundary')
parser.add_argument('outlet_btag', help='Boundary tag for outlet boundary')
parser.add_argument('manning', help='Manning Roughness of landscape')
parser.add_argument('q', help='Discharge at inlet (cms)')
parser.add_argument('slope', help='Average slope at inlet')
args = parser.parse_args()


def MetaData(topo_in,bound_in,res,t_to_print,duration,inlet_btag,outlet_btag,manning,q,slope):
    #Create list of parameters for metadata file for each model result
    Parameter_names = ["topo_in:","bound_in:","res:","t_to_print:","duration:","inlet_btag:","outlet_btag:","manning:","q:","slope:"]
    Parameter_vals = [topo_in,bound_in,res,t_to_print,duration,inlet_btag,outlet_btag,manning,q,slope]
    Parameters = ["","","","","","","","","",""]
    for i in numpy.arange(0,10):
        Parameters[i] = str(Parameter_names[i])+" "+str(Parameter_vals[i])

    #Save Parameter List as Metadata
    parlist = topo_in+'_'+str(q)+'cms_'+str(res)+'m'
    with open(parlist+'_par.txt', 'w') as f:
    	for item in Parameters:
            f.write(str(item) + "\n")
    print('Metadata file created')
    return parlist

def CreateDomain(topo_in,bound_in,datadir,inlet_btag,outlet_btag,res,parlist):
    # Import Channel topography DEM
    # Create DEM from asc data
    anuga.asc2dem(topo_in + '.asc', use_cache=False, verbose=True)

    # Create DEM from asc data
    anuga.dem2pts(topo_in + '.dem', use_cache=False, verbose=True)

    # Define boundaries for mesh
    # Read in coordinate file
    bounding_polygon = anuga.read_polygon(bound_in + '.csv')
    print(inlet_btag)
    #Determine tag to assign as reflective exterior
    for i in numpy.arange(1,4):
        print(i)
        if int(i) == int(inlet_btag):
            print(str(inlet_btag)+' already taken')
        else:
            print(str(i)+' used as exterior tag')
            break
    # Create mesh with bounding polygon
    domain = anuga.create_domain_from_regions(bounding_polygon,
                                        boundary_tags={'inlet': [inlet_btag],'exterior': [i],'outlet': [outlet_btag]},
                                        maximum_triangle_area=res,
                                        mesh_filename=parlist+'.msh',
                                        use_cache=False,
                                        verbose=True)

    # Name domain and decide where to save
    domain.set_name(parlist)
    domain.set_datadir('.')
    return domain

def SetQuant(topo_in,manning,domain):
    print('Setting Quantities')
    # Set quantities for domain
    # Set elevation from topography file
    domain.set_quantity('elevation', filename=topo_in+'.pts')
    # Set Manning Roughness of bed
    domain.set_quantity('friction',manning)
    # Set initial stage as dry bed
    domain.set_quantity('stage',expression='elevation')
    return domain

def SetBC(q,slope,domain):
    print('Setting BCs')
    # Define and set boundaries
    # Define transmissive boundary for downstream outlet
    Bt = anuga.Transmissive_boundary(domain)
    #Define reflective boundary for channel edges
    Br = anuga.Reflective_boundary(domain)
    # Define Dirichlet boundary for upstream inlet flow
    Bi = Inflow_boundary(domain,q,slope)
    domain.set_boundary({'inlet': Bi,'exterior': Br,'outlet': Bt})
    return domain

def Evolve(t_to_print,duration,domain):
    for t in domain.evolve(yieldstep=t_to_print, finaltime=duration):
        print domain.timestepping_statistics()

def Main(topo_in,bound_in,datadir,res,t_to_print,duration,inlet_btag,outlet_btag,manning,q,slope):
    parlist = MetaData(topo_in,bound_in,res,t_to_print,duration,inlet_btag,outlet_btag,manning,q,slope)
    res = int(res)
    t_to_print = int(t_to_print)
    duration = int(duration)
    inlet_btag = int(inlet_btag)
    outlet_btag = int(outlet_btag)
    q = float(q)
    slope = float(slope)
    domain = CreateDomain(topo_in,bound_in,datadir,inlet_btag,outlet_btag,res,parlist)
    domain = SetQuant(topo_in,manning,domain)
    domain = SetBC(q,slope,domain)
    print('Running Simulation')
    Evolve(t_to_print,duration,domain)

Main(args.topo_in,args.bound_in,args.datadir,args.res,args.t_to_print,args.duration,args.inlet_btag,args.outlet_btag,args.manning,args.q,args.slope)


