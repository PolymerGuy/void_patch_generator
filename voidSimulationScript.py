# -*- coding: mbcs -*-
# Do not delete the following import lines
import numpy as np
from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior




# Various helper functions
def read_voids(path):
    return np.genfromtxt(path,skip_header=1,delimiter=",").tolist()


def void_is_within_bounds(void,domain_shape):
    xmin = domain_shape[0][0]
    ymin = domain_shape[0][1]

    xmax = domain_shape[1][0]
    ymax = domain_shape[1][1]

    void_x_min = void[0]-void[2]
    void_x_max = void[0]+void[2]
    void_y_min = void[1]-void[2]
    void_y_max = void[1]+void[2]

    if void_x_min<xmin or void_x_max>xmax or void_y_min<ymin or void_y_max>ymax:
        return False
    else:
        return True

def void_pressure_start_ygrad(void,domain_shape):
    ymin = domain_shape[0][1]
    ymax = domain_shape[1][1]
    heigth = ymax-ymin
    # Start from top
    return 1.-void[1]/heigth
    
#==============================================================================
#                                User defined settings
#==============================================================================
# Tolerance for picking surfaces
tol = 1e-5

# Void data
path_to_voids = "/home/sindreno/voidcloud/void_patch_generator/voids.csv"

# Job name
job_file_name = "propagation"

# Material settings
modulus = 1000.
pois_ratio = 0.45
yield_stress = 50.
hard_mod = 30.
locking_strain = 1.0
locking_modulus = 600

# Mesh settings
global_mesh_size = 5
void_mesh_size = 1

# Mass scaling
target_time_inc = 5.e-5

# Offset between void and particle
void_particle_offset = 1.e-3


# Domain size
domain_corners = ((-20,30),(1000,600))

# Pressure and deformation settings
# First external pressure is applied
external_pressure = 100.
pressurize_step_time = 1.0
# Then a global strain is applied
global_strain = 0.1
deformation_step_time = 1.0
# Then, the pressure is propagated into the specimen within a given region
void_pressure_mag = 100.
void_pressure_domain = ((420,0),(580,600))
void_pressure_rise_time = 0.05
propagate_step_time = 1.0
#==============================================================================
#                                Initial calculations
#==============================================================================
voids = read_voids(path_to_voids)
voids_with_pressure = [void_is_within_bounds(void,void_pressure_domain) for void in voids]
void_start_times = [void_pressure_start_ygrad(void,void_pressure_domain)*(propagate_step_time-void_pressure_rise_time) for void in voids]

#==============================================================================
#                                Make model
#==============================================================================
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

# Draw domain
s.rectangle(point1=domain_corners[0], point2=domain_corners[1])

# Draw the voids
for void in voids:
    x,y,radius = void
    # Outer surface
    s.CircleByCenterPerimeter(center=(x, y), point1=(x, y+radius))
    # Inner surface
    s.CircleByCenterPerimeter(center=(x, y), point1=(x, y+radius-void_particle_offset))

p = mdb.models['Model-1'].Part(name='Domain', dimensionality=TWO_D_PLANAR, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Domain']
p.BaseShell(sketch=s)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Domain']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']

#==============================================================================
#                         Define and assign material
#==============================================================================

# Material model
mdb.models['Model-1'].Material(name='SimplePoly')
mdb.models['Model-1'].materials['SimplePoly'].Density(table=((1e-09, ), ))
mdb.models['Model-1'].materials['SimplePoly'].Elastic(table=((modulus, pois_ratio), ))
mdb.models['Model-1'].materials['SimplePoly'].Plastic(table=((yield_stress, 0.0), (
    yield_stress+hard_mod*locking_strain, locking_strain), (yield_stress+hard_mod*locking_strain + locking_modulus*0.1, locking_strain+0.1)))

mdb.models['Model-1'].HomogeneousSolidSection(name='PolySection', 
    material='SimplePoly', thickness=None)

p = mdb.models['Model-1'].parts['Domain']
f = p.faces
all_faces = f.getByBoundingBox(domain_corners[0][0],domain_corners[0][1],0,domain_corners[1][0],domain_corners[1][1],0)
region = p.Set(faces=all_faces, name='all')
p = mdb.models['Model-1'].parts['Domain']
p.SectionAssignment(region=region, sectionName='PolySection', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

#==============================================================================
#                                   Assembly
#==============================================================================

a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Domain']
a.Instance(name='Domain-1', part=p, dependent=ON)

#==============================================================================
#                             Step definitions
#==============================================================================
mdb.models['Model-1'].ExplicitDynamicsStep(name='Pressurise', 
    previous='Initial', massScaling=((SEMI_AUTOMATIC, MODEL, AT_BEGINNING, 
    0.0, target_time_inc, BELOW_MIN, 0, 0, 0.0, 0.0, 0, None), ),timePeriod=pressurize_step_time)
mdb.models['Model-1'].ExplicitDynamicsStep(name='Deform', 
    previous='Pressurise',timePeriod=deformation_step_time)
mdb.models['Model-1'].ExplicitDynamicsStep(name='Pressure propagation', 
    previous='Deform',timePeriod=propagate_step_time)

#==============================================================================
#                             Assign loads and contact 
#==============================================================================

# Smooth step amplitude
mdb.models['Model-1'].SmoothStepAmplitude(name='InitialPressure', timeSpan=STEP, 
        data=((0.0, 0.0), (pressurize_step_time, 1.0)))

# Contact properties
mdb.models['Model-1'].ContactProperty('ContactProperties')
mdb.models['Model-1'].interactionProperties['ContactProperties'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
a = mdb.models['Model-1'].rootAssembly



# Loop over all voids to assign pressure and contact properties
for i,void in enumerate(voids):
    x,y,radius = void

    void_outer_name = 'void_outer%i'%(i)
    void_inner_name = 'void_inner%i'%(i)

    void_sat_time = void_start_times[i] + void_pressure_rise_time

# Smooth step amplitude
    mdb.models['Model-1'].SmoothStepAmplitude(name='void%ipressure'%i, timeSpan=STEP, 
        data=((void_start_times[i], 0.0), (void_sat_time, 1.0)))

    # Define inner and outer surfaces and store as surfaces
    mdb.models['Model-1'].rootAssembly.Surface(name=void_outer_name, side1Edges=
        mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
        ((x+radius,y,0), )))   

    mdb.models['Model-1'].rootAssembly.Surface(name=void_inner_name, side1Edges=
        mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
        ((x+radius-void_particle_offset,y,0), )))

    outer_region = mdb.models['Model-1'].rootAssembly.surfaces[void_outer_name]
    inner_region = mdb.models['Model-1'].rootAssembly.surfaces[void_inner_name]


    if void_pressure_mag>0. and voids_with_pressure[i]:
        # Apply the pressure
        mdb.models['Model-1'].Pressure(name=void_outer_name, 
            createStepName='Pressure propagation', region=outer_region, distributionType=UNIFORM, 
            field='', magnitude=void_pressure_mag, amplitude='void%ipressure'%i)

        mdb.models['Model-1'].Pressure(name=void_inner_name, 
            createStepName='Pressure propagation', region=inner_region, distributionType=UNIFORM, 
            field='', magnitude=void_pressure_mag, amplitude='void%ipressure'%i)

    # Add contact conditions
    mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='Int%i'%i, 
    createStepName='Pressurise', master = outer_region, slave = inner_region, 
    mechanicalConstraint=KINEMATIC, sliding=FINITE, 
    interactionProperty='ContactProperties', initialClearance=OMIT, 
    datumAxis=None, clearanceRegion=None)


# Add pressure to top surface
mdb.models['Model-1'].rootAssembly.Surface(name="top surface", side1Edges=
    mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
    ((domain_corners[1][0]-tol,domain_corners[1][1],0), )))   
top_surface = mdb.models['Model-1'].rootAssembly.surfaces["top surface"]
mdb.models['Model-1'].Pressure(name="Top surface pressure", 
    createStepName='Pressurise', region=top_surface, distributionType=UNIFORM, 
    field='', magnitude=external_pressure, amplitude='InitialPressure')


#==============================================================================
#                             Assign constraints
#==============================================================================
# Left edge is fixed
mdb.models['Model-1'].rootAssembly.Set(name="left edge", edges=
    mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
    ((domain_corners[0][0],domain_corners[0][1]+tol,0), )))   
region = mdb.models['Model-1'].rootAssembly.sets["left edge"]

mdb.models['Model-1'].DisplacementBC(name='FixedLeft', 
    createStepName='Pressurise', region=region, u1=0.0, u2=UNSET, 
    ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
    fieldName='', localCsys=None)

# Bottom edge is fixed
mdb.models['Model-1'].rootAssembly.Set(name="bottom edge", edges=
    mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
    ((domain_corners[0][0]+tol,domain_corners[0][1],0), )))   
region = mdb.models['Model-1'].rootAssembly.sets["bottom edge"]

mdb.models['Model-1'].DisplacementBC(name='Fixedbottom', 
    createStepName='Pressurise', region=region,u1=UNSET, u2=0.0, 
    ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
    fieldName='', localCsys=None)

# Rigth edge is fixed
mdb.models['Model-1'].rootAssembly.Set(name="rigth edge", edges=
    mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
    ((domain_corners[1][0],domain_corners[1][1]-tol,0), )))   
region = mdb.models['Model-1'].rootAssembly.sets["rigth edge"]

mdb.models['Model-1'].DisplacementBC(name='FixedRigth', 
    createStepName='Pressurise', region=region, u1=0.0, u2=UNSET, 
    ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
    fieldName='', localCsys=None)

# Rigth edge is moved to impose global strain
# Smooth step amplitude
mdb.models['Model-1'].SmoothStepAmplitude(name='DeformStep', timeSpan=STEP, 
        data=((0.0, 0.0), (deformation_step_time, 1.0)))

mdb.models['Model-1'].boundaryConditions['FixedRigth'].setValuesInStep(
    stepName='Deform', u1=global_strain*(domain_corners[1][0]-domain_corners[0][0]), amplitude='DeformStep')


if True:

    #==============================================================================
    #                             Generate mesh
    #==============================================================================

    session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
        bcs=OFF, predefinedFields=OFF, connectors=OFF)
    session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
        meshTechnique=ON)
    p = mdb.models['Model-1'].parts['Domain']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=OFF, 
        engineeringFeatures=OFF, mesh=ON)
    session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(
        meshTechnique=ON)
    session.viewports['Viewport: 1'].view.setValues(width=883.68, height=294.115, 
        viewOffsetX=27.3271, viewOffsetY=-3.16406)

    for i,void in enumerate(voids):
        x,y,radius = void

        # Outer surface
        p = mdb.models['Model-1'].parts['Domain']
        e = p.edges
        pickedEdges = mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
            (((x+radius,y,0)), ))
        p.seedEdgeBySize(edges=pickedEdges, size=void_mesh_size, deviationFactor=0.1, 
            constraint=FINER)

        # Inner surface
        pickedEdges = mdb.models['Model-1'].rootAssembly.instances['Domain-1'].edges.findAt(
            (((x+radius-void_particle_offset,y,0)), ))
        p.seedEdgeBySize(edges=pickedEdges, size=void_mesh_size, deviationFactor=0.1, 
            constraint=FINER)

        p = mdb.models['Model-1'].parts['Domain']


    # Global seed
    p = mdb.models['Model-1'].parts['Domain']
    p.seedPart(size=global_mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    elemType1 = mesh.ElemType(elemCode=CPE4R, elemLibrary=EXPLICIT, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT)

    p = mdb.models['Model-1'].parts['Domain']
    pickedRegions =(all_faces, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    a1 = mdb.models['Model-1'].rootAssembly
    a1.regenerate()
    a = mdb.models['Model-1'].rootAssembly

    #==============================================================================
    #                             Make job and run
    #==============================================================================
    mdb.Job(name='FirstTest', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE, 
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=8, 
        activateLoadBalancing=False, multiprocessingMode=THREADS, numCpus=8)

    # mdb.jobs['FirstTest'].submit(consistencyChecking=OFF)

    mdb.saveAs(
        pathName='/home/sindreno/localizationStudy/cleanscript/VoidTemplate')
