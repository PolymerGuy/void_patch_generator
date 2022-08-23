from matplotlib import use
from matplotlib.pyplot import axis
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
    return np.genfromtxt(path, skip_header=1, delimiter=",").tolist()


def void_is_within_bounds(void, domain_shape):
    xmin = domain_shape[0][0]
    ymin = domain_shape[0][1]

    xmax = domain_shape[1][0]
    ymax = domain_shape[1][1]

    void_x_min = void[0] - void[2]
    void_x_max = void[0] + void[2]
    void_y_min = void[1] - void[2]
    void_y_max = void[1] + void[2]

    if void_x_min < xmin or void_x_max > xmax or void_y_min < ymin or void_y_max > ymax:
        return False
    else:
        return True


def void_pressure_start_ygrad(void, domain_shape):
    ymin = domain_shape[0][1]
    ymax = domain_shape[1][1]
    heigth = ymax - ymin
    # Start from top
    return 1.0 - void[1] / heigth


def edit_material_keywords(abq_model):
    abq_model.keywordBlock.synchVersions(storeNodesAndElements=False)
    model_keyword_blocks = abq_model.keywordBlock.sieBlocks
    for index, block in enumerate(model_keyword_blocks):
        print(index)
        print(block)
        if "User Material, constants=13" in block:
            abq_model.keywordBlock.insert(
                index - 1, "\n*INCLUDE,INPUT=./spm_vumat/depvar_spm.inc"
            )
            break
    return


# ==============================================================================
#                                User defined settings
# ==============================================================================
# Tolerance for picking surfaces
tol = 1e-5

# Void data
if axisymmetric_model:
    path_to_voids = "/home/jonas/git/void_patch_generator/voids_axisymmetric.csv"
else:
    path_to_voids = "/home/jonas/git/void_patch_generator/voids.csv"

# Model and job name
abq_model_name = "pvdf_SPM"
job_file_name = abq_model_name

# Model type (default=shell)
axisymmetric_model = False

# Material settings
use_vumat = True
if use_vumat:
    # define VUMAT material parameters
    mat_params = {
        # part A of the model
        "youngs_modulus": 3000.0,  # Young's modulus [MPa]
        "poissons_ratio": 0.3,  # Poisson's ratio [-]
        "eps_0": 1e-3,  # reference strain rate [1/s]
        "C_t": 7e-2,  # temperature-dependent strain-rate-sensitivity parameter [-]
        "sigma_t0": 46.8,  # initial yield stress in uniaxial tension [Mpa]
        "sigma_ts": 37.8,  # saturation yield stress in uniaxial tension [Mpa]
        "alpha": 1.0,  # ratio of uniaxial tensile and compressive yield stresses [-]
        "beta": 1.0,  # parameter controlling the plastic dilation [-]
        "H": 15.0,  # hardening modulus / ramping parameter [-]
        # part B of the model
        "C_r": 5.5,  # initial elastic modulus of part B of the model [MPa]
        "kappa": 0.0,  # bulk modulus for network contribution in part B of the model [Mpa]
        "lambda_l": 2.0,  # locking stretch [-]
        # ABAQUS related
        "tsf": 1,  # Zero time step factor [-]
    }
else:
    modulus = 1000.0
    pois_ratio = 0.45
    yield_stress = 50.0
    hard_mod = 30.0
    locking_strain = 1.0
    locking_modulus = 600


# Mesh settings
# global_mesh_size = 5
global_mesh_size = 50
# void_mesh_size = 1
void_mesh_size = 8

# Mass scaling
target_time_inc = 5.0e-5

# Offset between void and particle
void_particle_offset = 1.0e-3


# Domain size
domain_corners = ((0, 0), (1000, 600))

# Pressure and deformation settings
# First external pressure is applied
external_pressure = 100.0
pressurize_step_time = 1.0
# Then a global strain is applied
global_strain = 0.1
deformation_step_time = 1.0
# Then, the pressure is propagated into the specimen within a given region
void_pressure_mag = 100.0
void_pressure_domain = ((420, 0), (580, 600))
void_pressure_rise_time = 0.05
propagate_step_time = 1.0
# ==============================================================================
#                                Initial calculations
# ==============================================================================
voids = read_voids(path_to_voids)
voids_with_pressure = [
    void_is_within_bounds(void, void_pressure_domain) for void in voids
]
void_start_times = [
    void_pressure_start_ygrad(void, void_pressure_domain)
    * (propagate_step_time - void_pressure_rise_time)
    for void in voids
]

# ==============================================================================
#                                Make model
# ==============================================================================
try:
    del mdb.models[abq_model_name]
except:
    pass

abq_model = mdb.Model(modelType=STANDARD_EXPLICIT, name=abq_model_name)
s = abq_model.ConstrainedSketch(name="__profile__", sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

# Draw domain
s.rectangle(point1=domain_corners[0], point2=domain_corners[1])

# Draw the voids
for void in voids:
    x, y, radius = void
    # Outer surface
    s.CircleByCenterPerimeter(center=(x, y), point1=(x, y + radius))
    # Inner surface
    s.CircleByCenterPerimeter(
        center=(x, y), point1=(x, y + radius - void_particle_offset)
    )


p = abq_model.Part(name="Domain", dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
p = abq_model.parts["Domain"]
p.BaseShell(sketch=s)
s.unsetPrimaryObject()

if axisymmetric_model:
    p.setValues(space=AXISYMMETRIC, type=DEFORMABLE_BODY)
    abq_model.rootAssembly.regenerate()
else:
    pass

p = abq_model.parts["Domain"]
session.viewports["Viewport: 1"].setValues(displayedObject=p)
del abq_model.sketches["__profile__"]

# ==============================================================================
#                         Define and assign material
# ==============================================================================

# Material model
if use_vumat:
    # define SPM_PVDF
    abq_model.Material(name="SPM_PVDF")
    abq_model.materials["SPM_PVDF"].Density(table=((1.780e-9,),))
    abq_model.materials["SPM_PVDF"].UserMaterial(
        mechanicalConstants=(
            mat_params["youngs_modulus"],  # props(1)
            mat_params["poissons_ratio"],  # props(2)
            mat_params["eps_0"],  # props(3)
            mat_params["C_t"],  # props(4)
            mat_params["sigma_t0"],  # props(5)
            mat_params["C_r"],  # props(6)
            mat_params["lambda_l"],  # props(7)
            mat_params["alpha"],  # props(8)
            mat_params["beta"],  # props(9)
            mat_params["kappa"],  # props(10)
            mat_params["sigma_ts"],  # props(11)
            mat_params["H"],  # props(12)
            mat_params["tsf"],  # props(13)
        ),
        thermalConstants=(),
        type=MECHANICAL,
        unsymm=OFF,
    )

    abq_model.HomogeneousSolidSection(
        name="SPM_PVDF", material="SPM_PVDF", thickness=None
    )
else:
    abq_model.Material(name="SimplePoly")
    abq_model.materials["SimplePoly"].Density(table=((1e-09,),))
    abq_model.materials["SimplePoly"].Elastic(table=((modulus, pois_ratio),))
    abq_model.materials["SimplePoly"].Plastic(
        table=(
            (yield_stress, 0.0),
            (yield_stress + hard_mod * locking_strain, locking_strain),
            (
                yield_stress + hard_mod * locking_strain + locking_modulus * 0.1,
                locking_strain + 0.1,
            ),
        )
    )

    abq_model.HomogeneousSolidSection(
        name="PolySection", material="SimplePoly", thickness=None
    )

p = abq_model.parts["Domain"]
f = p.faces
all_faces = f.getByBoundingBox(
    domain_corners[0][0],
    domain_corners[0][1],
    0,
    domain_corners[1][0],
    domain_corners[1][1],
    0,
)
region = p.Set(faces=all_faces, name="all")
p = abq_model.parts["Domain"]

if use_vumat:
    p.SectionAssignment(
        region=region,
        sectionName="SPM_PVDF",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )
else:
    p.SectionAssignment(
        region=region,
        sectionName="PolySection",
        offset=0.0,
        offsetType=MIDDLE_SURFACE,
        offsetField="",
        thicknessAssignment=FROM_SECTION,
    )
# ==============================================================================
#                                   Assembly
# ==============================================================================

a = abq_model.rootAssembly
session.viewports["Viewport: 1"].setValues(displayedObject=a)
session.viewports["Viewport: 1"].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF
)
a = abq_model.rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = abq_model.parts["Domain"]
a.Instance(name="Domain-1", part=p, dependent=ON)

# ==============================================================================
#                             Step definitions
# ==============================================================================
abq_model.ExplicitDynamicsStep(
    name="Pressurise",
    previous="Initial",
    massScaling=(
        (
            SEMI_AUTOMATIC,
            MODEL,
            AT_BEGINNING,
            0.0,
            target_time_inc,
            BELOW_MIN,
            0,
            0,
            0.0,
            0.0,
            0,
            None,
        ),
    ),
    timePeriod=pressurize_step_time,
)
abq_model.ExplicitDynamicsStep(
    name="Deform", previous="Pressurise", timePeriod=deformation_step_time
)
abq_model.ExplicitDynamicsStep(
    name="Pressure propagation", previous="Deform", timePeriod=propagate_step_time
)

# ==============================================================================
#                             Field Output request
# ==============================================================================

abq_model.FieldOutputRequest(
    createStepName="Pressurise",
    name="F-Output-1",
    variables=(
        "S",
        "SVAVG",
        "PE",
        "PEVAVG",
        "PEEQ",
        "PEEQVAVG",
        "LE",
        "U",
        "V",
        "A",
        "RF",
        "CSTRESS",
        "EVF",
        "SDV",
    ),
)

# ==============================================================================
#                             History Output request
# ==============================================================================

abq_model.HistoryOutputRequest(
    createStepName="Pressurise", name="H-Output-1", variables=PRESELECT
)

# ==============================================================================
#                             Assign loads and contact
# ==============================================================================

# Smooth step amplitude
abq_model.SmoothStepAmplitude(
    name="InitialPressure",
    timeSpan=STEP,
    data=((0.0, 0.0), (pressurize_step_time, 1.0)),
)

# Contact properties
abq_model.ContactProperty("ContactProperties")
abq_model.interactionProperties["ContactProperties"].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT
)
a = abq_model.rootAssembly


# Loop over all voids to assign pressure and contact properties
for i, void in enumerate(voids):
    x, y, radius = void

    void_outer_name = "void_outer%i" % (i)
    void_inner_name = "void_inner%i" % (i)

    void_sat_time = void_start_times[i] + void_pressure_rise_time

    # Smooth step amplitude
    abq_model.SmoothStepAmplitude(
        name="void%ipressure" % i,
        timeSpan=STEP,
        data=((void_start_times[i], 0.0), (void_sat_time, 1.0)),
    )

    # Define inner and outer surfaces and store as surfaces
    abq_model.rootAssembly.Surface(
        name=void_outer_name,
        side1Edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
            ((x + radius, y, 0),)
        ),
    )

    abq_model.rootAssembly.Surface(
        name=void_inner_name,
        side1Edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
            ((x + radius - void_particle_offset, y, 0),)
        ),
    )

    outer_region = abq_model.rootAssembly.surfaces[void_outer_name]
    inner_region = abq_model.rootAssembly.surfaces[void_inner_name]

    if void_pressure_mag > 0.0 and voids_with_pressure[i]:
        # Apply the pressure
        abq_model.Pressure(
            name=void_outer_name,
            createStepName="Pressure propagation",
            region=outer_region,
            distributionType=UNIFORM,
            field="",
            magnitude=void_pressure_mag,
            amplitude="void%ipressure" % i,
        )

        abq_model.Pressure(
            name=void_inner_name,
            createStepName="Pressure propagation",
            region=inner_region,
            distributionType=UNIFORM,
            field="",
            magnitude=void_pressure_mag,
            amplitude="void%ipressure" % i,
        )

    # Add contact conditions
    abq_model.SurfaceToSurfaceContactExp(
        name="Int%i" % i,
        createStepName="Pressurise",
        main=outer_region,
        secondary=inner_region,
        mechanicalConstraint=KINEMATIC,
        sliding=FINITE,
        interactionProperty="ContactProperties",
        initialClearance=OMIT,
        datumAxis=None,
        clearanceRegion=None,
    )


# Add pressure to top surface
abq_model.rootAssembly.Surface(
    name="top surface",
    side1Edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
        ((domain_corners[1][0] - tol, domain_corners[1][1], 0),)
    ),
)
top_surface = abq_model.rootAssembly.surfaces["top surface"]
abq_model.Pressure(
    name="Top surface pressure",
    createStepName="Pressurise",
    region=top_surface,
    distributionType=UNIFORM,
    field="",
    magnitude=external_pressure,
    amplitude="InitialPressure",
)


# ==============================================================================
#                             Assign constraints
# ==============================================================================
# Left edge is fixed
abq_model.rootAssembly.Set(
    name="left edge",
    edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
        ((domain_corners[0][0], domain_corners[0][1] + tol, 0),)
    ),
)
region = abq_model.rootAssembly.sets["left edge"]

abq_model.DisplacementBC(
    name="FixedLeft",
    createStepName="Pressurise",
    region=region,
    u1=0.0,
    u2=UNSET,
    ur3=UNSET,
    amplitude=UNSET,
    fixed=OFF,
    distributionType=UNIFORM,
    fieldName="",
    localCsys=None,
)

# Bottom edge is fixed
abq_model.rootAssembly.Set(
    name="bottom edge",
    edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
        ((domain_corners[0][0] + tol, domain_corners[0][1], 0),)
    ),
)
region = abq_model.rootAssembly.sets["bottom edge"]

abq_model.DisplacementBC(
    name="Fixedbottom",
    createStepName="Pressurise",
    region=region,
    u1=UNSET,
    u2=0.0,
    ur3=UNSET,
    amplitude=UNSET,
    fixed=OFF,
    distributionType=UNIFORM,
    fieldName="",
    localCsys=None,
)

# Rigth edge is fixed
abq_model.rootAssembly.Set(
    name="rigth edge",
    edges=abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
        ((domain_corners[1][0], domain_corners[1][1] - tol, 0),)
    ),
)
region = abq_model.rootAssembly.sets["rigth edge"]

abq_model.DisplacementBC(
    name="FixedRigth",
    createStepName="Pressurise",
    region=region,
    u1=0.0,
    u2=UNSET,
    ur3=UNSET,
    amplitude=UNSET,
    fixed=OFF,
    distributionType=UNIFORM,
    fieldName="",
    localCsys=None,
)

# Rigth edge is moved to impose global strain
# Smooth step amplitude
abq_model.SmoothStepAmplitude(
    name="DeformStep", timeSpan=STEP, data=((0.0, 0.0), (deformation_step_time, 1.0))
)

abq_model.boundaryConditions["FixedRigth"].setValuesInStep(
    stepName="Deform",
    u1=global_strain * (domain_corners[1][0] - domain_corners[0][0]),
    amplitude="DeformStep",
)


if True:

    # ==============================================================================
    #                             Generate mesh
    # ==============================================================================

    session.viewports["Viewport: 1"].assemblyDisplay.setValues(
        mesh=ON, loads=OFF, bcs=OFF, predefinedFields=OFF, connectors=OFF
    )
    session.viewports["Viewport: 1"].assemblyDisplay.meshOptions.setValues(
        meshTechnique=ON
    )
    p = abq_model.parts["Domain"]
    session.viewports["Viewport: 1"].setValues(displayedObject=p)
    session.viewports["Viewport: 1"].partDisplay.setValues(
        sectionAssignments=OFF, engineeringFeatures=OFF, mesh=ON
    )
    session.viewports["Viewport: 1"].partDisplay.meshOptions.setValues(meshTechnique=ON)
    session.viewports["Viewport: 1"].view.setValues(
        width=883.68, height=294.115, viewOffsetX=27.3271, viewOffsetY=-3.16406
    )

    for i, void in enumerate(voids):
        x, y, radius = void

        # Outer surface
        p = abq_model.parts["Domain"]
        e = p.edges
        pickedEdges = abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
            (((x + radius, y, 0)),)
        )
        p.seedEdgeBySize(
            edges=pickedEdges,
            size=void_mesh_size,
            deviationFactor=0.1,
            constraint=FINER,
        )

        # Inner surface
        pickedEdges = abq_model.rootAssembly.instances["Domain-1"].edges.findAt(
            (((x + radius - void_particle_offset, y, 0)),)
        )
        p.seedEdgeBySize(
            edges=pickedEdges,
            size=void_mesh_size,
            deviationFactor=0.1,
            constraint=FINER,
        )

        p = abq_model.parts["Domain"]

    # Global seed
    p = abq_model.parts["Domain"]
    p.setMeshControls(elemShape=QUAD, regions=all_faces)
    p.seedPart(size=global_mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    if axisymmetric_model:
        elemType1 = mesh.ElemType(
            elemCode=CAX4R,
            elemLibrary=STANDARD,
            secondOrderAccuracy=OFF,
            hourglassControl=DEFAULT,
            distortionControl=DEFAULT,
        )
        elemType2 = mesh.ElemType(elemCode=CAX3, elemLibrary=STANDARD)
    else:
        elemType1 = mesh.ElemType(
            elemCode=CPE4R,
            elemLibrary=EXPLICIT,
            secondOrderAccuracy=OFF,
            hourglassControl=DEFAULT,
            distortionControl=DEFAULT,
        )
        elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=EXPLICIT)

    p = abq_model.parts["Domain"]
    pickedRegions = (all_faces,)
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    a1 = abq_model.rootAssembly
    a1.regenerate()
    a = abq_model.rootAssembly

    # ==============================================================================
    #                             Make job and run
    # ==============================================================================
    if use_vumat:
        edit_material_keywords(abq_model)
        mdb.Job(
            name=job_file_name,
            model=abq_model.name,
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine="./spm_vumat/spm.f",
            scratch="",
            resultsFormat=ODB,
            parallelizationMethodExplicit=DOMAIN,
            numDomains=8,
            activateLoadBalancing=False,
            multiprocessingMode=THREADS,
            numCpus=8,
        )
    else:
        mdb.Job(
            name=job_file_name,
            model=abq_model.name,
            description="",
            type=ANALYSIS,
            atTime=None,
            waitMinutes=0,
            waitHours=0,
            queue=None,
            memory=90,
            memoryUnits=PERCENTAGE,
            explicitPrecision=SINGLE,
            nodalOutputPrecision=SINGLE,
            echoPrint=OFF,
            modelPrint=OFF,
            contactPrint=OFF,
            historyPrint=OFF,
            userSubroutine="",
            scratch="",
            resultsFormat=ODB,
            parallelizationMethodExplicit=DOMAIN,
            numDomains=8,
            activateLoadBalancing=False,
            multiprocessingMode=THREADS,
            numCpus=8,
        )
    # mdb.jobs['FirstTest'].submit(consistencyChecking=OFF)

    mdb.saveAs(
        pathName="/home/jonas/abaqus_work_dir/pressurised_voids/pressurised_voids"
    )
