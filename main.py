# -*- coding: utf-8 -*-
"""
Cold Spray Impact Simulation - Abaqus Python Script
====================================================
Simulates a spherical aluminum particle impacting an aluminum substrate.

Model Parameters:
  - Particle radius   : 40 µm
  - Substrate size    : 250 µm x 250 µm (cylinder, revolved)
  - Impact velocity   : -700 m/s (Y-direction)
  - Initial temperature: 298 K
  - Step time period  : 60 ns (6e-8 s)
  - Material          : Aluminum (Johnson-Cook plasticity + Mie-Gruneisen EOS)

Usage:
  Run from Abaqus CAE: File > Run Script, or via command line:
    abaqus cae noGUI=cold_spray_simulation.py
"""

# ── Standard Abaqus imports ───────────────────────────────────────────────────
from abaqus import *
from abaqusConstants import *
import section
import regionToolset
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch

# =============================================================================
# 0.  FRESH MODEL DATABASE
# =============================================================================
Mdb()
mdl = mdb.models['Model-1']

# =============================================================================
# 1.  GEOMETRY – PARTICLE  (sphere, R = 40 µm, revolved half-circle)
# =============================================================================
PARTICLE_RADIUS = 4.0e-05   # 40 µm

s = mdl.ConstrainedSketch(name='__profile__', sheetSize=0.001)
s.sketchOptions.setValues(decimalPlaces=5)
s.setPrimaryObject(option=STANDALONE)

# Axis of revolution (construction line along Y-axis)
s.ConstructionLine(point1=(0.0, -0.0005), point2=(0.0, 0.0005))
g = s.geometry
s.FixedConstraint(entity=g[2])

# Full circle then trim to right-half arc
s.CircleByCenterPerimeter(
    center=(0.0, 0.0),
    point1=(PARTICLE_RADIUS * 0.75, -PARTICLE_RADIUS * 0.75)
)
s.CoincidentConstraint(entity1=s.vertices[0], entity2=g[2], addUndoState=False)
s.RadialDimension(
    curve=g[3],
    textPoint=(1.3e-04, -6.9e-05),
    radius=PARTICLE_RADIUS
)

# Diameter line (axis chord) for the revolve profile
s.Line(point1=(0.0,  PARTICLE_RADIUS), point2=(0.0, -PARTICLE_RADIUS))
s.VerticalConstraint(entity=g[4], addUndoState=False)
s.ParallelConstraint(entity1=g[2], entity2=g[4], addUndoState=False)
s.CoincidentConstraint(entity1=s.vertices[2], entity2=g[2], addUndoState=False)
s.CoincidentConstraint(entity1=s.vertices[3], entity2=g[2], addUndoState=False)

# Trim left half of the circle to leave a semicircle
s.autoTrimCurve(
    curve1=g[3],
    point1=(-PARTICLE_RADIUS * 0.99, 0.0)
)

# Create the particle part (3-D solid of revolution, 360°)
p_particle = mdl.Part(
    name='Particle',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY
)
p_particle.BaseSolidRevolve(sketch=s, angle=360.0, flipRevolveDirection=OFF)
s.unsetPrimaryObject()
del mdl.sketches['__profile__']

# =============================================================================
# 2.  GEOMETRY – SUBSTRATE  (annular cylinder revolved from rectangle)
# =============================================================================
# Outer dimensions: 250 µm radius × 250 µm depth
SUBSTRATE_RADIUS = 2.5e-04
SUBSTRATE_DEPTH  = 2.5e-04

s1 = mdl.ConstrainedSketch(name='__profile__', sheetSize=0.001)
s1.sketchOptions.setValues(decimalPlaces=5)
s1.setPrimaryObject(option=STANDALONE)

s1.ConstructionLine(point1=(0.0, -0.0005), point2=(0.0, 0.0005))
g1 = s1.geometry
s1.FixedConstraint(entity=g1[2])

# Rectangle: from axis (x=0) outward
s1.rectangle(point1=(0.0, 0.0), point2=(SUBSTRATE_RADIUS, -SUBSTRATE_DEPTH))
v1 = s1.vertices
s1.CoincidentConstraint(entity1=v1[0], entity2=g1[2], addUndoState=False)
s1.ObliqueDimension(
    vertex1=v1[3], vertex2=v1[0],
    textPoint=(1.0e-04, 5.0e-05), value=SUBSTRATE_RADIUS
)
s1.ObliqueDimension(
    vertex1=v1[2], vertex2=v1[3],
    textPoint=(3.2e-04, -6.8e-05), value=SUBSTRATE_DEPTH
)

p_substrate = mdl.Part(
    name='substrate',
    dimensionality=THREE_D,
    type=DEFORMABLE_BODY
)
p_substrate.BaseSolidRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
s1.unsetPrimaryObject()
del mdl.sketches['__profile__']

# =============================================================================
# 3.  MATERIAL – ALUMINUM (Johnson-Cook + Mie-Gruneisen EOS)
# =============================================================================
mat = mdl.Material(name='Aluminum')

mat.Conductivity(table=((237.2,),))
mat.Density(table=((2700.0,),))
mat.Elastic(type=SHEAR, table=((27.0e9,),))
mat.Eos(
    type=USUP,
    table=((5386.0, 1.339, 1.97),)   # c0, s, Gamma_0
)
mat.InelasticHeatFraction()

# Johnson-Cook: A, B, n, m, T_melt, T_ref
mat.Plastic(
    hardening=JOHNSON_COOK,
    scaleStress=None,
    table=((148.4e6, 345.5e6, 0.183, 0.895, 916.0, 298.0),)
)
# Johnson-Cook rate dependency: C, epsilon_dot_ref
mat.plastic.RateDependent(
    type=JOHNSON_COOK,
    table=((0.001, 1.0),)
)
mat.SpecificHeat(table=((898.2,),))

# Global constants required for thermal calculations
mdl.setValues(absoluteZero=0.0, stefanBoltzmann=6.67e-6)

# =============================================================================
# 4.  SECTION ASSIGNMENT
# =============================================================================
mdl.HomogeneousSolidSection(
    name='Section-1', material='Aluminum', thickness=None
)

for prt in (p_particle, p_substrate):
    cells = prt.cells.getSequenceFromMask(mask=('[#1 ]',),)
    region = prt.Set(cells=cells, name='Set-1')
    prt.SectionAssignment(
        region=region, sectionName='Section-1',
        offset=0.0, offsetType=MIDDLE_SURFACE,
        offsetField='', thicknessAssignment=FROM_SECTION
    )

# =============================================================================
# 5.  ASSEMBLY  (translate particle so it just touches the substrate top)
# =============================================================================
asm = mdl.rootAssembly
asm.DatumCsysByDefault(CARTESIAN)

inst_particle  = asm.Instance(name='Particle-1',  part=p_particle,  dependent=OFF)
inst_substrate = asm.Instance(name='substrate-1', part=p_substrate, dependent=OFF)

# Lift particle above substrate top face by one radius (40 µm)
asm.translate(instanceList=('Particle-1',), vector=(0.0, PARTICLE_RADIUS, 0.0))

# Make instances dependent for meshing from part mesh
asm.makeDependent(instances=(inst_particle, inst_substrate))

# =============================================================================
# 6.  ANALYSIS STEP  (Explicit coupled temp-displacement, 60 ns)
# =============================================================================
mdl.TempDisplacementDynamicsStep(
    name='Step-1',
    previous='Initial',
    timePeriod=6.0e-08,
    improvedDtMethod=ON
)

# Field output: stress, plastic strain, temperature, velocity, contact, status
mdl.fieldOutputRequests['F-Output-1'].setValues(
    variables=(
        'S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG',
        'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
        'NT', 'TEMP', 'HFL', 'RFL', 'EVF', 'STATUS'
    ),
    numIntervals=30
)

# =============================================================================
# 7.  MESH – PARTICLE
#     Partitioned into 8 octants for structured hex meshing
# =============================================================================
p_particle.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
p_particle.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
p_particle.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)

datums = p_particle.datums

# Three sequential partitions (XY, YZ, XZ planes)
p_particle.PartitionCellByDatumPlane(
    datumPlane=datums[3],
    cells=p_particle.cells.getSequenceFromMask(mask=('[#1 ]',),)
)
p_particle.PartitionCellByDatumPlane(
    datumPlane=datums[4],
    cells=p_particle.cells.getSequenceFromMask(mask=('[#3 ]',),)
)
p_particle.PartitionCellByDatumPlane(
    datumPlane=datums[5],
    cells=p_particle.cells.getSequenceFromMask(mask=('[#f ]',),)
)

p_particle.seedPart(size=2.5e-06, deviationFactor=0.1, minSizeFactor=0.1)
p_particle.generateMesh()

# Thermomechanical hex elements (Explicit)
eT1 = mesh.ElemType(
    elemCode=C3D8RT, elemLibrary=EXPLICIT,
    kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
    hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=ON
)
eT2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=EXPLICIT)
eT3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=EXPLICIT)

p_particle.setElementType(
    regions=(p_particle.cells.getSequenceFromMask(mask=('[#ff ]',),),),
    elemTypes=(eT1, eT2, eT3)
)

# =============================================================================
# 8.  MESH – SUBSTRATE
#     Partition: circular refinement zone under particle + octant planes
# =============================================================================
# --- Circular partition on top face ---
f_sub  = p_substrate.faces
e_sub  = p_substrate.edges
t_sub  = p_substrate.MakeSketchTransform(
    sketchPlane=f_sub[2], sketchUpEdge=e_sub[2],
    sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0)
)
s2 = mdl.ConstrainedSketch(
    name='__profile__', sheetSize=0.00141,
    gridSpacing=3e-05, transform=t_sub
)
s2.sketchOptions.setValues(decimalPlaces=5)
s2.setPrimaryObject(option=SUPERIMPOSE)
p_substrate.projectReferencesOntoSketch(sketch=s2, filter=COPLANAR_EDGES)
s2.CircleByCenterPerimeter(
    center=(0.0, 0.0),
    point1=(5.25e-05, -6.75e-05)   # ~85 µm radius refinement zone
)
p_substrate.PartitionFaceBySketch(
    sketchUpEdge=e_sub[2],
    faces=f_sub.getSequenceFromMask(mask=('[#4 ]',),),
    sketch=s2
)
s2.unsetPrimaryObject()
del mdl.sketches['__profile__']

# Extrude the circle edge downward to partition the solid
p_substrate.PartitionCellByExtrudeEdge(
    line=e_sub[2],
    cells=p_substrate.cells.getSequenceFromMask(mask=('[#1 ]',),),
    edges=(p_substrate.edges[1],),
    sense=FORWARD
)

# --- Octant datum plane partitions ---
p_substrate.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
p_substrate.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)

d_sub = p_substrate.datums
p_substrate.PartitionCellByDatumPlane(
    datumPlane=d_sub[5],
    cells=p_substrate.cells.getSequenceFromMask(mask=('[#3 ]',),)
)
p_substrate.PartitionCellByDatumPlane(
    datumPlane=d_sub[6],
    cells=p_substrate.cells.getSequenceFromMask(mask=('[#f ]',),)
)

# Coarse global seed, fine seed on impact-zone edges
p_substrate.seedPart(size=3.5e-05, deviationFactor=0.1, minSizeFactor=0.1)
p_substrate.seedEdgeBySize(
    edges=p_substrate.edges.getSequenceFromMask(
        mask=('[#ea1903f1 #63f ]',),
    ),
    size=3.5e-06, deviationFactor=0.1,
    minSizeFactor=0.1, constraint=FINER
)
p_substrate.generateMesh()

# Thermomechanical hex elements (Explicit – must match the step type)
eS1 = mesh.ElemType(
    elemCode=C3D8RT, elemLibrary=EXPLICIT,
    kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF,
    hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=ON
)
eS2 = mesh.ElemType(elemCode=C3D6T, elemLibrary=EXPLICIT)
eS3 = mesh.ElemType(elemCode=C3D4T, elemLibrary=EXPLICIT)

p_substrate.setElementType(
    regions=(p_substrate.cells.getSequenceFromMask(mask=('[#ff ]',),),),
    elemTypes=(eS1, eS2, eS3)
)

asm.regenerate()

# =============================================================================
# 9.  ALE ADAPTIVE MESHING  (particle + inner substrate cells)
# =============================================================================
cells_p  = asm.instances['Particle-1'].cells.getSequenceFromMask(
    mask=('[#ff ]',),
)
cells_s  = asm.instances['substrate-1'].cells.getSequenceFromMask(
    mask=('[#d2 ]',),
)
ale_region = regionToolset.Region(cells=cells_p + cells_s)

mdl.steps['Step-1'].AdaptiveMeshDomain(
    region=ale_region,
    controls=None,
    frequency=1,
    meshSweeps=3
)

# =============================================================================
# 10.  CONTACT INTERACTION
# =============================================================================
mdl.ContactProperty('Friction')
mdl.interactionProperties['Friction'].TangentialBehavior(
    formulation=PENALTY,
    directionality=ISOTROPIC,
    slipRateDependency=OFF,
    pressureDependency=OFF,
    temperatureDependency=OFF,
    dependencies=0,
    table=((0.3,),),                 # friction coefficient µ = 0.3
    shearStressLimit=None,
    maximumElasticSlip=FRACTION,
    fraction=0.005,
    elasticSlipStiffness=None
)
mdl.interactionProperties['Friction'].NormalBehavior(
    pressureOverclosure=HARD,
    allowSeparation=OFF,
    constraintEnforcementMethod=DEFAULT
)

# Contact surfaces
surf_particle = asm.Surface(
    side1Faces=asm.instances['Particle-1'].faces.getSequenceFromMask(
        mask=('[#cc4b0 ]',),
    ),
    name='m_Surf-1'
)
surf_substrate = asm.Surface(
    side1Faces=asm.instances['substrate-1'].faces.getSequenceFromMask(
        mask=('[#80808400 ]',),
    ),
    name='s_Surf-1'
)

mdl.SurfaceToSurfaceContactExp(
    name='Int-1',
    createStepName='Step-1',
    main=surf_particle,
    secondary=surf_substrate,
    mechanicalConstraint=KINEMATIC,
    sliding=FINITE,
    interactionProperty='Friction',
    initialClearance=OMIT,
    datumAxis=None,
    clearanceRegion=None
)

# =============================================================================
# 11.  BOUNDARY CONDITIONS
# =============================================================================
# Encastre: bottom and outer cylindrical faces of substrate
fixed_faces = asm.instances['substrate-1'].faces.getSequenceFromMask(
    mask=('[#20400090 ]',),
)
fixed_region = asm.Set(faces=fixed_faces, name='Set-2')
mdl.EncastreBC(
    name='BC-1',
    createStepName='Initial',
    region=fixed_region,
    localCsys=None
)

# =============================================================================
# 12.  PREDEFINED FIELDS (initial velocity + temperature)
# =============================================================================
particle_set = asm.instances['Particle-1'].sets['Set-1']
substrate_set = asm.instances['substrate-1'].sets['Set-1']

# Particle impact velocity (negative Y = downward)
mdl.Velocity(
    name='Predefined Field-1',
    region=particle_set,
    field='',
    distributionType=MAGNITUDE,
    velocity2=-700.0,           # m/s
    omega=0.0
)

# Initial temperatures – 298 K for both parts
mdl.Temperature(
    name='Predefined Field-2',
    createStepName='Initial',
    region=particle_set,
    distributionType=UNIFORM,
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
    magnitudes=(298.0,)
)
mdl.Temperature(
    name='Predefined Field-3',
    createStepName='Initial',
    region=substrate_set,
    distributionType=UNIFORM,
    crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
    magnitudes=(298.0,)
)

# =============================================================================
# 13.  JOB CREATION AND SUBMISSION
# =============================================================================
mdb.Job(
    name='cs1',
    model='Model-1',
    description='Cold spray Al particle on Al substrate, v=-700 m/s',
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
    userSubroutine='',
    scratch='',
    resultsFormat=ODB,
    numDomains=4,
    activateLoadBalancing=False,
    numThreadsPerMpiProcess=0,
    numCpus=4
)

mdb.jobs['cs1'].submit(consistencyChecking=OFF)
mdb.jobs['cs1'].waitForCompletion()

print("\n====================================")
print("  Job cs1 submitted successfully.")
print("  Output file: cs1.odb")
print("====================================\n")
