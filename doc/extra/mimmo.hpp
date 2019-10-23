/*! @mainpage mimmo

<B>mimmo</B> is a C++ library for Manipulation and Morphing of surface/volume meshes (Open source). \n
\n
The main aim of <B>mimmo</B> is to provide a framework of tools able to handle a tessellated geometry and quickly modify it. It relies on
Radial Basis Function (<B>RBF</B>) and extended Free Form Deformation (<B>eFFD</B>) techniques to achieve it. \n
\n
Many efforts are spent to make <B>mimmo</B> a flexible platform, able to adapt to any User problem/need, and open to any customization. Its inner structure
is made by essential executable blocks, each one providing a basic feature, such as I/O handling, manipulating, selecting
sub-portion of target geometries, etc... \n
\n
Blocks are connectable between each other, so that they can describe a morphing workflow, easily
controllable by the User. \n

# Quick contents description

## core
<B>core</B> provides the basic structures and low level tools to manage them.
Main executable blocks are:

- <B>mimmo::BaseManipulation</B> is the base class and interface for each executable block
- <B>mimmo::Chain</B> provides tools to manage and execute a set of manipulation blocks in an automatic work-flow
- <B>mimmo::MimmoObject</B> is the geometry container used to store and perform operation on a target geometry
- <B>mimmo::MimmoPiercedVector</B> is the container for data attached to a MimmoObject geometry
- <B>mimmo::Lattice</B> is the executable block to create structured mesh of points
- <B> miscellanous </B> utilities integrating CG,SkdTree and I/O VTU packages of bitpit


## iogeneric
<B>iogeneric</B> provides all executable blocks to read/write geometry from/on external files, as well as raw data in generic formats.
Main executable blocks are:

- <B>mimmo::MimmoGeometry</B> is the I/O handling block for geometry files.
- <B>mimmo::GenericInput</B> provides input interfaces to unformatted/CSV files
- <B>mimmo::GenericOutput</B> provides output interfaces to unformatted/CSV files
- <B>mimmo::GenericInputMPVData</B> provides input interfaces for MimmoPiercedVector data files
- <B>mimmo::GenericOutputMPVData</B> provides output interfaces for MimmoPiercedVector data files
- <B>mimmo::GenericDispls</B> interfaces for FFDLattice input dof displacements from file
- <B>mimmo::IOCloudPoints</B> interfaces for MRBF input nodes and dof displacements from file


## iocgns
<B>iocgns</B> collects all interfaces to volume mesh CGNS format input/output.

- <B>mimmo::IOCGNS</B> input/output block for CGNS meshes

## iovtk
<B>iovtk</B> collects all interfaces to volume mesh VTK format input/output.

- <B>mimmo::IOVTKScalar</B> input/output block for VTK surface meshes with scalar field

## ioofoam
<B>ioofoam</B> collects all interfaces to volume mesh OpenFoam format input/output.

- <B>mimmo::IOOFOAM</B> input/output block for meshes in OpenFOAM format (VTK surface and volume Point CLoud)

## geohandlers
<B>geohandlers</B> provides tools to handle, transform, select a target geometry mesh and attached data fields.
Main executable blocks are:

- <B>mimmo::...Selection...</B> different selection blocks from primitive-based(Box, Cylinder, etc..) inclusion selections
                              to PID or proximity Mapping.
- <B>mimmo::ClipGeometry</B> cutting-by-plane a target geometry
- <B>mimmo::Extract...</B> extract field subportion from a target MimmoPiercedVector field.
- <B>mimmo::Reconstruct...</B> reconstruct field on a target geometry from various MimmoPiercedVector fields defined on its subportions.
- <B>mimmo::RefineGeometry</B> refine a surface geometry into a finer triangular mesh.
- <B>mimmo::SurfaceTriangulator</B> triangulate a tessellated surface geometry
- <B>mimmo::StitchGeometry</B> stitch multiple geometries together

## manipulators
<B>manipulators</B> provides core engines to morph a target mesh directly or with <B>eFFD</B> and <B>RBF</B> techniques.
Main executable blocks are:

- <B>mimmo::FFDLattice</B> provides extended Free Form Deformation based manipulation tools
- <B>mimmo::MRBF</B> provides Radial Basis Functions based manipulation tools
- <B>mimmo::...Geometry</B> provide global direct transformation such as rotation, scaling, translation, twisting and bending
- <B>mimmo::Apply</B> applier of a deformation field to a target geometry

## propagators
<B>propagators</B> executable blocks for mesh morphing: provides laplacian
                   propagation of boundary information into a volume mesh


## utils
<B>utils</B> provides miscellaneous features useful to morphing workflow

*/
