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

# Modules

## core 
<B>core</B> provides the basic structures and low level tools to manage them.
Main executable blocks are:

- <B>mimmo::BaseManipulation</B> is the base class and interface for each executable block
- <B>mimmo::MimmoObject</B> is the geometry container used to store and perform operation on a target geometry
- <B>mimmo::Chain</B> provides tools to manage and execute a set of manipulation blocks in an automatic work-flow

## iogeneric
<B>iogeneric</B> provides all executable blocks to read/write geometry from/on external files, as well as raw data in generic formats.
Main executable blocks are:

- <B>mimmo::MimmoGeometry</B> is the input/output block for geometry files. Supported format are:
  - <B>STL</B> Ascii/Binary triangulation
  - <B>VTU</B> Surface/Volume/Point Cloud/Curve
  - <B>NAS</B> Nastran triangulation
  - Ascii <B>OpenFoam</B> point cloud
- <B>mimmo::GenericInput</B> provides input interfaces to unformatted/CSV files
- <B>mimmo::GenericOutput</B> provides output interfaces to unformatted/CSV files

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

- <B>mimmo::GenericSelection</B> is the base block that provides interfaces to specific selection blocks:
  - <B>mimmo::SelectionByBox</B>
  - <B>mimmo::SelectionBySphere</B>
  - <B>mimmo::SelectionByCylinder</B>
  - <B>mimmo::SelectionByPID</B>
  - <B>mimmo::SelectionByMapping</B>
- <B>mimmo::ClipGeometry</B> and <B>mimmo::StitchGeometry</B> are executable blocks to cut a target geometry and stitch multiple geometries

## manipulators
<B>manipulators</B> provides core engines to morph a target mesh directly or with <B>eFFD</B> and <B>RBF</B> techniques.
Main executable blocks are:

- <B>mimmo::FFDLattice</B> provides extended Free Form Deformation based manipulation tools
- <B>mimmo::MRBF</B> provides Radial Basis Functions based manipulation tools
- <B>mimmo::RotationGeometry</B>, <B>mimmo::ScaleGeometry</B>, <B>mimmo::TranslationGeometry</B> and <B>mimmo::TwistGeometry</B> provide global direct transformation

## utils
<B>utils</B> provides miscellaneous features for morphing workflow

*/

