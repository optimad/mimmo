/*! @mainpage mimmo

<B>mimmo</B> is a C++ library for Manipulation and Morphing of surface/volume meshes (Open source).
The main aim of mimmo is to provide a framework of tools able to handle a tessellated geometry and quickly modify it. It relies on 
Radial Basis Function (RBF) and extended Free Form Deformation (eFFD) techniques to achieve it. 
Many efforts are spent to make mimmo a flexible platform, able to adapt to any User problem/need, and open to any customization. Its inner structure 
is made by essential executable blocks (EB), each one providing a basic feature, such as I/O handling, manipulating, selecting 
sub-portion of target geometries, etc... . Blocks are connectable between each other, so that they can describe a morphing workflow, easily
controllable by the User.

# Modules

## core 
<B>core</B> provides the basic EB structure and low level tools to manage them.

## iogeneric
<B>iogeneric</B> provides all EBs to read/write geometry from/on external files, as well as raw data in generic formats.

## iocgns
<B>iocgns</B> collects all interfaces to volume mesh CGNS format input/output

## iovtk
<B>iovtk</B> collects all interfaces to volume mesh VTK format input/output.

## ioofoam
<B>ioofoam</B> collects all interfaces to volume mesh OpenFoam format input/output.

## geohandlers
<B>geohandlers</B> provides tools to handle, transform, select a target geometry mesh and scalar/vector field attached to it, if any

## manipulators
<B>manipulators</B> provides core engines to manipulate the target mesh with eFFD or RBF techniques

## utils
<B>utils</B> provides miscellaneous features for morphing workflow

*/

