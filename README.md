# About mimmo
**mimmo** is a C++ library for Manipulation and Morphing of surface/volume meshes (Open source).

The main aim of **mimmo** is to provide a framework of tools able to handle a tessellated geometry and quickly modify it. 
It relies on Radial Basis Function (RBF) and extended Free Form Deformation (eFFD) techniques to achieve it. 

Many efforts are spent to make **mimmo** a flexible platform, able to adapt to any User problem/need, and open to any customization. Its inner structure 
is made by essential executable blocks (EB), each one providing a basic feature, such as I/O handling, manipulating, selecting 
sub-portion of target geometries, etc... . 

Blocks are connectable between each other through input/output ports, so that they can describe a morphing workflow, easily controllable by the User.

Further info available @ [mimmo website](http://optimad.github.io/mimmo/)
 
**mimmo** is developed and maintained by **Optimad engineering srl** and is being distributed under the GNU LESSER GENERAL PUBLIC LICENSE version 3.
