# Change Log
All notable changes to **mimmo** project will be documented in this file.

This library _tries_ to adhere to [Semantic Versioning](http://semver.org/).

## Unreleased
### Fixed
- fixed cmake retracking dependancies (as in case of bitpit)
- fixed BasicShape: XML interface
- fixed BasicShape: selections method
- fixed IOCGNS: prism element mapping from/to cgns data structure.
- fixed w/ CGNS compilation-> bump to 3.3.1 version
- fixed compilation with Openfoam search: find LABEL_SIZE, ARCH_OPTION and PRECISION directly from FOAM env variables.

### Added
- MRBF class: added heaviside functions to Mimmo Radial Basis Function object
- added variable radius inner limits for cylindrical and spherical BasicShapes and Lattices.
- added new elementary basic shape and lattice WEDGE (triangular prism)
- added Linear System Solver Class based on PETSc
- added Volume Mesh handler class MimmoFvMesh
- added new modules "system" and "propagators"
- added new selectors FVSelectionBy (Box,Sphere,Cylinder) suitable for MimmoFvMesh objects.
- added reading and writing STL solid names in multiSolid STL.

### Changed
- PropagateField classes: employing direct solution of Laplacian system as well as smoothing. Changed
                          also User Interface and parameters controlling field propagation. Moved into the
                          new module "propagators".
- MimmoPiercedVector: changes to class interface
- IOOFOAM class : OpenFoam class interface now depending on OpenFoam native libraries, instead on VTK.
                  User interface has changed. Compatibility is ensured fro 2.x, 3.x 4.x and 5.x OpenFoam
                  Foundation versions, and ESI-OpenFoam v1606+, v1806+.
- SkdTreeUtils : extractTarget method changed in algorithm and interface.



 ## [1.2.0] - 2017-10-24
### Added
- added new geohandlers:Reconstruct classes
- added new geohandlers:Extract classes
- added new geohandlers:Switch classes
- added new manipulators:PropagateField class. Class that propagates a scalar/vector field from boundaries to interior of a domain. Mesh Morphing note : it may directly apply a propagated deformation (vector) field on points of a volume/surface mesh.
- added new MimmoPiercedVector container for data fields attached to a MimmoObject
- added new class for I/O of MimmoPiercedVector iogeneric::GenericInputMPVData and iogeneric::GenericInputMPVData
- added mimmo native dump/restore format *.geomimmo (MimmoGeometry)
- added new option to assign a reference PID to a geometry by MimmoGeometry execution
- added new option to write pidded 3D surfaces in MultiSolid ASCII (MimmoGeometry)
- added new features to OBBox utility (multi geometry in input, option to switch from Oriented to Axis Aligned box calculation)
- added global expert mode MIMMO_EXPERT to override mandatory ports checking in execution of chains
- Ports management: created new static singleton PortManager for on-the fly ports registration
- Ports management: added a C macro REGISTER_PORT for compile-time ports registration
- added full support to non-homogeneous unstructured meshes in MimmoObject

### Fixed
- fixed minor bugs in FFDLattice manipulator
- fixed bug in CGNSPIDExtractor (PID assignment)
- fixed bugs in MeshSelection classes
- fixed warnings

### Changed
- update CreateSeedOnSurface block: added sensitivity map option to drive the seeding.
- modified ports : added new ports for MimmoPiercedVector / Dismissed some ports no more useful (eg. pair MimmoObject & data attached).
- iocgns_example_00001 now applies volume mesh morphing by propagation of a deformation field defined on boundaries.
- update bitpit objects and functions with bitipit 1.4 release
- Changed types of file(FileType enum) handled in I/O by MimmoGeometry

### Removed
- removed geohandlers:Overlap classes, wrapped in new added geohandlers: Reconstruct classes
- removed geohandlers:Split classes, substituted by new added geohandlers: Extract classes
- removed iogeneric:MultipleMimmoGeometries
- removed manipulators:MultiApply
- Ports management: enum PortType, containerTAG and datatypeTAG are dissmissed

## [1.1.0] - 2017-06-09
### Added
- This CHANGELOG file.
- Unique logger for **mimmo** added.
- Added warnings and errors for mandatory ports and connections not linked.
- Added directly apply for manipulator blocks.
- Added useful ports and options to global manipulators.

### Changed
- MimmoGeometry now works in three operational modes (read/write/converter).
- Changed the name of **mimmo** executable in **mimmo++**
- Changed command options of **mimmo++**. For further info launch mimmo++ --help
- Changed layout of mimmo XML-TUI dictionaries. See examples in samples folder.
- ScalingGeometry renamed in ScaleGeometry.
- README updated with new modifications.

### Fixed
- Fixed plot in execution options for blocks.
- Fixed absorb/flush xml parameters.
- Bug fixing in global manipulators.
- Bug fixing in BendGeometry manipulator (local coordinate).
- Bug fixed in read CGNS format files.
- Cmake now activates retracking of optional modules
- Fixed bad and missing documentation.

### Removed
- Removed forced bvTree and kdTree refresh.
- Removed redundant setGeometry methods.

## [1.0.0] - 2017-05-11
### Added
- First **mimmo** release!
