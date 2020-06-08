# Change Log
All notable changes to **mimmo** project will be documented in this file.

This library _tries_ to adhere to [Semantic Versioning](http://semver.org/).

## Unreleased
### Fixed
- fixed copy during compilation of binary samples additional files
- various bug fixes
### Added
- restored name in MimmoPiercedVector objects
- new methods to write geometries with several and different kind of fields (scalar/vector, defined on points/cells) added to BaseManipulation class
- added geometric tolerance to mimmo objects
- added use of geometric tolerance and cleaning flag during reading with MimmoGeometry objects
- Added RefineGeometry class to handle refinement of surface mesh: ternary and red-green engines provided for triangular based meshes
- Added IOWaveFrontOBJ to handle input/output of surface mesh in Wavefront OBJ format. Enabled class ManipulateWFOBJData to manipulate data attached to the OBJ mesh (texture, normal fields, cell groups)
- Added new feature to MRBF manipulator class, i.e. the possibility to handle RBF node set with variable support radii.
### Changed
- modified plotOptionalResults methods in some classes to use the new write function of base class
- Point Cloud geometries now build the cell too. The cells are defined as bitpit::VERTEX type elements.
- introducing SelectField blocks, instead of SwitchField blocks for selecting a field in a block chain
- reworked IOCLoudPoints class: interface enhanced
- general cmake enhanced to retrack more robustly mimmo external dependencies.
### Removed
- SwitchField classes removed.

 ## [1.3.0] - 2019-10-25
### Fixed
- fixed cmake retracking dependancies (as in case of bitpit v.1.6.0).
- fixed BasicShape: XML interface.
- fixed BasicShape: selections method.
- fixed IOCGNS: prism element mapping from/to cgns data structure.
- fixed w/ CGNS compilation-> bump to 3.3.1 version.
- fixed compilation with Openfoam search: find LABEL_SIZE, ARCH_OPTION and PRECISION directly from FOAM env variables.
- fixed ExtractFields: fixed execution of classes workflow.
- fixed SwitchFields: fixed execution of classes workflow.
- fixed SelectionByMapping: fixed execution of class workflow.

### Added
- MRBF class: added heaviside functions to Mimmo Radial Basis Function object.
- BasicShapes/Lattice classes: added variable radius inner limits for cylindrical and spherical shapes.
- BasicShapes/Lattice classes: added new elementary basic shape and lattice WEDGE (triangular prism).
- SurfaceTriangulator class: added block to handle mixed type 3D surface tessellation as homogeneous triangulations ("geohandlers" module).
- RefineGeometry class: added new class to refine surface meshes ("geohandlers" module).
- MimmoFVMesh class: prototype of Volume Mesh handler. Can be used also as stand alone executable block.
- FVSelectionBy (Box,Sphere,Cylinder) classes: added new geometry handlers/selectors suitable for MimmoFvMesh objects (geohandlers module).
- added new modules "propagators".
- added new kit of functions for handling Stencils of linear systems (propagators module).
- added solving of linear systems through LA module classes of bitpit (v1.6.0).
- PropagateField classes: handling propagation of bc scalar/vector field inside bulk volume mesh employing solution of Laplacian systems.
- GenericInput/OutputMPVData: added classes to read/write MimmoPiercedVector fields to file (ascii/binary).
- MeshChecker class : added new class to check Volume mesh quality ("utils" module).
- Project(Segment/3DCurve/Patch)OnSurface classes: added classes to project a segment/3DCurve/SurfacePatch onto a target surface mesh.
- Module class: calculate magnitude scalar field from a generic input vector field.

### Changed
- MimmoGeometry/MimmoObject class: added reading and writing STL solid names in multiSolid STL.
- MimmoObject: structural changes to class interface
- MimmoPiercedVector class: structural changes to class interface
- IOOFOAM class : OpenFoam class interface now depending on OpenFoam native libraries, instead on VTK.
                  User interface has changed. Compatibility is ensured fro 2.x, 3.x 4.x and 5.x OpenFoam
                  Foundation versions, and ESI-OpenFoam v1606+, v1806+.
- IOCGNS class : rework of I/O of the class. Reading multizone unstructured mesh.
- OBBox class : rework of basic workflow.
- SkdTreeUtils namespace : extractTarget method changed in algorithm and interface.
- Scalar/Vector MimmoPiercedVector field PORTS are now exchanging data by structure pointer (M_GDISPLS, M_FILTER, M_SCALARFIELD, M_VECTORFIELD).
- Coefficients PORT M_BCOEFFS in class manipulators::BendGeometry are now exchanging data by structure pointer.
- XML executable mimmo++ : tiny fixes to User interface.

### Removed
- CGNSPidExtractor class: Removed. Use SelectionByPID and SurfaceTriangulator instead.



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
