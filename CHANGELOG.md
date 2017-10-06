# Change Log
All notable changes to **mimmo** project will be documented in this file.

This library _tries_ to adhere to [Semantic Versioning](http://semver.org/).


## [Unreleased]
### Added
- added new option to assign a reference PID to a geometry by MimmoGeometry execution
- added new option to write pidded 3D surfaces in MultiSolid ASCII (MimmoGeometry)
- added mimmo native dump/restore format *.geomimmo (MimmoGeometry)
- added new features to OBBox utility (multi geometry in input, option to switch from Oriented to Axis Aligned box calculation)
- added global expert mode MIMMO_EXPERT to override mandatory ports checking in execution of chains
- added new MimmoPiercedVector container for data fields attached to a MimmoObject.
- added new ports for MimmoPiercedVector / Dismissed some ports no more useful.
- added new geohandlers:Switch classes
- Ports management: created new static singleton PortManager for on-the fly ports registration.
- Ports management: added a C macro REGISTER_PORT for compile-time ports registration.
- added new class for I/O of MimmoPiercedVector iogeneric::GenericInputMPVData and iogeneric::GenericInputMPVData

### Fixed
- fixed minor bugs in FFDLattice manipulator
- fixed coverity issues

### Changed
- update CreateSeedOnSurface block: added sensitivity map option to drive the seeding. 

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

