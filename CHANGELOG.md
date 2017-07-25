# Change Log
All notable changes to **mimmo** project will be documented in this file.

This library _tries_ to adhere to [Semantic Versioning](http://semver.org/).

## [Unreleased]
- added new option to assign a reference PID to a geometry by MimmoGeometry execution
- added new option to write pidded 3D surfaces in MultiSolid ASCII (MimmoGeometry)
- added mimmo native dump/restore format *.geomimmo (MimmoGeometry)
- fixed minor bugs in FFDLattice manipulator
- added new features to OBBox utility (multi geometry in input, option to switch from Oriented to Axis Aligned box calculation)
- fixed coverity issues
- update CreateSeedOnSurface block: added sensitivity map option to drive the seeding. 
- added global expert mode MIMMO_EXPERT to override mandatory ports checking in execution of chains
- introducing new MimmoPiercedVector container for data fields attached to a MimmoObject.
- Added new ports for MimmoPiercedVector / Dismissed some ports no more useful.
- Dismissed geohandlers:Overlap classes, wrapped in new geohandlers:Reconstruct classes
- Added new geohandlers:Switch classes
- Dismissed geohandlers:Split classes, substituted by new geohandlers: Extract classes
- Dismissed iogeneric:MultipleMimmoGeometries
- Dismissed manipulators:MultiApply


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

