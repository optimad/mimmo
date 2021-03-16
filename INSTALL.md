# mimmo installation

mimmo currently runs on Linux platforms. Compatibility on Windows is enabled via MSYS2/MinGW and relative instructions can be found in README_WINDOWS.md.
Full compliance with MacOS systems is currently under investigation by developers.

## Dependencies
mimmo depends on
* c++ compiler supporting `-std=c++11`. It has been tested with g++ = 4.8.5
* cmake. Tested with cmake >= 3.13.2
* bitpit library. It has been tested with bitpit 1.7.1. Visit www.optimad.it/products/bitpit/ for further information.
* (optionally) cgns. It has been tested with cgns = 3.3.1, with sub-dependency hdf5 = 1.10.4.
* (optionally) OpenFoam. It has been tested with OpenFoam Foundation 7 and ESI-OpenCFD v1906+ versions.
   Retro-compatibility with older versions has be not been verified yet and may be not guaranteed.
   Anyway, share with us your own experience, will be glad to integrate this current guide.
* (optionally) MPI implementation. It has been tested with OpenMPI >= 4.0.0.
* (optionally) Metis/parMetis. Needed by MPI version, it has been tested with Metis = 5.1.0, parMetis = 4.0.3         


## Configuring mimmo
mimmo uses ccmake as building tool.
In the mimmo's root folder make a building folder, e.g. build
```bash
    mimmo$ mkdir build
```
Enter the `build` folder
```bash
    mimmo$ cd build
```
 In order to configure it, run:
```bash
    mimmo/build$ ccmake ../
```

 By this way, mimmo can be configured for production and installation.
Setting some variable in ccmake interface you can customize a bit your configuration.

The `CMAKE_BUILD_TYPE` variable has to be used to set the type of build. The possible options are : `None`, the environment compiler flags are used; `Release`, using compiler optimization flag `-O2`; `Debug`, related to compiler flags `-O0 -fmessage-length=0`, `RelWithDebInfo`, that uses compilation flags `-O2 -g` and `MinSizeRel` to have the smallest binary size.

In addition the `ENABLE_PROFILING` variable can be set to `ON` in order to add profiling flag `-pg` during the compilation.

The `ENABLE_MPI` variable can be used to compile the parallel implementation of the mimmo packages and to allow the dependency on MPI libraries.

The `BUILD_EXAMPLES` can be used to compile examples sources in `mimmo/examples`. Note that the tests sources in `mimmo/test`are necessarily compiled and successively available at `mimmo/build/test/` as well as the compiled examples are available at `mimmo/build/examples/`.

The module variables  can be used to compile each module singularly by setting the related varible `ON/OFF`. Some modules are always compiled (as for core, manipulators), while for `MIMMO_MODULE_GEOHANDLERS`, `MIMMO_MODULE_IOCGNS`, `MIMMO_MODULE_IOOFOAM`, `MIMMO_MODULE_PROPAGATORS` and `MIMMO_MODULE_UTILS` the compilation can be toggled. Possible dependencies between mimmo modules are automatically resolved.
When possible, dependencies on external libraries are automatically resolved. Otherwise cmake will ask to specify the installation info of the missing packages.
In particular:
1) `MIMMO_MODULE_IOCGNS` will require cgns libraries and will expose a variable `CGNS_DIR` to specify manually the path to to cgns installation in case of missing or non-compliant package. If the package is regularly found, a `CGNS_DIR_FOUND` variable will be filled accordingly.
2) `MPI versions` (ENABLE_MPI on) will require `MIMMO_MODULE_PARALLEL` to be active and Metis and parMetis libraries to be available on your system. Two variables `METIS_DIR` and `PARMETIS_DIR` will be exposed  to specify manually the path to to the two installations, in case of missing or non-compliant packages. If packages are regularly found, a `METIS_DIR_FOUND` and `PARMETIS_DIR_FOUND` variables will be filled accordingly.
3) `MIMMO_MODULE_OPENFOAM` will require one of OpenFOAM distributions of OpenFOAM Foundation or ESI-OpenCFD  regularly installed on your system (and environment variables loaded too), with devel libraries and include directories available for the User. An internal search function will expose a variable `OPENFOAM_DISTRO` just to let the User specify its distribution type, if an OpenFoam Foundation or ESIOpenCFD one. All other needed information are automatically retrieved. If the package is regularly found a `OPENFOAM_DIR` and `OPENFOAM_API` variable will be filled accordingly.  


The `BUILD_XMLTUI` variable defines if the XML mimmo interpreter has to be compiled. The compiled executable `mimmo++` is available at `mimmo/build/binaries/`.

Finally, you can choose the installation folder setting the cmake variable `CMAKE_INSTALL_PREFIX`. The default installation folder is `/usr/local/`.

Remember that if you choose the default installation path or another path without write permission you will need administration privileges to install mimmo in.

## Building and Installing
Once cmake has configured mimmo's building, just do
```bash
    mimmo/build$ make   
```
to build and
```bash
    mimmo/build$ make install   
```
to install.

If you have just built mimmo, its headers will be available at `mimmo/include/` folder and a static library `libmimmo.a` will be available at `mimmo/build/lib/` folder.

<!--- (or `libmimmo_MPI.a` in case of parallel compilation) -->

If you have also installed mimmo, its headers will be available at `/my/installation/folder/mimmo/include/` folder and a static library `libmimmo.a` will be available at `/my/installation/folder/lib/` folder. The `mimmo++` executable will be installed in `/my/installation/folder/bin/`.

A shared version of the library is provided setting the cmake variable `BUILD_SHARED_LIBS` to ON, during the ccmake settings.

You can test the compilation by run the command `ctest` from the building folder.

For a complete guide to installation of mimmo please visit
<a href="http://www.optimad.github.io/mimmo/documentation/installation.html">mimmo installation webpage</a>.

## Building Documentation
In order to build properly the documentation Doxygen and Graphviz are needed.
Doxygen and Graphviz versions currently employed to test documentation are 1.8.16 and 2.20.2 respectively.
In the ccmake interface the variable `BUILD_DOCUMENTATION` can be set to `ON` in order to build the documentation during the library compilation.

<!-- If turned on the new variable `DOC_EXTRACT_PRIVATE` can be used to include all the private class members in the documentation. -->

After the `make` or `make install` the doxygen documentation will be built. You can chose to compile only the documentation with command
```bash
    mimmo/build$ make doc   
```
You can now browse the html documentation with your favorite browser by opening 'html/index.html'.

## Help
For any problem, please contact <a href="http://www.optimad.it">Optimad engineering srl</a> at info@optimad.it.
