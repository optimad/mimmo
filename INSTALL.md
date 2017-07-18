# mimmo installation

mimmo currently runs on Linux platforms. Compatibility on Windows and MacOSX is currently being tested by developers.

## Dependencies
mimmo depends on
* c++ compiler supporting `-std=c++11`. It has been tested with g++ >= 4.7.3
* cmake >= 2.8
* lapacke/lapack libraries. It has been tested with Lapack >= 3.5.0
* xml2 libraries. (should be provided by default on Linux system). Tested with LibXml2 >= 2.9.1
* bitpit library. It has been tested with bitpit >= 1.4.0. Visit www.optimad.it/products/bitpit/ for further information.
* (optionally) MPI implementation. It has been tested with OpenMPI >= 1.6.5. 
* (optionally) vtk. It has been tested with vtk >= 6.3. 
* (optionally) cgns. It has been tested with cgns = 3.2.1. 
* (optionally) hdf5. It has been tested with hdf5 = 1.8.15. 

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

The module variables (available in the advanced mode) can be used to compile each module singularly by setting the related varible `ON/OFF`. `MIMMO_MODULE_CORE` is always compiled, while for `MIMMO_MODULE_GEOHANDLERS`, `MIMMO_MODULE_IOCGNS`, `MIMMO_MODULE_IOOFOAM`, `MIMMO_MODULE_IOVTK` and `MIMMO_MODULE_UTILS` the compilation can be toggled. Possible dependencies between mimmo modules are automatically resolved.
Dependencies on external libraries when possible are automatically resolved  through find package command.

The `BUILD_XMLTUI` variable defines if the executable binary has to be compiled. The compiled executable `mimmo++` is available at `mimmo/build/binaries/`.

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

If you have just built mimmo, its headers will be available at `mimmo/include/` folder and a static library `libmimmo.a` (or `libmimmo_MPI.a` in case of parallel compilation) will be available at `mimmo/build/lib/` folder.

If you have also installed mimmo, its headers will be available at `/my/installation/folder/mimmo/include/` folder and a static library `libmimmo.a` will be available at `/my/installation/folder/lib/` folder. The `mimmo++` executable will be installed in `/my/installation/folder/bin/`.

A shared version of the library is provided setting the cmake variable `BUILD_SHARED_LIBS` to ON, during the ccmake settings.

You can test the compilation by run the command `ctest` from the building folder.

## Building Documentation
In order to build properly the documentation Doxygen (>=1.8.6) and Graphviz (>=2.20.2) are needed.

In the ccmake interface the variable `BUILD_DOCUMENTATION` can be set to `ON` in order to build the documentation during the library compilation. 
If turned on the new variable `DOC_EXTRACT_PRIVATE` can be used to include all the private class members in the documentation.
  
After the `make` or `make install` the doxygen documentation will be built. You can chose to compile only the documentation with command 
```bash
    mimmo/build$ make doc   
```
You can now browse the html documentation with your favorite browser by opening 'html/index.html'.

## Help
For any problem, please contact <a href="http://www.optimad.it">Optimad engineering srl</a> at info@optimad.it. 
