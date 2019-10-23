#!/bin/bash

MIMMO_SOURCE_DIR="$1"
BITPIT_DIR="$2"
LOCAL_BUILD="$3"

MANUAL_ROOT_DIR="${PWD}"
JEKYLL_ROOT_DIR="${MANUAL_ROOT_DIR}/.."

MANUAL_BUILD_DIR="${MANUAL_ROOT_DIR}/build"

# Generate build directory
rm -rf "${MANUAL_BUILD_DIR}"
mkdir "${MANUAL_BUILD_DIR}"
cd "${MANUAL_BUILD_DIR}"

# Link files needed by jekyll
cp -r "${MANUAL_ROOT_DIR}/_includes" "${MANUAL_BUILD_DIR}"
cp -r "${MANUAL_ROOT_DIR}/_layouts" "${MANUAL_BUILD_DIR}"
cp -r "${MANUAL_ROOT_DIR}/stylesheets" "${MANUAL_BUILD_DIR}"
cp "${JEKYLL_ROOT_DIR}/_includes/"* "${MANUAL_BUILD_DIR}/_includes/"
cp "${JEKYLL_ROOT_DIR}/_layouts/"* "${MANUAL_BUILD_DIR}/_layouts/"
cp "${MANUAL_ROOT_DIR}/doxygen_footer.html" "${MANUAL_BUILD_DIR}"
cp "${MANUAL_ROOT_DIR}/doxygen_header.html" "${MANUAL_BUILD_DIR}"

# Configration files needed by jekyll
JEKYLL_CONFIG_GLOBAL="${JEKYLL_ROOT_DIR}/_config.yml"
JEKYLL_CONFIG_LOCAL="${JEKYLL_ROOT_DIR}/_config.local.yml"
JEKYLL_CONFIG_MANUAL="${MANUAL_BUILD_DIR}/_config.manual.yml"

JEKYLL_CONFIG_MANUAL_CONTENTS="""
layouts_dir:  ${MANUAL_BUILD_DIR}/_layouts

exclude: [generate.sh]

"""
echo "${JEKYLL_CONFIG_MANUAL_CONTENTS}" > "${JEKYLL_CONFIG_MANUAL}"

JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_GLOBAL}"
if [ -f "${JEKYLL_CONFIG_LOCAL}" -a ${LOCAL_BUILD} -eq 1 ]; then
    JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_LIST},${JEKYLL_CONFIG_LOCAL}"
fi
JEKYLL_CONFIG_LIST="${JEKYLL_CONFIG_LIST},${JEKYLL_CONFIG_MANUAL}"

# Generate templates needed by doxygen
cd ${JEKYLL_ROOT_DIR}

bundle exec jekyll build \
    --config "${JEKYLL_CONFIG_LIST}" \
    --source ${MANUAL_BUILD_DIR} \
    --destination ${MANUAL_BUILD_DIR}/templates

# Generate doxygen files
cd ${MANUAL_BUILD_DIR}
mkdir mimmo
cd mimmo

cmake ${MIMMO_SOURCE_DIR} \
    -DCMAKE_INSTALL_PREFIX:PATH="" \
    -DBITPIT_DIR="${BITPIT_DIR}" \
    -DBUILD_DOCUMENTATION=ON \
    -DENABLE_MPI=OFF \
    -DMIMMO_MODULE_GEOHANDLERS=ON \
    -DMIMMO_MODULE_IOCGNS=ON \
    -DMIMMO_MODULE_IOOFOAM=ON \
    -DMIMMO_MODULE_IOVTK=ON \
    -DMIMMO_MODULE_PROPAGATORS=ON \
    -DMIMMO_MODULE_PARALLEL=OFF \
    -DMIMMO_MODULE_UTILS=ON \
    -DDOXY_HTML_HEADER="${MANUAL_BUILD_DIR}/templates/doxygen_header.html" \
    -DDOXY_HTML_FOOTER="${MANUAL_BUILD_DIR}/templates/doxygen_footer.html" \
    -DDOXY_HTML_EXTRA_STYLESHEET="${MANUAL_BUILD_DIR}/stylesheets/manual.css" \
    -DDOXY_DISABLE_INDEX="NO" \
    -DDOXY_GENERATE_TREEVIEW="NO" \
    -DPETSC_COMPILER:FILEPATH="mpicxx" \
    -DPETSC_CURRENT:BOOL=OFF \
    -DPETSC_DEFINITIONS:STRING="-D__INSDIR__=" \
    -DPETSC_DIR:PATH="/usr/lib/petsc" \
    -DPETSC_INCLUDES:STRING="/usr/lib/petscdir/3.7.7/x86_64-linux-gnu-real/include;/usr/include/hypre;/usr/include/superlu;/usr/include/suitesparse;/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include;/usr/include/scotch;/usr/lib/x86_64-linux-gnu/openmpi/include" \
    -DPETSC_INCLUDE_CONF:PATH="/usr/lib/petsc/include" \
    -DPETSC_INCLUDE_DIR:PATH="/usr/lib/petsc/include" \
    -DPETSC_LIBRARIES:STRING="/usr/lib/petscdir/3.7.7/x86_64-linux-gnu-real/lib/libpetsc_real.so;/usr/lib/x86_64-linux-gnu/libdmumps.so;/usr/lib/x86_64-linux-gnu/libzmumps.so;/usr/lib/x86_64-linux-gnu/libsmumps.so;/usr/lib/x86_64-linux-gnu/libcmumps.so;/usr/lib/x86_64-linux-gnu/libmumps_common.so;/usr/lib/x86_64-linux-gnu/libpord.so;/usr/lib/x86_64-linux-gnu/libHYPRE_IJ_mv.so;/usr/lib/x86_64-linux-gnu/libHYPRE_parcsr_ls.so;/usr/lib/x86_64-linux-gnu/libHYPRE_sstruct_ls.so;/usr/lib/x86_64-linux-gnu/libHYPRE_sstruct_mv.so;/usr/lib/x86_64-linux-gnu/libHYPRE_struct_ls.so;/usr/lib/x86_64-linux-gnu/libHYPRE_struct_mv.so;/usr/lib/x86_64-linux-gnu/libHYPRE_utilities.so;/usr/lib/x86_64-linux-gnu/libsuperlu.so;/usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so;/usr/lib/x86_64-linux-gnu/libfftw3.so;/usr/lib/x86_64-linux-gnu/libfftw3_mpi.so;/usr/lib/x86_64-linux-gnu/libumfpack.so;/usr/lib/x86_64-linux-gnu/libamd.so;/usr/lib/x86_64-linux-gnu/libcholmod.so;/usr/lib/x86_64-linux-gnu/libklu.so;/usr/lib/x86_64-linux-gnu/liblapack.so;/usr/lib/x86_64-linux-gnu/libblas.so;/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/libhdf5hl_fortran.so;/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/libhdf5_fortran.so;/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/libhdf5_hl.so;/usr/lib/x86_64-linux-gnu/hdf5/openmpi/lib/libhdf5.so;/usr/lib/x86_64-linux-gnu/libhwloc.so;/usr/lib/x86_64-linux-gnu/libptesmumps.so;/usr/lib/x86_64-linux-gnu/libptscotch.so;/usr/lib/x86_64-linux-gnu/libptscotcherr.so;/usr/lib/x86_64-linux-gnu/libssl.so;/usr/lib/x86_64-linux-gnu/libcrypto.so;/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempif08.so;/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_usempi_ignore_tkr.so;/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_mpifh.so;/usr/lib/gcc/x86_64-linux-gnu/7/libgfortran.so;/usr/lib/gcc/x86_64-linux-gnu/7/libquadmath.so;/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so;/usr/lib/gcc/x86_64-linux-gnu/7/libstdc++.so;/usr/lib/x86_64-linux-gnu/libm.so;/usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so;/usr/lib/gcc/x86_64-linux-gnu/7/libgcc_s.so;/usr/lib/x86_64-linux-gnu/libpthread.so;/usr/lib/x86_64-linux-gnu/libdl.so"
    -DPETSC_LIBRARY_SINGLE:FILEPATH="/usr/lib/petscdir/3.7.7/x86_64-linux-gnu-real/lib/libpetsc_real.so" \
    -DPETSC_MPIEXEC:FILEPATH="mpiexec" \

    
cd doc

make

make DESTDIR=${MANUAL_BUILD_DIR} install/local

# Copy the generated manual into the website
cd ${MANUAL_BUILD_DIR}

MIMMO_VERSION=`ls doc`
MIMMO_VERSION=${MIMMO_VERSION#"mimmo-"}
MIMMO_VERSION="${MIMMO_VERSION}/"

MANUAL_INSTALL_DIR="${JEKYLL_ROOT_DIR}/documentation/manual/${MIMMO_VERSION}"

rm -rf "${MANUAL_INSTALL_DIR}"
mkdir -p "${MANUAL_INSTALL_DIR}"
cp -r "doc/mimmo-${MIMMO_VERSION}/html/"* "${MANUAL_INSTALL_DIR}"
