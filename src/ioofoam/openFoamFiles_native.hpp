/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#ifndef FOAM_FILES_NATIVE_H
#define FOAM_FILES_NATIVE_H

#include "MimmoPiercedVector.hpp"
#include <bitpit_patchkernel.hpp>
#include <array>
#include <vector>
#include <fvCFD.H>

namespace mimmo{

/*!
 *
 * \brief Utilities employing native OpenFOAM libraries to read an OpenFOAM mesh
 * and relative fields defined onto it.
 * \ingroup ioofoam
 */
namespace foamUtilsNative{
    int countFieldComponents(const char *rootPath, const char *fileName);

    int getFieldSize(const char *rootPath, const char *fileName, int patchIdx);
    void readScalarField(const char *rootPath, const char *fileName, int patchIdx, std::size_t &size, std::vector<double>&field);
    void writeScalarField(const char *rootPath, const char *fileName, int patchIdx, std::vector<double> &fieldd);
    void readVectorField(const char *rootPath, const char *fileName, int patchIdx,  std::size_t &size, std::vector<std::array<double,3> >&field);
    void writeVectorField(const char *rootPath, const char *fileName, int patchIdx, std::vector<std::array<double,3> > &field);

    std::string getPatchName(const char *rootPath, int patchIdx);
    int countPatches(const char *rootPath);
    int getPatchIndex(const char *rootPath, const std::string & patchName);

    void initializeCase(const char *rootPath, Foam::Time **runTime, Foam::fvMesh **mesh, int processorCount = 1);
    const word getFieldClass(const char *rootPath, const char *fileName);
    bool writePointsOnCase(const char *rootPath, std::vector<std::array<double,3> > &points, bool overwriteStart = false);

    livector1D mapEleVConnectivity(const livector1D &, const bitpit::ElementType &);

    bool interpolateFaceToPoint(dmpvector1D & facefield, dmpvector1D & pointfield);
    bool interpolateFaceToPoint(dmpvecarr3E & facefield, dmpvecarr3E & pointfield);


};//end namespace foamUtilsNative

}// end namespace mimmo.
#endif
