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

//=========================================================================
//  Description:   Handling of OpenFoam Files
//  Author:        Andrea Iob
//=========================================================================

#ifndef FOAM_FILES_NATIVE_H
#define FOAM_FILES_NATIVE_H

#include "fvCFD.H"
#include <array>
#include <vector>
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
    
    char * getPatchName(const char *rootPath, int patchIdx);
    int countPatches(const char *rootPath);
    int getPatchIndex(const char *rootPath, char *patchName);

    void initializeCase(const char *rootPath, Foam::Time **runTime, Foam::fvMesh **mesh);
    const word getFieldClass(const char *rootPath, const char *fileName);
    bool writePointsOnCase(const char *rootPath, std::vector<std::array<double,3> > &points, bool overwriteStart = false);
};//end namespace foamUtilsNative

}// end namespace mimmo.
#endif
