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
 \ *-----------------------------------------------------------------------*----*/

 //=========================================================================
 //  Description:   Handling of OpenFoam Files
 //  Author:        Andrea Iob
 //=========================================================================

#include "openFoamFiles_native.hpp"
#include <string>
#include <IOobject.H>
#include <IFstream.H>

namespace mimmo{

namespace foamUtilsNative{

using namespace Foam;
/*!
 * \return number of components of a data field in file fileName associated to mesh
 * contained in rootPath. Return 3 for vectorfield, 1 for scalarfield, -1 for unknown or
 * not supported data field.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 */
int countFieldComponents(const char *rootPath, const char *fileName)
{

    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volScalarField::typeName) {
        return 1;
    } else if (headerClassName == volVectorField::typeName) {
        return 3;
    } else {
        return -1;
    }
}

/*!
 * \return dimension of a target field associated to an openFoam mesh
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 * \param[in] patchIdx flag, if < 0 identify your field as an internal field. Otherwise, it identifies
 *                     the the field as relative to the pathIdx-th boundary patch.
 */
int getFieldSize(const char *rootPath, const char *fileName, int patchIdx){
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Read field size
    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volScalarField::typeName) {
        volScalarField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            return foamField.internalField().size();
        } else {
            return foamField.boundaryField()[patchIdx].size();
        }
    } else if (headerClassName == volVectorField::typeName) {
        volVectorField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            return foamField.internalField().size();
        } else {
            return foamField.boundaryField()[patchIdx].size();
        }
    }

    return -1;
}


/*!
 * Read a scalar field from an OpenFOAM mesh.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 * \param[in] patchIdx flag, if < 0 identify your field as an internal field. Otherwise, it identifies
 *                     the the field as relative to the pathIdx-th boundary patch.
 * \param[out] size    size of the written field
 * \param[out] field   data saved in a std::vector structure.
 */
void readScalarField(const char *rootPath, const char *fileName, int patchIdx, std::size_t &size, std::vector<double>&field){
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Read field
    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volScalarField::typeName) {
        volScalarField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            size = std::size_t(foamField.internalField().size());
            field.resize(size);
            forAll(foamField.internalField(), n) {
                field[n] = foamField[n];
            }
        } else {
            size = std::size_t(foamField.boundaryField()[patchIdx].size());
            field.resize(size);
            forAll(foamField.boundaryField()[patchIdx], n) {
                field[n] = foamField.boundaryField()[patchIdx][n];
            }
        }
    }
}

/*!
 * Read a vector field from an OpenFOAM mesh.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 * \param[in] patchIdx flag, if < 0 identify your field as an internal field. Otherwise, it identifies
 *                     the the field as relative to the pathIdx-th boundary patch.
 * \param[out] size    size of the written field
 * \param[out] field   data saved in a std::vector structure.
 */
void readVectorField(const char *rootPath, const char *fileName, int patchIdx, std::size_t &size, std::vector<std::array<double,3> >&field)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Read field
    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volVectorField::typeName) {
        volVectorField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            size = std::size_t(foamField.internalField().size());
            field.resize(size);
            forAll(foamField.internalField(), n) {
                for (int k = 0; k < 3; k++) {
                    field[n][k] = foamField[n][k];
                }
            }
        } else {
            size = std::size_t(foamField.boundaryField()[patchIdx].size());
            field.resize(size);
            forAll(foamField.boundaryField()[patchIdx], n) {
                for (int k = 0; k < 3; k++) {
                    field[n][k] = foamField.boundaryField()[patchIdx][n][k];
                }
            }
        }
    }

}

/*!
 *Overwrite/Update scalar field data on a pre-existent OpenFOAM mesh scalar field file.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 * \param[in] patchIdx flag, if < 0 identify your field as an internal field. Otherwise, it identifies
 *                     the the field as relative to the pathIdx-th boundary patch.
 * \param[out] field   data saved in a std::vector structure.
 */
void writeScalarField(const char *rootPath, const char *fileName, int patchIdx, std::vector<double> &field)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Write field
    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volScalarField::typeName) {

        volScalarField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            field.resize(foamField.internalField().size(),0.0);
            forAll(foamField.internalField(), n) {
                foamField[n] = field[n];
            }
        } else {
            field.resize(foamField.boundaryField()[patchIdx].size(),0.0);
            forAll(foamField.boundaryField()[patchIdx], n) {
#if OPENFOAM_OLDVER
                foamField.boundaryField()[patchIdx][n] = field[n];
#else
                foamField.boundaryFieldRef()[patchIdx][n] = field[n];
#endif
            }
        }

        foamField.write();
    }
}

/*!
 * Overwrite/Update vector field data on a pre-existent OpenFOAM mesh vector field file.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 * \param[in] patchIdx flag, if < 0 identify your field as an internal field. Otherwise, it identifies
 *                     the the field as relative to the pathIdx-th boundary patch.
 * \param[out] field   data saved in a std::vector structure.
 */
void writeVectorField(const char *rootPath, const char *fileName, int patchIdx,  std::vector<std::array<double,3> > & field)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Write field
    const word& headerClassName = getFieldClass(rootPath, fileName);
    if (headerClassName == volVectorField::typeName) {
        volVectorField foamField
        (
            IOobject
            (
                fileName,
             foamRunTime->timeName(),
             *foamMesh,
             IOobject::NO_READ,
             IOobject::NO_WRITE
            ),
         *foamMesh
        );

        if (patchIdx < 0) {
            field.resize(foamField.internalField().size(),{{0.0,0.0,0.0}});
            forAll(foamField.internalField(), n) {
                for (int k = 0; k < 3; k++) {
                    foamField[n][k] = field[n][k];
                }
            }
        } else {
            field.resize(foamField.boundaryField()[patchIdx].size(),{{0.0,0.0,0.0}});
            forAll(foamField.boundaryField()[patchIdx], n) {
                for (int k = 0; k < 3; k++) {
#if OPENFOAM_OLDVER
                    foamField.boundaryField()[patchIdx][n][k] = field[n][k];
#else
                    foamField.boundaryFieldRef()[patchIdx][n][k] = field[n][k];
#endif
                }
            }
        }
        foamField.write();
    }
}

/*!
 * \return the total number of boundary patches hold by the target OpenFOAM mesh.
 * \param[in] rootPath path to openfoam mesh directory
 */
int countPatches(const char *rootPath)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Patch count
    return foamMesh->boundaryMesh().size();
}


/*!
 * \return the name of the i-th boundary patch of the target OpenFOAM mesh.
 * If index is not found return an empty char string.
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] patchIdx index of the boundary patch.
 */
std::string getPatchName(const char *rootPath, int patchIdx)
{

    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Patch name
    if (patchIdx > foamMesh->boundaryMesh().size() || patchIdx < 0) {
         return std::string("");
    }

    return std::string(foamMesh->boundaryMesh()[patchIdx].name().c_str());

}


/*!
 * \return the index associated to a name of a boundary patch of the target OpenFOAM mesh.
 * If patch name is not found return a negative index -1.
 * \param[in] rootPath path to openfoam mesh directory.
 * \param[in] patchName name of the boundary patch.
 */
int getPatchIndex(const char *rootPath, const std::string &patchName)
{

    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Patch index
    forAll(foamMesh->boundaryMesh(), n) {
        if (strcmp(patchName.c_str(),foamMesh->boundaryMesh()[n].name().c_str()) == 0) {
            return n;
        }
    }

    return -1;

}

/*!
 * Read the mesh and initialize the OpenFOAM case.
 * \param[in] rootPath path to openfoam mesh directory.
 * \param[in] foamRunTime_retPtr pointer to foam Time object.
 * \param[in] foamMesh_retPtr pointer to foam fvMesh object.
 */
void initializeCase(const char *rootPath, Foam::Time **foamRunTime_retPtr, Foam::fvMesh **foamMesh_retPtr)
{
    // Arguments and root case
    int nArguments = 3;
    char **arguments = new char*[nArguments];

    char executable[] = "foamReader";
    arguments[0] = executable;

    char caseOption[] = "-case";
    arguments[1] = caseOption;

    char caseRoot[strlen(rootPath) + 1];
    strcpy(caseRoot, rootPath);
    arguments[2] = caseRoot;


    static Foam::argList *foamArgs = 0;
    if (!foamArgs) {
        foamArgs = new Foam::argList(nArguments, arguments);
        if (!foamArgs->checkRootCase()) {
            Foam::FatalError.exit();
        }
    }

    // Time
    static Foam::Time *foamRunTime = 0;
    if (!foamRunTime) {
         foamRunTime = new Foam::Time(Foam::Time::controlDictName, *foamArgs);
    }

    *foamRunTime_retPtr = foamRunTime;

    // Mesh
    static Foam::fvMesh *foamMesh = 0;
    if (!foamMesh) {
        foamMesh = new Foam::fvMesh(
            Foam::IOobject(
                Foam::fvMesh::defaultRegion,
                foamRunTime->timeName(),
                *foamRunTime,
                Foam::IOobject::MUST_READ
            )
        );
    }

    *foamMesh_retPtr = foamMesh;

}

/*!
 * \return the type of data associated to a field
 * \param[in] rootPath path to openfoam mesh directory
 * \param[in] fileName name of the file where the field is contained
 */
const word getFieldClass(const char *rootPath, const char *fileName)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    // Read field
    IOobject fieldObject
    (
        fileName,
        foamRunTime->timeName(),
        *foamMesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    Istream* headerStream = new IFstream(fieldObject.objectPath());
    fieldObject.readHeader(*headerStream);

    return fieldObject.headerClassName();

}


/*!
 * Write a pointField to an OpenFoam pre-existent case. The method use the Foam::fvMesh method movePoints
 * to update position of mesh points in the current case. Connectivity of the mesh elements remains the same.
 * Dimension of pointField argument is checked out. If it does not correspond to actual size of points in the
 * target case do nothing and return false;
 * \param[in] rootPath path to openfoam mesh directory.
 * \param[in] points new pointField to substitute
 * \param[in] overwriteStart if true new points set will be overwritten in the current case time, otherwise
 *                           a new time case incremented by one w.r.t the current one will be created.
 */
bool writePointsOnCase(const char *rootPath, std::vector<std::array<double,3> > & points, bool overwriteStart)
{
    // Initialize case
    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    initializeCase(rootPath, &foamRunTime, &foamMesh);

    functionObjectList & objfunctions = foamRunTime->functionObjects();
    objfunctions.off();

    if(std::size_t(foamMesh->points().size()) != points.size()){
        return false;
    }

    Foam::pointField movedPoints(foamMesh->points());

    forAll(movedPoints, i){
        for(Foam::label j=0; j<3; ++j){
            movedPoints[i][j] = points[i][j];
        }
    }

    foamMesh->movePoints(movedPoints);

    if(overwriteStart){
#if OPENFOAM_OLDVER
        foamMesh->writeObjects(Foam::IOstream::streamFormat::BINARY, foamRunTime->writeVersion(), foamRunTime->writeCompression());
#else
        foamMesh->writeObject(Foam::IOstream::streamFormat::BINARY, foamRunTime->writeVersion(), foamRunTime->writeCompression(), true);
#endif
    }else{

        Foam::dimensionedScalar fakeEndTime
        (
            "fakeEndTime",
         dimensionSet(0,0,1,0,0,0,0),
         scalar(1)
        );
        foamRunTime->setEndTime(foamRunTime->startTime()+fakeEndTime);

        while(foamRunTime->loop()){
#if OPENFOAM_OLDVER
            foamMesh->writeObjects(Foam::IOstream::streamFormat::BINARY, foamRunTime->writeVersion(), foamRunTime->writeCompression());
#else
            foamMesh->writeObject(Foam::IOstream::streamFormat::BINARY, foamRunTime->writeVersion(), foamRunTime->writeCompression(), true);
#endif
        }

    }
    return true;
}


/*!
 * Reorder OpenFoam vertex-cell connectivity of an elementary shape in
 * a suitable one for corrispondent bitpit elementary shape. The reordering
 * preserve the local enumeration of faces from OpenFOAM to bitpit.
 * \param[in] FoamConn ordered vertex connectivity of a shape cell
 * \param[in] eltype   reference Bitpit::element type for reordering
 */
livector1D
mapEleVConnectivity(const livector1D & FoamConn, const bitpit::ElementType & eltype){
	livector1D conn;

	switch(eltype){
	case bitpit::ElementType::TETRA:
		conn.resize(4);
		conn[0] = FoamConn[2];
		conn[1] = FoamConn[1];
		conn[2] = FoamConn[3];
		conn[3] = FoamConn[0];
		break;
	case bitpit::ElementType::HEXAHEDRON:
		conn.resize(8);
		conn[0] = FoamConn[0];
		conn[1] = FoamConn[3];
		conn[2] = FoamConn[7];
		conn[3] = FoamConn[4];
		conn[4] = FoamConn[1];
		conn[5] = FoamConn[2];
		conn[6] = FoamConn[6];
		conn[7] = FoamConn[5];
		break;
	case bitpit::ElementType::WEDGE:
		conn.resize(6);
		conn[0] = FoamConn[0];
		conn[1] = FoamConn[2];
		conn[2] = FoamConn[1];
		conn[3] = FoamConn[3];
		conn[4] = FoamConn[5];
		conn[5] = FoamConn[4];
		break;
	case bitpit::ElementType::PYRAMID:
		conn.resize(5);
		conn[0] = FoamConn[3];
		conn[1] = FoamConn[2];
		conn[2] = FoamConn[1];
		conn[3] = FoamConn[0];
		conn[4] = FoamConn[4];
		break;
	default:
		//do nothing
		break;
	}

	return conn;
}

/*!
 * It interpolates a surface scalar field defined on cell centers of a surface mesh on the vertices of the same mesh.
 * It uses the angles related to each vertex to define weight for each face value.
 * \return false if errors occurred during the interpolation.
 */
bool interpolateFaceToPoint(dmpvector1D & facefield, dmpvector1D & pointfield){

	if (facefield.getGeometry()->getType() != 1){
		throw std::runtime_error ("Error interpolating data on surface mesh field: geometry is not a surface mesh");

	}

	dmpvector1D sumAngles;
	pointfield.clear();
	for (bitpit::Vertex & vertex : facefield.getGeometry()->getVertices()){
		long idV = vertex.getId();
		sumAngles.insert(idV, 0.);
		pointfield.insert(idV, 0.);
	}

	for (bitpit::Cell & cell : facefield.getGeometry()->getCells()){

		//Cell id
		long id = cell.getId();

		// Cell data
		double cellData = facefield[id];

		// Get vertex id
		bitpit::ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

		for (int vertexLocalId=0; vertexLocalId<cell.getVertexCount(); vertexLocalId++){

			// Cell angle at vertex
			double cellVertexAngle = static_cast<bitpit::SurfaceKernel*>(facefield.getGeometry()->getPatch())->evalAngleAtVertex(id, vertexLocalId);

			//Update data vertex and sum angles
			pointfield[cellVertexIds[vertexLocalId]] += cellData * cellVertexAngle;
			sumAngles[cellVertexIds[vertexLocalId]] += cellVertexAngle;

		}

	}

	//Normalize data with sum angles
	for (bitpit::Vertex & vertex : facefield.getGeometry()->getVertices()){
		long idV = vertex.getId();
		pointfield[idV] /= sumAngles[idV];
	}

	//Assign geometry and data location to pointfield
	pointfield.setGeometry(facefield.getGeometry());
	pointfield.setDataLocation(1);

	return true;

}



/*!
 * It interpolates a vector surface field defined on cell centers of a surface mesh on the vertices of the same mesh.
 * It uses the angles related to each vertex to define weight for each face value.
 * \return false if errors occurred during the interpolation.
 */
bool interpolateFaceToPoint(dmpvecarr3E & facefield, dmpvecarr3E & pointfield){

	if (facefield.getGeometry()->getType() != 1){
		throw std::runtime_error ("Error interpolating data on surface mesh field: geometry is not a surface mesh");

	}

	dmpvecarr3E sumAngles;
	pointfield.clear();
	for (bitpit::Vertex & vertex : facefield.getGeometry()->getVertices()){
		long idV = vertex.getId();
		sumAngles.insert(idV, {{0.,0.,0.}});
		pointfield.insert(idV, {{0.,0.,0.}});
	}

	for (bitpit::Cell & cell : facefield.getGeometry()->getCells()){

		//Cell id
		long id = cell.getId();

		// Cell data
		std::array<double,3> cellData = facefield[id];

		// Get vertex id
		bitpit::ConstProxyVector<long> cellVertexIds = cell.getVertexIds();

		for (int vertexLocalId=0; vertexLocalId<cell.getVertexCount(); vertexLocalId++){

			// Cell angle at vertex
			double cellVertexAngle = static_cast<bitpit::SurfaceKernel*>(facefield.getGeometry()->getPatch())->evalAngleAtVertex(id, vertexLocalId);

			//Update data vertex and sum angles
			pointfield[cellVertexIds[vertexLocalId]] += cellData * cellVertexAngle;
			sumAngles[cellVertexIds[vertexLocalId]] += cellVertexAngle;

		}

	}

	//Normalize data with sum angles
	for (bitpit::Vertex & vertex : facefield.getGeometry()->getVertices()){
		long idV = vertex.getId();
		pointfield [idV] /= sumAngles[idV];
	}

	//Assign geometry and data location to pointfield
	pointfield.setGeometry(facefield.getGeometry());
	pointfield.setDataLocation(1);

	return true;

}



} //end of namespace foamUtilsNative

} //end of namespace mimmo
