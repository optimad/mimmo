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
\*---------------------------------------------------------------------------*/
#include "IOOFOAM.hpp"
#include "openFoamFiles_native.hpp"

namespace mimmo{

/*!
 * Default constructor of IOOFOAM.
 * \param[in] type int from enum IOOFMode, for class mode: READ, WRITE, WRITEPOINTSONLY. Default is READ.
 */
IOOFOAM::IOOFOAM(int type):MimmoFvMesh(){
	m_name          = "mimmo.IOOFOAM";
	auto maybeIOMode = IOOFMode::_from_integral_nothrow(type);
	if(maybeIOMode){
		m_type = type;
	}else{
		m_type = IOOFMode::READ;
	}
	setDefaults();
}

/*!
 * Custom constructor. Pass bulk and boundary of the mesh.
 * The ownership of the argument will be transferred to the internal 
 * members of the class.
 * Please note, The class admits only the following combinations as bulk/boundary pair:
 *
 *  - bulk MimmoObject of type 2-Volume and boundary MimmoObject of type 1-Surface
 *  - bulk MimmoObject of type 1-Surface and boundary MimmoObject of type 4-3DCurve.
 * 
 * The class is meant for writing modes WRITE and WRITEPOINTSONLY. If a READ mode is passed 
 * in argument, it throws an error.
 * 
 *\param[in] bulk unique pointer to the bulk MimmoObject
 *\param[in] boundary unique pointer to the boundary MimmoObject
 *\param[in] type int from enum IOOFMode, for class mode: WRITE, WRITEPOINTSONLY. Default is WRITE.
 *
 */
IOOFOAM::IOOFOAM(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type): MimmoFvMesh(bulk,boundary){
	auto maybeIOMode = IOOFMode::_from_integral_nothrow(type);
	if(maybeIOMode){
		if(maybeIOMode->_to_integral() == IOOFMode::READ){
			throw std::runtime_error("Error in IOOFOAM constructor. Forcing READ mode in a constructor meant for writing modes only");
		}
		m_type = type;
	}else{
		m_type = IOOFMode::WRITE;
	}
	setDefaults();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAM::IOOFOAM(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.IOOFOAM";
	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);

	std::string fallback_type = "IOModeNONE";
	std::string input_type = rootXML.get("IOMode", fallback_type);
	input_type = bitpit::utils::string::trim(input_type);

	auto maybeIOOFMode = IOOFMode::_from_string_nothrow(input_type.c_str());

	if(input == "mimmo.IOOFOAM" && maybeIOOFMode){
		m_type = maybeIOOFMode->_to_integral();
		setDefaults();
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of IOOFOAM.
 */
IOOFOAM::~IOOFOAM(){};

/*!Copy constructor of IOOFOAM. Internal volume and boundary mesh are not copied.
 */
IOOFOAM::IOOFOAM(const IOOFOAM & other):MimmoFvMesh(other){
	m_type = other.m_type;
	m_path = other.m_path;
	m_overwrite = other.m_overwrite;
	m_OFbitpitmapfaces = other.m_OFbitpitmapfaces;
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
};

/*!
 * Assignment operator. Internal volume and boundary mesh are not copied.
 */
IOOFOAM & IOOFOAM::operator=(IOOFOAM other){
	swap(other);
	return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void IOOFOAM::swap(IOOFOAM & x) noexcept
{
	std::swap(m_type, x.m_type);
	std::swap(m_path, x.m_path);
	std::swap(m_overwrite, x.m_overwrite);
	std::swap(m_OFbitpitmapfaces, x.m_OFbitpitmapfaces);

	MimmoFvMesh::swap(x);
};

/*!
 * Default values for IOOFOAM.
 */
void
IOOFOAM::setDefaults(){

	m_path   = ".";
	m_overwrite = false;
	m_OFE_supp.clear();
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAM::buildPorts(){
	MimmoFvMesh::buildPorts();
	bool built = m_arePortsBuilt;
	// creating input ports
	built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAM>(this, &IOOFOAM::setFacesMap, M_UMAPIDS));
	// creating output ports
	built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &MimmoFvMesh::getGeometry, M_GEOMOFOAM));
	built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &MimmoFvMesh::getBoundaryGeometry, M_GEOMOFOAM2));
	built = (built && createPortOut<std::unordered_map<long,long>, IOOFOAM>(this, &IOOFOAM::getFacesMap, M_UMAPIDS));

	m_arePortsBuilt = built;
};

/*!
 * It sets the name of directory to read/write the OpenFOAM mesh.
 * \param[in] dir mesh input directory.
 */
void
IOOFOAM::setDir(std::string dir){
	m_path = dir;
}

/*!
 * Set overwrite parameter. This option is valid only in WRITEPOINTSONLY mode.
 * If true overwrite points in the current OpenFoam case time of the mesh at WriteDir. 
 * If false save them in a newly created case time at current time + 1;
 * \param[in] flag activation flag.
 */
void
IOOFOAM::setOverwrite(bool flag){
	m_overwrite = flag;
}

/*!
 * Get Overwrite parameter.  See setOverwrite method.
 * \return overwrite parameter content
 */
bool
IOOFOAM::getOverwrite(){
	return m_overwrite;
}

/*!
 * Set current bulk geometry to be written. Method meant for writing class modes only.
 * Any previous mesh internally allocated will be destroyed. See MimmoFvMesh::setGeometry
 * \param[in] bulk pointer to external bulk mesh MimmoObject.
 */
void
IOOFOAM::setGeometry(MimmoObject * bulk){
	if(m_type == IOOFMode::READ) return;
	MimmoFvMesh::setGeometry(bulk);
}

/*!
 * Set current boundary geometry to be written. Method meant for writing class modes only.
 * Any previous mesh internally allocated will be destroyed. See MimmoFvMesh::setBoundaryGeometry
 * \param[in] bulk pointer to external boundary mesh MimmoObject.
 */
void
IOOFOAM::setBoundaryGeometry(MimmoObject * boundary){
	if(m_type == IOOFMode::READ) return;
	MimmoFvMesh::setBoundaryGeometry(boundary);
}


/*!
 * Get current map between OFOAM faces -> bitpit interfaces.
 * \return reference to interfaces map.
 */
std::unordered_map<long,long>
IOOFOAM::getFacesMap(){
	return m_OFbitpitmapfaces;
}

/*!
 * Set current map between OFOAM faces -> bitpit interfaces.
 * \param[in] reference to interfaces map.
 */
void
IOOFOAM::setFacesMap(std::unordered_map<long,long> mapFaces){
	m_OFbitpitmapfaces = mapFaces;
}

/*!Execution command.
 * It reads the geometry in reading mode, otherwise it writes, according to the write mode.
 */
void
IOOFOAM::execute(){
	bool check = true;
	switch (m_type){
	case IOOFMode::READ :
		check = read();
		break;
	case IOOFMode::WRITE :
		//             check = write(); //TODO coding complete mesh writing from scratch.
		//             break;
	case IOOFMode::WRITEPOINTSONLY :
		check = writePointsOnly();
		break;
	default:
		//never been reached
		break;
	}
	if (!check){
		throw std::runtime_error (m_name + ": an error occured while reading from/writing to files");
	}
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	std::string input;

	BaseManipulation::absorbSectionXML(slotXML, name);

	if(slotXML.hasOption("Dir")){
		input = slotXML.get("Dir");
		input = bitpit::utils::string::trim(input);
		if(input.empty())   input = ".";
		setDir(input);
	};

	if(slotXML.hasOption("Overwrite")){
		input = slotXML.get("Overwrite");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
		}
		setOverwrite(value);
	};
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

	slotXML.set("IOMode", IOOFMode::_from_integral(m_type)._to_string());
	slotXML.set("Dir", m_path);
	slotXML.set("Overwrite", std::to_string(m_overwrite));
};

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
 * It reads the OpenFOAM mesh from input file and store in a the class structures m_bulk and m_boundary.
 * \return false if errors occured during the reading.
 */
bool
IOOFOAM::read(){

	Foam::Time *foamRunTime = 0;
	Foam::fvMesh *foamMesh = 0;

	//read mesh from OpenFoam case directory
	foamUtilsNative::initializeCase(m_path.c_str(), &foamRunTime, &foamMesh);

	//prepare my bulk geometry container
	std::unique_ptr<bitpit::PatchKernel> mesh(new mimmo::MimmoVolUnstructured(3));
	mesh->reserveVertices(std::size_t(foamMesh->nPoints()));
	mesh->reserveCells(std::size_t(foamMesh->nCells()));

	//start absorbing mesh nodes/points.
	Foam::pointField nodes = foamMesh->points();
	darray3E coords;
	forAll(nodes, in){
		for (int k = 0; k < 3; k++) {
			coords[k] = nodes[in][k];
		}
		mesh->addVertex(coords, long(in));
	}

	//absorbing cells.
	const Foam::cellList & cells           = foamMesh->cells();
	const Foam::cellShapeList & cellShapes = foamMesh->cellShapes();
	const Foam::faceList & faces           = foamMesh->faces();
	const Foam::labelList & faceOwner      = foamMesh->faceOwner();
	const Foam::labelList & faceNeighbour  = foamMesh->faceNeighbour();

	Foam::label sizeNeighbours = faceNeighbour.size();

	std::string eleshape;
	bitpit::ElementType eltype;
	long iDC;
	short PID = 0;
	livector1D conn, temp;
	livector2D adjacency;

	forAll(cells, iC){

		iDC = long(iC);
		eleshape = std::string(cellShapes[iC].model().name());
		//first step verify the model in twin cellShapes list.
		conn.clear();
		adjacency.clear();
		adjacency.resize(std::size_t(cells[iC].size()), livector1D(1, bitpit::Cell::NULL_ID));

		if(m_OFE_supp.count(eleshape) > 0){

			eltype = m_OFE_supp[eleshape];
			temp.resize(cellShapes[iC].size());
			forAll(cellShapes[iC], loc){
				temp[loc] = (long)cellShapes[iC][loc];
			}
			conn = mapEleVConnectivity(temp, eltype);

			auto ordFaceList = cellShapes[iC].meshFaces(faces, cells[iC]);
			Foam::label refFace;
			forAll(ordFaceList, ofcount){
				refFace = ordFaceList[ofcount];
				if(refFace >= sizeNeighbours) continue;
				if(iC == faceOwner[refFace]){
					adjacency[int(ofcount)][0] = faceNeighbour[refFace];
				}else{
					adjacency[int(ofcount)][0] = faceOwner[refFace];
				}
			}

		}else{

			eltype = bitpit::ElementType::POLYHEDRON;
			//manually calculate connectivity.
			conn.push_back((long)cells[iC].size()); //total number of faces on the top.
			forAll(cells[iC], locC){
				Foam::label iFace = cells[iC][locC];
				long faceNVertex = (long)faces[iFace].size();
				temp.resize(faceNVertex);
				forAll(faces[iFace], locF){
					temp[locF] = (long)faces[iFace][locF];
				}

				conn.push_back(faceNVertex);
				if(iFace >= sizeNeighbours) {
					//border face, normal outwards, take as it is
					conn.insert(conn.end(), temp.begin(), temp.end());
					continue;
				}
				bool normalIsOut;
				//recover right adjacency and check if face normal pointing outwards.
				//OpenFoam policy wants the face normal between cell pointing towards
				// the cell with greater id.
				if(iC == faceOwner[iFace]){
					adjacency[int(locC)][0] = faceNeighbour[iFace];
					normalIsOut = faceNeighbour[iFace] > iC;
				}else{
					adjacency[int(locC)][0] = faceOwner[iFace];
					normalIsOut = faceOwner[iFace] > iC;
				}

				if(normalIsOut){
					conn.insert(conn.end(), temp.begin(), temp.end());
				}else{
					conn.insert(conn.end(), temp.rbegin(), temp.rend());
				}
			}
		}
		bitpit::PatchKernel::CellIterator it = mesh->addCell(eltype, true, conn, iDC);
		it->setPID(int(PID));
		it->setAdjacencies(adjacency);
	}
	//build as bitpit the Interfaces -> adjacency is automatically recoverd from openFoam info.
	mesh->buildInterfaces();

	forAll(cells, iC){
		iDC = long(iC);
		eleshape = std::string(cellShapes[iC].model().name());
		Foam::labelList ofFaceList;
		if(m_OFE_supp.count(eleshape) > 0){
			eltype = m_OFE_supp[eleshape];
			ofFaceList = cellShapes[iC].meshFaces(faces, cells[iC]);
		}else{
			eltype = bitpit::ElementType::POLYHEDRON;
			ofFaceList = cells[iC];
		}
		long * bitFaceList = mesh->getCell(iC).getInterfaces();
		std::size_t sizeFList = mesh->getCell(iC).getInterfaceCount();
		for(std::size_t i = 0; i<sizeFList; ++i){
			m_OFbitpitmapfaces[ofFaceList[i]] = bitFaceList[i];
		}
	}

	//finally store bulk mesh in the internal bulk member of the class (from MimmoFvMesh);
	m_bulk = std::move(std::unique_ptr<MimmoObject>(new MimmoObject(2, mesh)));
	m_internalBulk = true;
	m_bulkext = NULL;

	//from MimmoFvMesh protected utilities, create the raw boundary mesh, storing it in m_boundary internal member.
	//PID will be every where 0. Need to be compiled to align with patch division of foamBoundaryMesh.
	//every cell of the boundary mesh will have the same id of the border interfaces of the bulk.
	createBoundaryMesh();

	//Once OFoam faces and bitpit Interfaces link is set,
	//Extract boundary patch info from foamBoundary and
	// get boundary patch division automatically pidding the class boundary mesh.

	const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
	long startIndex;
	long endIndex;

	forAll(foamBMesh, iBoundary){
		PID = short(iBoundary+1);
		startIndex = foamBMesh[iBoundary].start();
		endIndex = startIndex + long(foamBMesh[iBoundary].size());
		for(long ind=startIndex; ind<endIndex; ++ind){
			m_boundary->setPIDCell(m_OFbitpitmapfaces[ind], PID);
		}
	}
	m_boundary->resyncPID();
	//minimo sindacale fatto.
	return true;

}

/*!
 * It writes the OpenFOAM mesh to an output file from internal structure m_bulk and m_boundary
 * \return false if errors occured during the writing.
 */
bool
IOOFOAM::write(){
	if(!checkMeshCoherence()) return false;

	//do nothing for now

	return true;
}

/*!
 * It writes only the set of points coordinates in m_bulk internal mesh on a pre-existent OpenFOAM mesh
 * \return false if errors occured during the writing.
 */
bool
IOOFOAM::writePointsOnly(){
	if(!checkMeshCoherence()) return false;

	dvecarr3E points = getGeometry()->getVertexCoords();

	return foamUtilsNative::writePointsOnCase(m_path.c_str(), points, m_overwrite);
}



/*
 * ========================================================================================================
 */


/*!
 * Default constructor of IOOFOAMField.
 * \param[in] type int from enum IOOFMode, for class mode: READ, WRITE. Default is READ.
 */
IOOFOAMField::IOOFOAMField(int type):MimmoFvMesh(){
	m_name           = "mimmo.IOOFOAMField";
	auto maybeIOMode = IOOFMode::_from_integral_nothrow(type);
	if(maybeIOMode){
		m_type = type;
	}else{
		m_type = IOOFMode::READ;
	}
	setDefaults();
}


/*!Default destructor of IOOFOAMField.
 */
IOOFOAMField::~IOOFOAMField(){};

///*!Copy constructor of IOOFOAMField. Internal volume and boundary mesh are not copied.
// */
//IOOFOAMField::IOOFOAMField(const IOOFOAMField & other):MimmoFvMesh(other){
//    m_type = other.m_type;
//    m_path = other.m_path;
//    m_path = other.m_path;
//    m_fieldname = other.m_fieldname;
//    m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
//    m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
//    m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
//    m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
//};
//
///*!
// * Assignment operator. Internal volume and boundary mesh are not copied.
// */
//IOOFOAMField & IOOFOAMField::operator=(IOOFOAMField other){
//    swap(other);
//    return *this;
//}
//
///*!
// * Swap function
// * \param[in] x object to be swapped
// */
//void IOOFOAMField::swap(IOOFOAMField & x) noexcept
//{
//   std::swap(m_type, x.m_type);
//   std::swap(m_path, x.m_path);
//   std::swap(m_overwrite, x.m_overwrite);
//   std::swap(m_fieldname, x.m_fieldname);
//
//   MimmoFvMesh::swap(x);
//};

/*!
 * Default values for IOOFOAMField.
 */
void
IOOFOAMField::setDefaults(){

	m_path   = ".";
	m_fieldname = "";
	m_overwrite = false;
	m_OFE_supp.clear();
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAMField::buildPorts(){
	MimmoFvMesh::buildPorts();
	bool built = m_arePortsBuilt;
	// creating input ports
	built = (built && createPortIn<MimmoObject*, IOOFOAMField>(this, &MimmoFvMesh::setGeometry, M_GEOMOFOAM ));
	built = (built && createPortIn<MimmoObject*, IOOFOAMField>(this, &MimmoFvMesh::setBoundaryGeometry, M_GEOMOFOAM2));
	built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAMField>(this, &IOOFOAMField::setFacesMap, M_UMAPIDS));
	// creating output ports
	built = (built && createPortOut<MimmoObject*, IOOFOAMField>(this, &MimmoFvMesh::getGeometry, M_GEOMOFOAM));
	built = (built && createPortOut<MimmoObject*, IOOFOAMField>(this, &MimmoFvMesh::getBoundaryGeometry, M_GEOMOFOAM2));
	m_arePortsBuilt = built;
};

/*!
 * It sets the name of directory to read/write the OpenFOAM mesh.
 * \param[in] dir mesh input directory.
 */
void
IOOFOAMField::setDir(std::string dir){
	m_path = dir;
}

/*!
 * It sets the name of the field to read/write.
 * \param[in] fieldname name of input/output field.
 */
void
IOOFOAMField::setField(std::string fieldname){
	m_fieldname = fieldname;
}

/*!
 * Set overwrite parameter.
 * If true overwrite field in the current OpenFoam case time of the mesh at WriteDir.
 * If false save them in a newly created case time at current time + 1;
 * \param[in] flag activation flag.
 */
void
IOOFOAMField::setOverwrite(bool flag){
	m_overwrite = flag;
}

/*!
 * Get Overwrite parameter.  See setOverwrite method.
 * \return overwrite parameter content
 */
bool
IOOFOAMField::getOverwrite(){
	return m_overwrite;
}

/*!
 * Set type parameter.
 * \param[in] type parameter content
 */
void
IOOFOAMField::setType(int type){
	m_type = type;
}

/*!
 * Get type parameter.
 * \return type parameter content
 */
int
IOOFOAMField::getType(){
	return m_type;
}

/*!
 * Set current map between OFOAM faces -> bitpit interfaces.
 * \param[in] reference to interfaces map.
 */
void
IOOFOAMField::setFacesMap(std::unordered_map<long,long> mapFaces){
	m_OFbitpitmapfaces = mapFaces;
}

/*!Execution command.
 * It reads the geometry in reading mode, otherwise it writes, according to the write mode.
 */
void
IOOFOAMField::execute(){
	bool check = true;
	switch (m_type){
	case IOOFMode::READ :
		check = read();
		break;
	case IOOFMode::WRITE :
		//             check = write(); //TODO coding complete mesh writing from scratch.
		break;
	default:
		//never been reached
		break;
	}
	if (!check){
		throw std::runtime_error (m_name + ": an error occured while reading from/writing to files");
	}
}


/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAMField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	std::string input;
	BaseManipulation::absorbSectionXML(slotXML, name);

	if(slotXML.hasOption("Dir")){
		input = slotXML.get("Dir");
		input = bitpit::utils::string::trim(input);
		if(input.empty())   input = ".";
		setDir(input);
	};

	if(slotXML.hasOption("Field")){
		input = slotXML.get("Field");
		input = bitpit::utils::string::trim(input);
		if(input.empty())   input = ".";
		setField(input);
	};

	if(slotXML.hasOption("Overwrite")){
		input = slotXML.get("Overwrite");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
		}
		setOverwrite(value);
	};
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAMField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	BaseManipulation::flushSectionXML(slotXML, name);

	slotXML.set("IOMode", IOOFMode::_from_integral(m_type)._to_string());
	slotXML.set("Dir", m_path);
	slotXML.set("Field", m_fieldname);
	slotXML.set("Overwrite", std::to_string(m_overwrite));
};



/*
 * ========================================================================================================
 */

/*!
 * Default constructor of IOOFOAMScalarField.
 * \param[in] type int from enum IOOFMode, for class mode: READ, WRITE. Default is READ.
 */
IOOFOAMScalarField::IOOFOAMScalarField(int type):IOOFOAMField(){
	m_name           = "mimmo.IOOFOAMScalarField";
	IOOFOAMField::setDefaults();
}


/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAMScalarField::IOOFOAMScalarField(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.IOOFOAMScalarField";
	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);

	std::string fallback_type = "IOModeNONE";
	std::string input_type = rootXML.get("IOMode", fallback_type);
	input_type = bitpit::utils::string::trim(input_type);

	auto maybeIOOFMode = IOOFMode::_from_string_nothrow(input_type.c_str());

	if(input == "mimmo.IOOFOAMScalarField" && maybeIOOFMode){
		m_type = maybeIOOFMode->_to_integral();
		IOOFOAMField::setDefaults();
		IOOFOAMField::absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of IOOFOAMScalarField.
 */
IOOFOAMScalarField::~IOOFOAMScalarField(){};

/*!Copy constructor of IOOFOAMScalarField. Internal volume and boundary mesh are not copied.
 */
IOOFOAMScalarField::IOOFOAMScalarField(const IOOFOAMScalarField & other):IOOFOAMField(other){
	m_type = other.m_type;
	m_path = other.m_path;
	m_fieldname = other.m_fieldname;
	m_field = other.m_field;
	m_boundaryField = other.m_boundaryField;
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
};

/*!
 * Assignment operator. Internal volume and boundary mesh are not copied.
 */
IOOFOAMScalarField & IOOFOAMScalarField::operator=(IOOFOAMScalarField other){
	swap(other);
	return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void IOOFOAMScalarField::swap(IOOFOAMScalarField & x) noexcept
{
	MimmoFvMesh::swap(x);
	std::swap(m_type, x.m_type);
	std::swap(m_path, x.m_path);
	std::swap(m_overwrite, x.m_overwrite);
	std::swap(m_fieldname, x.m_fieldname);
	std::swap(m_field, x.m_field);
	std::swap(m_boundaryField, x.m_boundaryField);
};

/*! It builds the input/output ports of the object
 */
void
IOOFOAMScalarField::buildPorts(){
	IOOFOAMField::buildPorts();
	bool built = m_arePortsBuilt;
	// creating output ports
	built = (built && createPortOut<dmpvector1D, IOOFOAMScalarField>(this, &IOOFOAMScalarField::getField, M_SCALARFIELD));
	built = (built && createPortOut<dmpvector1D, IOOFOAMScalarField>(this, &IOOFOAMScalarField::getBoundaryField, M_SCALARFIELD2));
	m_arePortsBuilt = built;
};

/*!
 * It gets the internal scalar field.
 * \return internal scalar field.
 */
dmpvector1D
IOOFOAMScalarField::getField(){
	return m_field;
}

/*!
 * It gets the boundary scalar field.
 * \return boundary scalar field.
 */
dmpvector1D
IOOFOAMScalarField::getBoundaryField(){
	return m_boundaryField;
}

/*!
 * It reads the OpenFOAM field from input file and store in related variables m_field and m_boundaryfield.
 * \return false if errors occured during the reading.
 */
bool
IOOFOAMScalarField::read(){

	//read mesh from OpenFoam case directory (initialize)
	Foam::Time *foamRunTime = 0;
	Foam::fvMesh *foamMesh = 0;
	foamUtilsNative::initializeCase(m_path.c_str(), &foamRunTime, &foamMesh);


	if ( getGeometry() != nullptr){
		std::size_t size;
		std::vector<double> field;
		foamUtilsNative::readScalarField(m_path.c_str(), m_fieldname.c_str(), -1, size, field);
		m_field.clear();
		m_field.reserve(size);

		auto itfield = field.begin();
		for (bitpit::Cell & cell : getGeometry()->getCells()){
			m_field.insert(cell.getId(), *itfield);
			itfield++;
		}
		m_field.setGeometry(getGeometry());
		//Data on cells
		m_field.setDataLocation(2);

		//TODO ADD CELL TO POINT INTERPOLATION

	}

	//One field stored.
	if ( getBoundaryGeometry() != nullptr ){

		const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
		long startIndex;
		long endIndex;

		std::unordered_set<short> pids = getBoundaryGeometry()->getPIDTypeList();
		dmpvector1D boundaryFieldOnFace;
		for (short pid : pids){
			if (pid > 0){
				std::size_t size = 0;
				std::vector<double> field;
				foamUtilsNative::readScalarField(m_path.c_str(), m_fieldname.c_str(), pid-1, size, field);
				boundaryFieldOnFace.reserve(boundaryFieldOnFace.size() + size);
				if (size > 0){
					long iBoundary = pid-1;
					startIndex = foamBMesh[iBoundary].start();
					endIndex = startIndex + long(foamBMesh[iBoundary].size());
					long ind = startIndex;
					for (double val : field){
						boundaryFieldOnFace.insert(m_OFbitpitmapfaces[ind], val);
						ind++;
					}
				}
				else{
					for (bitpit::Cell cell : getBoundaryGeometry()->getCells()){
						if (cell.getPID() == pid)
							boundaryFieldOnFace.insert(cell.getId(), 0.);
					}
				}
			}
			else{
				for (bitpit::Cell cell : getBoundaryGeometry()->getCells()){
					if (cell.getPID() == pid){
						if (!boundaryFieldOnFace.exists(cell.getId()))
							boundaryFieldOnFace.insert(cell.getId(), 0.);
					}
				}
			}
		}
		boundaryFieldOnFace.setGeometry(getBoundaryGeometry());
		boundaryFieldOnFace.setDataLocation(1);
		interpolateFaceToPoint(boundaryFieldOnFace, m_boundaryField);
	}

	//TODO exception for null or empty geometries and return false for error during reading
	return true;

}

/*!
 * It writes the OpenFOAM field to an output file from internal structure m_field and m_boundaryfield
 * \return false if errors occured during the writing.
 */
bool
IOOFOAMScalarField::write(){
	if(!checkMeshCoherence()) return false;

	//do nothing for now

	return true;
}



/*
 * ========================================================================================================
 */

/*!
 * Default constructor of IOOFOAMVectorField.
 * \param[in] type int from enum IOOFMode, for class mode: READ, WRITE. Default is READ.
 */
IOOFOAMVectorField::IOOFOAMVectorField(int type):IOOFOAMField(){
	m_name           = "mimmo.IOOFOAMVectorField";
	IOOFOAMField::setDefaults();
}


/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAMVectorField::IOOFOAMVectorField(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.IOOFOAMVectorField";
	std::string fallback_name = "ClassNONE";
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);

	std::string fallback_type = "IOModeNONE";
	std::string input_type = rootXML.get("IOMode", fallback_type);
	input_type = bitpit::utils::string::trim(input_type);

	auto maybeIOOFMode = IOOFMode::_from_string_nothrow(input_type.c_str());

	if(input == "mimmo.IOOFOAMVectorField" && maybeIOOFMode){
		m_type = maybeIOOFMode->_to_integral();
		IOOFOAMField::setDefaults();
		IOOFOAMField::absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of IOOFOAMVectorField.
 */
IOOFOAMVectorField::~IOOFOAMVectorField(){};

/*!Copy constructor of IOOFOAMVectorField. Internal volume and boundary mesh are not copied.
 */
IOOFOAMVectorField::IOOFOAMVectorField(const IOOFOAMVectorField & other):IOOFOAMField(other){
	m_type = other.m_type;
	m_path = other.m_path;
	m_fieldname = other.m_fieldname;
	m_field = other.m_field;
	m_boundaryField = other.m_boundaryField;
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
};

/*!
 * Assignment operator. Internal volume and boundary mesh are not copied.
 */
IOOFOAMVectorField & IOOFOAMVectorField::operator=(IOOFOAMVectorField other){
	swap(other);
	return *this;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void IOOFOAMVectorField::swap(IOOFOAMVectorField & x) noexcept
{
	MimmoFvMesh::swap(x);
	std::swap(m_type, x.m_type);
	std::swap(m_path, x.m_path);
	std::swap(m_overwrite, x.m_overwrite);
	std::swap(m_fieldname, x.m_fieldname);
	std::swap(m_field, x.m_field);
	std::swap(m_boundaryField, x.m_boundaryField);
};

/*! It builds the input/output ports of the object
 */
void
IOOFOAMVectorField::buildPorts(){
	IOOFOAMField::buildPorts();
	bool built = m_arePortsBuilt;
	// creating output ports
	built = (built && createPortOut<dmpvecarr3E, IOOFOAMVectorField>(this, &IOOFOAMVectorField::getField, M_VECTORFIELD));
	built = (built && createPortOut<dmpvecarr3E, IOOFOAMVectorField>(this, &IOOFOAMVectorField::getBoundaryField, M_VECTORFIELD2));
	m_arePortsBuilt = built;
};

/*!
 * It gets the internal vector field.
 * \return internal vector field.
 */
dmpvecarr3E
IOOFOAMVectorField::getField(){
	return m_field;
}

/*!
 * It gets the boundary vector field.
 * \return boundary vector field.
 */
dmpvecarr3E
IOOFOAMVectorField::getBoundaryField(){
	return m_boundaryField;
}

/*!
 * It reads the OpenFOAM field from input file and store in related variables m_field and m_boundaryfield.
 * \return false if errors occured during the reading.
 */
bool
IOOFOAMVectorField::read(){

	Foam::Time *foamRunTime = 0;
	Foam::fvMesh *foamMesh = 0;

	//read mesh from OpenFoam case directory
	foamUtilsNative::initializeCase(m_path.c_str(), &foamRunTime, &foamMesh);

	if ( getGeometry() != nullptr){
		std::size_t size;
		std::vector<std::array<double,3>> field;
		foamUtilsNative::readVectorField(m_path.c_str(), m_fieldname.c_str(), -1, size, field);
		m_field.clear();
		m_field.reserve(size);

		auto itfield = field.begin();
		for (bitpit::Cell & cell : getGeometry()->getCells()){
			m_field.insert(cell.getId(), *itfield);
			itfield++;
		}
		m_field.setGeometry(getGeometry());
		m_field.setDataLocation(2);
	}

	//One field stored.
	if ( getBoundaryGeometry() != nullptr ){

		const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
		long startIndex;
		long endIndex;

		std::unordered_set<short> pids = getBoundaryGeometry()->getPIDTypeList();
		dmpvecarr3E boundaryFieldOnFace;
		for (short pid : pids){
			if (pid > 0){
				std::size_t size = 0;
				std::vector<std::array<double,3>> field;
				foamUtilsNative::readVectorField(m_path.c_str(), m_fieldname.c_str(), pid-1, size, field);
				boundaryFieldOnFace.reserve(boundaryFieldOnFace.size() + size);
				if (size > 0){
					long iBoundary = pid-1;
					startIndex = foamBMesh[iBoundary].start();
					endIndex = startIndex + long(foamBMesh[iBoundary].size());
					long ind = startIndex;
					for (std::array<double,3> val : field){
						boundaryFieldOnFace.insert(m_OFbitpitmapfaces[ind], val);
						ind++;
					}
				}
				else{
					for (bitpit::Cell cell : getBoundaryGeometry()->getCells()){
						if (cell.getPID() == pid)
							boundaryFieldOnFace.insert(cell.getId(), {{0.,0.,0.}});
					}
				}
			}
		}
		boundaryFieldOnFace.setGeometry(getBoundaryGeometry());
		boundaryFieldOnFace.setDataLocation(1);
		interpolateFaceToPoint(boundaryFieldOnFace, m_boundaryField);
	}

	//TODO exception for null or empty geometries and return false for error during reading
	return true;

}

/*!
 * It writes the OpenFOAM field to an output file from internal structure m_field and m_boundaryfield
 * \return false if errors occured during the writing.
 */
bool
IOOFOAMVectorField::write(){
	if(!checkMeshCoherence()) return false;

	//do nothing for now

	return true;
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




}
