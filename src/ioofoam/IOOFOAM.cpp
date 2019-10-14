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
 * Default constructor of IOOFOAM_Kernel.
 * \param[in] type int from enum IOOFMode. Default is READ.
 */
IOOFOAM_Kernel::IOOFOAM_Kernel(int type):MimmoFvMesh(){
	m_type = IOOFMode::READ;
    auto maybeIOMode = IOOFMode::_from_integral_nothrow(type);
    if(maybeIOMode){
		m_type = type;
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
 *\param[in] type int from enum IOOFMode. Default is WRITE.
 *
 */
IOOFOAM_Kernel::IOOFOAM_Kernel(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type): MimmoFvMesh(bulk,boundary){
	auto maybeIOMode = IOOFMode::_from_integral_nothrow(type);
    m_type = IOOFMode::WRITE;
	if(maybeIOMode){
		if(maybeIOMode->_to_integral() == IOOFMode::READ){
			throw std::runtime_error("Error in IOOFOAM constructor. Forcing READ mode in a constructor meant for writing modes only");
		}
		m_type = type;
	}
	setDefaults();
}

/*!Default destructor of IOOFOAM_Kernel.
 */
IOOFOAM_Kernel::~IOOFOAM_Kernel(){};


/*!
    Convenience copy constructor
*/
IOOFOAM_Kernel::IOOFOAM_Kernel(const IOOFOAM_Kernel & other): MimmoFvMesh(other){

    m_OFE_supp.clear();
    m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
    m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
    m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
    m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;

    m_type = other.m_type;
    m_name = other.m_name;
    m_path = other.m_path;
    m_fieldname = other.m_fieldname;
    m_OFbitpitmapfaces = other.m_OFbitpitmapfaces;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void IOOFOAM_Kernel::swap(IOOFOAM_Kernel & x) noexcept
{
	std::swap(m_type, x.m_type);
	std::swap(m_path, x.m_path);
	std::swap(m_fieldname, x.m_fieldname);
	std::swap(m_OFbitpitmapfaces, x.m_OFbitpitmapfaces);

	MimmoFvMesh::swap(x);
};

/*!
 * Default values for IOOFOAM_Kernel.
 */
void
IOOFOAM_Kernel::setDefaults(){

	m_path   = ".";
    m_fieldname = "";
	m_OFE_supp.clear();
	m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
	m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
	m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
	m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAM_Kernel::buildPorts(){
	bool built = true;
	// creating input ports
	built = (built && createPortIn<MimmoObject*, IOOFOAM_Kernel>(this, &IOOFOAM_Kernel::setGeometry, M_GEOMOFOAM));
	built = (built && createPortIn<MimmoObject*, IOOFOAM_Kernel>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOMOFOAM2));
    built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAM_Kernel>(this, &IOOFOAM_Kernel::setFacesMap, M_UMAPIDS));

	// creating output ports
	built = (built && createPortOut<MimmoObject*, IOOFOAM_Kernel>(this, &MimmoFvMesh::getGeometry, M_GEOMOFOAM));
	built = (built && createPortOut<MimmoObject*, IOOFOAM_Kernel>(this, &MimmoFvMesh::getBoundaryGeometry, M_GEOMOFOAM2));
	built = (built && createPortOut<std::unordered_map<long,long>, IOOFOAM_Kernel>(this, &IOOFOAM_Kernel::getFacesMap, M_UMAPIDS));

	m_arePortsBuilt = built;
};

/*!
 * Get current map between OFOAM faces -> bitpit interfaces.
 * \return reference to interfaces map.
 */
std::unordered_map<long,long>
IOOFOAM_Kernel::getFacesMap(){
	return m_OFbitpitmapfaces;
}

/*!
 * Get type parameter - int casting of IOOFMode enum.
 * \return type parameter content
 */
int
IOOFOAM_Kernel::getType(){
	return m_type;
}

/*!
 * It sets the name of directory to read/write the OpenFOAM mesh.
 * \param[in] dir mesh input directory.
 */
void
IOOFOAM_Kernel::setDir(const std::string &dir){
	m_path = dir;
}

/*!
 * Set current bulk geometry MimmoObject mesh.
 * The mesh must be of volume type. Internal check to know if interfaces are built is done,
 * since without them can be not coherent link with OpenFoam mesh from file (setDir method).
 * Any previous mesh internally allocated will be destroyed. See MimmoFvMesh::setGeometry.
 * \param[in] bulk pointer to external bulk mesh MimmoObject.
 */
void
IOOFOAM_Kernel::setGeometry(MimmoObject * bulk){
    if(!bulk) return;
    if(!bulk->areInterfacesBuilt()) {
        *(m_log)<<"Warning IOOFOAM_Kernel:: linked MimmoObject bulk mesh cannot be coherent with an OpenFoam mesh"<<std::endl;
    }
	MimmoFvMesh::setGeometry(bulk);
}

/*!
 * Set current boundary geometry -> mesh coherent with bulk set with setGeometry();
 * Any previous mesh internally allocated will be destroyed. See MimmoFvMesh::setBoundaryGeometry
 * \param[in] boundary pointer to external boundary mesh MimmoObject.
 */
void
IOOFOAM_Kernel::setBoundaryGeometry(MimmoObject * boundary){
    if(!boundary) return;
	MimmoFvMesh::setBoundaryGeometry(boundary);
}

/*!
 * Set current map between OFOAM faces -> bitpit interfaces.
 * \param[in] mapFaces reference to interfaces map.
 */
void
IOOFOAM_Kernel::setFacesMap(std::unordered_map<long,long> mapFaces){
    m_OFbitpitmapfaces = mapFaces;
}

/*!
 * It sets the name of the field to read/write.
 * \param[in] fieldname name of input/output field.
 */
void
IOOFOAM_Kernel::setFieldName(const std::string & fieldname){
	m_fieldname = fieldname;
}

/*!
 * Set mode of the class, IOOFMode enum.
 * \param[in] type IOOFMode enum
 */
void
IOOFOAM_Kernel::setType(IOOFMode type){
	m_type = type;
}

/*!
 * Set mode of the class, IOOFMode enum.
 * \param[in] type int casting to IOOFMode enum
 */
void
IOOFOAM_Kernel::setType(int type){
    auto maybe = IOOFMode::_from_integral_nothrow(type);
    if(maybe)   setType(*maybe);
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM_Kernel::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	std::string input;

	BaseManipulation::absorbSectionXML(slotXML, name);

	if(slotXML.hasOption("Dir")){
		input = slotXML.get("Dir");
		input = bitpit::utils::string::trim(input);
		if(input.empty())   input = ".";
		setDir(input);
	};
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM_Kernel::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	BaseManipulation::flushSectionXML(slotXML, name);

	slotXML.set("IOMode", IOOFMode::_from_integral(m_type)._to_string());
	slotXML.set("Dir", m_path);
};

/*
 * ========================================================================================================
 */

/*!
    Default constructor
*/
IOOFOAM::IOOFOAM(int type):IOOFOAM_Kernel(type){
    m_name = "mimmo.IOOFOAM";
    m_overwrite=false;
};

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
 *\param[in] type int from enum IOOFMode. Default is WRITE.
 *
 */

IOOFOAM::IOOFOAM(std::unique_ptr<MimmoObject> & bulk, std::unique_ptr<MimmoObject> &boundary, int type):
         IOOFOAM_Kernel(bulk,boundary, type)
{
    m_name = "mimmo.IOOFOAM";
    m_overwrite= false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAM::IOOFOAM(const bitpit::Config::Section & rootXML){

	m_name = "mimmo.IOOFOAM";
    m_overwrite= false;
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

/*!
   Destructor of IOOFOAM
 */
IOOFOAM::~IOOFOAM(){};

/*!
Copy constructor of IOOFOAM. Internal volume and boundary mesh are not copied.
 */
IOOFOAM::IOOFOAM(const IOOFOAM & other):IOOFOAM_Kernel(other){
    setDefaults();
	m_type = other.m_type;
	m_path = other.m_path;
	m_overwrite = other.m_overwrite;
	m_OFbitpitmapfaces = other.m_OFbitpitmapfaces;
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
	std::swap(m_overwrite, x.m_overwrite);
	IOOFOAM_Kernel::swap(x);
};

/*!
 * Default values for IOOFOAM_Kernel.
 */
void
IOOFOAM::setDefaults(){
	m_overwrite = false;
    IOOFOAM_Kernel::setDefaults();
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
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAM::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);

	std::string input;

	IOOFOAM_Kernel::absorbSectionXML(slotXML, name);
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

	IOOFOAM_Kernel::flushSectionXML(slotXML, name);
	slotXML.set("Overwrite", std::to_string(m_overwrite));
};

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
	long PID = 0;
	livector1D conn, temp;

	forAll(cells, iC){

		iDC = long(iC);
		eleshape = std::string(cellShapes[iC].model().name());
		//first step verify the model in twin cellShapes list.
		conn.clear();

		if(m_OFE_supp.count(eleshape) > 0){

			eltype = m_OFE_supp[eleshape];
			temp.resize(cellShapes[iC].size());
			forAll(cellShapes[iC], loc){
				temp[loc] = (long)cellShapes[iC][loc];
			}
			conn = foamUtilsNative::mapEleVConnectivity(temp, eltype);

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
					normalIsOut = faceNeighbour[iFace] > iC;
				}else{
					normalIsOut = faceOwner[iFace] > iC;
				}

				if(normalIsOut){
					conn.insert(conn.end(), temp.begin(), temp.end());
				}else{
					conn.insert(conn.end(), temp.rbegin(), temp.rend());
				}
			}
		}
		bitpit::PatchKernel::CellIterator it = mesh->addCell(eltype, conn, iDC);
		it->setPID(int(PID));
	}

    mesh->buildAdjacencies();
    mesh->buildInterfaces();

    bitpit::PiercedVector<bitpit::Cell> & bitCells = mesh->getCells();
    bitpit::PiercedVector<bitpit::Interface> & bitInterfaces = mesh->getInterfaces();

	forAll(faces, iOF){

        std::vector<long> vListOF(faces[iOF].size());
        forAll(faces[iOF], index){
            vListOF[index] = long(faces[iOF][index]);
        }
        std::sort(vListOF.begin(), vListOF.end());

        long iDC = long(faceOwner[iOF]);
		long * bitFaceList = bitCells[iDC].getInterfaces();
		std::size_t sizeFList = bitCells[iDC].getInterfaceCount();

        long iBIT = bitpit::Interface::NULL_ID;
        std::size_t j(0);
        while(iBIT < 0 && j<sizeFList){
            long * vconn = bitInterfaces[bitFaceList[j]].getConnect();
            std::size_t vconnsize = bitInterfaces[bitFaceList[j]].getConnectSize();
            std::vector<long> vListBIT(vconn, vconn+vconnsize);
            std::sort(vListBIT.begin(), vListBIT.end());

            if(std::equal(vListBIT.begin(), vListBIT.end(), vListOF.begin()) ){
                iBIT = bitFaceList[j];
            }else{
                ++j;
            }
        }
		m_OFbitpitmapfaces.insert(std::make_pair(long(iOF),iBIT) );
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
		PID = long(iBoundary+1);
		startIndex = foamBMesh[iBoundary].patch().start();
		endIndex = startIndex + long(foamBMesh[iBoundary].patch().size());
        for(long ind=startIndex; ind<endIndex; ++ind){
			m_boundary->setPIDCell(m_OFbitpitmapfaces[ind], PID);
		}
        m_boundary->setPIDName(PID,foamBMesh[iBoundary].name().c_str());
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

	dvecarr3E points = getGeometry()->getVerticesCoords();

	return foamUtilsNative::writePointsOnCase(m_path.c_str(), points, m_overwrite);
}

/*
 * ========================================================================================================
 */

/*!
 * Default constructor of IOOFOAMScalarField.
 * \param[in] type int from enum IOOFMode, for class mode: READ, WRITE. Default is READ.
 */
IOOFOAMScalarField::IOOFOAMScalarField(int type):IOOFOAM_Kernel(type){
	m_name           = "mimmo.IOOFOAMScalarField";
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
		setDefaults();
		absorbSectionXML(rootXML);
	}else{
		warningXML(m_log, m_name);
	};
}

/*!Default destructor of IOOFOAMScalarField.
 */
IOOFOAMScalarField::~IOOFOAMScalarField(){};

/*!Copy constructor of IOOFOAMScalarField. Internal volume and boundary mesh are not copied.
 */
IOOFOAMScalarField::IOOFOAMScalarField(const IOOFOAMScalarField & other):IOOFOAM_Kernel(other){
    m_field = other.m_field;
    m_boundaryField = other.m_boundaryField;
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
	m_field.swap(x.m_field);
	m_boundaryField.swap(x.m_boundaryField);
    IOOFOAM_Kernel::swap(x);
};

/*! It builds the input/output ports of the object
 */
void
IOOFOAMScalarField::buildPorts(){
	IOOFOAM_Kernel::buildPorts();
	bool built = m_arePortsBuilt;
	// creating output ports
	built = (built && createPortOut<dmpvector1D*, IOOFOAMScalarField>(this, &IOOFOAMScalarField::getField, M_SCALARFIELD));
	built = (built && createPortOut<dmpvector1D*, IOOFOAMScalarField>(this, &IOOFOAMScalarField::getBoundaryField, M_SCALARFIELD2));
	m_arePortsBuilt = built;
};

/*!
 * It gets the internal scalar field.
 * \return internal scalar field.
 */
dmpvector1D *
IOOFOAMScalarField::getField(){
	return &m_field;
}

/*!
 * It gets the boundary scalar field.
 * \return boundary scalar field.
 */
dmpvector1D *
IOOFOAMScalarField::getBoundaryField(){
	return &m_boundaryField;
}

/*!Execution command.
 * It reads the geometry in reading mode, otherwise it writes, according to the write mode.
 */
void
IOOFOAMScalarField::execute(){
	bool check = true;
	switch (m_type){
	case IOOFMode::READ :
		check = read();
		break;
	case IOOFMode::WRITE :
		//             check = write(); //TODO coding complete mesh writing from scratch.
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
IOOFOAMScalarField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	std::string input;
	IOOFOAM_Kernel::absorbSectionXML(slotXML, name);

	if(slotXML.hasOption("FieldName")){
		input = slotXML.get("FieldName");
		input = bitpit::utils::string::trim(input);
		if(input.empty())   input = ".";
		setFieldName(input);
	};
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOOFOAMScalarField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	IOOFOAM_Kernel::flushSectionXML(slotXML, name);

	slotXML.set("FieldName", m_fieldname);
};

/*!
 * It reads the OpenFOAM field from input file and store in related variables m_field and m_boundaryField.
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

		std::unordered_set<long> pids = getBoundaryGeometry()->getPIDTypeList();
		dmpvector1D boundaryFieldOnFace;
		for (long pid : pids){
			if (pid > 0){
				std::size_t size = 0;
				std::vector<double> field;
				foamUtilsNative::readScalarField(m_path.c_str(), m_fieldname.c_str(), pid-1, size, field);
				boundaryFieldOnFace.reserve(boundaryFieldOnFace.size() + size);
				if (size > 0){
					long iBoundary = pid-1;
					startIndex = foamBMesh[iBoundary].start();
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
		foamUtilsNative::interpolateFaceToPoint(boundaryFieldOnFace, m_boundaryField);
	}

	//TODO exception for null or empty geometries and return false for error during reading
	return true;

}

/*!
 * It writes the OpenFOAM field to an output file from internal structure m_field and m_boundaryField
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
 IOOFOAMVectorField::IOOFOAMVectorField(int type):IOOFOAM_Kernel(type){
 	m_name           = "mimmo.IOOFOAMVectorField";
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
 		setDefaults();
 		absorbSectionXML(rootXML);
 	}else{
 		warningXML(m_log, m_name);
 	};
 }

 /*!Default destructor of IOOFOAMVectorField.
  */
 IOOFOAMVectorField::~IOOFOAMVectorField(){};

 /*!Copy constructor of IOOFOAMVectorField. Internal volume and boundary mesh are not copied.
  */
 IOOFOAMVectorField::IOOFOAMVectorField(const IOOFOAMVectorField & other):IOOFOAM_Kernel(other){
     m_field = other.m_field;
     m_boundaryField = other.m_boundaryField;
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
 	m_field.swap(x.m_field);
 	m_boundaryField.swap(x.m_boundaryField);
     IOOFOAM_Kernel::swap(x);
 };

 /*! It builds the input/output ports of the object
  */
 void
 IOOFOAMVectorField::buildPorts(){
 	IOOFOAM_Kernel::buildPorts();
 	bool built = m_arePortsBuilt;
 	// creating output ports
 	built = (built && createPortOut<dmpvecarr3E*, IOOFOAMVectorField>(this, &IOOFOAMVectorField::getField, M_VECTORFIELD));
 	built = (built && createPortOut<dmpvecarr3E*, IOOFOAMVectorField>(this, &IOOFOAMVectorField::getBoundaryField, M_VECTORFIELD2));
 	m_arePortsBuilt = built;
 };

 /*!
  * It gets the internal vector field.
  * \return internal vector field.
  */
 dmpvecarr3E *
 IOOFOAMVectorField::getField(){
 	return &m_field;
 }

 /*!
  * It gets the boundary scalar field.
  * \return boundary scalar field.
  */
 dmpvecarr3E *
 IOOFOAMVectorField::getBoundaryField(){
 	return &m_boundaryField;
 }

 /*!Execution command.
  * It reads the geometry in reading mode, otherwise it writes, according to the write mode.
  */
 void
 IOOFOAMVectorField::execute(){
 	bool check = true;
 	switch (m_type){
 	case IOOFMode::READ :
 		check = read();
 		break;
 	case IOOFMode::WRITE :
 		//             check = write(); //TODO coding complete mesh writing from scratch.
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
 IOOFOAMVectorField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

 	BITPIT_UNUSED(name);
 	std::string input;
 	IOOFOAM_Kernel::absorbSectionXML(slotXML, name);

 	if(slotXML.hasOption("FieldName")){
 		input = slotXML.get("FieldName");
 		input = bitpit::utils::string::trim(input);
 		if(input.empty())   input = ".";
 		setFieldName(input);
 	};
 };

 /*!
  * It sets infos from class members in a XML bitpit::Config::section.
  * \param[in] slotXML bitpit::Config::Section of XML file
  * \param[in] name   name associated to the slot
  */
 void
 IOOFOAMVectorField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

 	BITPIT_UNUSED(name);
 	IOOFOAM_Kernel::flushSectionXML(slotXML, name);

 	slotXML.set("FieldName", m_fieldname);
 };

/*!
 * It reads the OpenFOAM field from input file and store in related variables m_field and m_boundaryField.
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

		std::unordered_set<long> pids = getBoundaryGeometry()->getPIDTypeList();
		dmpvecarr3E boundaryFieldOnFace;
		for (long pid : pids){
			if (pid > 0){
				std::size_t size = 0;
				std::vector<std::array<double,3>> field;
				foamUtilsNative::readVectorField(m_path.c_str(), m_fieldname.c_str(), pid-1, size, field);
				boundaryFieldOnFace.reserve(boundaryFieldOnFace.size() + size);
				if (size > 0){
					long iBoundary = pid-1;
					startIndex = foamBMesh[iBoundary].start();
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
		foamUtilsNative::interpolateFaceToPoint(boundaryFieldOnFace, m_boundaryField);
	}

	//TODO exception for null or empty geometries and return false for error during reading
	return true;

}

/*!
 * It writes the OpenFOAM field to an output file from internal structure m_field and m_boundaryField
 * \return false if errors occured during the writing.
 */
bool
IOOFOAMVectorField::write(){
	if(!checkMeshCoherence()) return false;

	//do nothing for now

	return true;
}








}
