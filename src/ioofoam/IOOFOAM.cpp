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
#include <processorCyclicFvPatch.H>

namespace mimmo{

/*!
 * Default constructor of IOOFOAM_Kernel.
 * \param[in] writemode if true, put the class in write mode, otherwise in read mode.
 */
IOOFOAM_Kernel::IOOFOAM_Kernel(bool writemode):BaseManipulation(){
	m_write = writemode;
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
 *
 *\param[in] bulk shared pointer to the bulk MimmoObject
 *\param[in] boundary shared pointer to the boundary MimmoObject
 *\param[in] writemode if true, put the class in write mode, otherwise in read mode.
 *
 */
IOOFOAM_Kernel::IOOFOAM_Kernel(MimmoSharedPointer<MimmoObject> bulk, MimmoSharedPointer<MimmoObject> boundary, bool writemode){
	m_write = writemode;
	setDefaults();

	// Check input geometries
	if(!bulk && !boundary){
	    throw std::runtime_error("Error in MimmoFvMesh Constructor: nullptr linked as bulk mesh or boundary mesh");
	}
	std::string key = std::to_string(bulk->getType()) + std::to_string(boundary->getType());

	if(key !="21" && key !="14"){
	    throw std::runtime_error("Error in MimmoFvMesh Constructor: not valid types in pair bulk-boundary meshes");
	}

    setGeometry(bulk);
    setBoundaryGeometry(boundary);


}

/*!Default destructor of IOOFOAM_Kernel.
 */
IOOFOAM_Kernel::~IOOFOAM_Kernel(){};


/*!
    Convenience copy constructor
*/
IOOFOAM_Kernel::IOOFOAM_Kernel(const IOOFOAM_Kernel & other): BaseManipulation(other){

    m_OFE_supp.clear();
    m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
    m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
    m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
    m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;

    m_write = other.m_write;
    m_name = other.m_name;
    m_path = other.m_path;
    m_fieldname = other.m_fieldname;
    m_OFbitpitmapfaces = other.m_OFbitpitmapfaces;
    m_boundary = other.m_boundary;
}

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void IOOFOAM_Kernel::swap(IOOFOAM_Kernel & x) noexcept
{
	std::swap(m_write, x.m_write);
	std::swap(m_path, x.m_path);
	std::swap(m_fieldname, x.m_fieldname);
    std::swap(m_OFbitpitmapfaces, x.m_OFbitpitmapfaces);
    std::swap(m_boundary, x.m_boundary);

	BaseManipulation::swap(x);
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
	m_boundary.reset();
}

/*!
 * Get current map between OFOAM faces -> bitpit interfaces.
 * \return reference to interfaces map.
 */
std::unordered_map<long,long>
IOOFOAM_Kernel::getFacesMap(){
	return m_OFbitpitmapfaces;
}

/*!
 * \return if the class is in writing or reading mode.
 */
bool
IOOFOAM_Kernel::isWriteMode(){
	return m_write;
}

/*!
 * Get current boundary geometry coherent with bulk.
 * \return boundary pointer to internal boundary mesh MimmoObject.
 */
MimmoSharedPointer<MimmoObject>
IOOFOAM_Kernel::getBoundaryGeometry(){
    return m_boundary;
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
 * Any previous mesh internally allocated will be destroyed. See BaseManipulation::setGeometry.
 * \param[in] bulk pointer to external bulk mesh MimmoObject.
 */
void
IOOFOAM_Kernel::setGeometry(MimmoSharedPointer<MimmoObject> bulk){
    if(!bulk) return;
    if(bulk->getInterfacesSyncStatus() != SyncStatus::SYNC) {
        *(m_log)<<"Warning IOOFOAM_Kernel:: linked MimmoObject bulk mesh can be not coherent with an OpenFoam mesh"<<std::endl;
    }
    if (bulk->getType() != 2 && bulk->getType() != 1){
        *(m_log)<<"Warning IOOFOAM_Kernel:: linked MimmoObject bulk mesh can be only volume or surface mesh. Continue without bulk geometry."<<std::endl;
        return;
    }
	BaseManipulation::setGeometry(bulk);
}

/*!
 * Set current boundary geometry -> mesh coherent with bulk set with setGeometry();
 * Any previous mesh internally allocated will be destroyed.
 * \param[in] boundary pointer to external boundary mesh MimmoObject.
 */
void
IOOFOAM_Kernel::setBoundaryGeometry(MimmoSharedPointer<MimmoObject> boundary){
    if(!boundary) return;
    if (boundary->getType() != 1 && boundary->getType() != 4){
        *(m_log)<<"Warning IOOFOAM_Kernel:: linked MimmoObject boundary mesh can be only surface or 3D-curve mesh. Continue without boundary geometry."<<std::endl;
        return;
    }
    m_boundary = boundary;
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

	slotXML.set("IOMode", std::to_string(m_write));
	slotXML.set("Dir", m_path);
};


/*!
 * Check coeherence between bulk and boundary mesh, i.e.:
 *
 * - bulk non-null and non-empty
 * - if boundary non null and non empty, it checks out:
 *   - bulk/boundary pair must be of MimmoObject types 2-Volume bulk and 1-Surface boundary or 1-Surface bulk and 4-3Dcurve boundary
 *   - ids tight matching between border interfaces on the bulk and cells of the boundary
 * - otherwise:
 *   - create internally a boundary patch, without pid segmentation informations, and return true.
 *
 * \return false if one of the check failed.
 */
bool
IOOFOAM_Kernel::checkMeshCoherence(){

    MimmoSharedPointer<MimmoObject> bulk = getGeometry();
    if(bulk == nullptr)  return false;

    MimmoSharedPointer<MimmoObject> boundary = getBoundaryGeometry();

    bool checkBoundaries = true;
    if(boundary == nullptr)  checkBoundaries = false;

    if(checkBoundaries){

        if(bulk->getInterfacesSyncStatus() != SyncStatus::SYNC)  return false;
        std::string key = std::to_string(bulk->getType()) + std::to_string(boundary->getType());
        if(key !="21" && key !="14")    return false;

        std::set<long> idBorderInterf;
        for(const auto & interf : bulk->getInterfaces()){
            if(interf.isBorder())   idBorderInterf.insert(interf.getId());
        };
        std::vector<long> idInterf;
        idInterf.insert(idInterf.end(), idBorderInterf.begin(), idBorderInterf.end());
        std::vector<long> idCell   = boundary->getCells().getIds(true);

        return idInterf == idCell;
    }else{
        m_geometry->extractBoundaryMesh();
        return true;
    }
}

/*
 * ========================================================================================================
 */



/*!
    Default constructor
    \param[in] writemode if true, put the class in write mode, otherwise in read mode.
*/
IOOFOAM::IOOFOAM(bool writemode):IOOFOAM_Kernel(writemode){
    m_name = "mimmo.IOOFOAM";
    m_overwrite = false;
    m_writepointsonly = true;
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
 * The class is meant for writing mode ONLY when built with this constructor.
 *
 *\param[in] bulk shared pointer to the bulk MimmoObject
 *\param[in] boundary shared pointer to the boundary MimmoObject
 *
 */

IOOFOAM::IOOFOAM(MimmoSharedPointer<MimmoObject> bulk, MimmoSharedPointer<MimmoObject> boundary):
         IOOFOAM_Kernel(bulk,boundary, true)
{
    m_name = "mimmo.IOOFOAM";
    m_overwrite= false;
    m_writepointsonly = true;
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

    std::string fallback_name2 = "0";
	std::string input2 = rootXML.get("IOMode", fallback_name2);
	input2 = bitpit::utils::string::trim(input2);
    {
        std::stringstream ss(input2);
        ss>>m_write;
    }

	if(input == "mimmo.IOOFOAM"){
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
	m_write = other.m_write;
    m_writepointsonly = other.m_writepointsonly;
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
    std::swap(m_writepointsonly, x.m_writepointsonly);
	IOOFOAM_Kernel::swap(x);
};

/*!
 * Default values for IOOFOAM_Kernel.
 */
void
IOOFOAM::setDefaults(){
	m_overwrite = false;
    m_writepointsonly = true;
    IOOFOAM_Kernel::setDefaults();
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAM::buildPorts(){
   bool built = true;
    // creating input ports
    // check if geometries are really mandatory (in case custom constructor with bulk and boundary is used.)
    bool m_mandbulk = (m_write && getGeometry() == nullptr);
    bool m_mandbndy = (m_write && getBoundaryGeometry() == nullptr);
    m_mandbndy = m_mandbndy && !m_writepointsonly; //boundary are not needed really if you write only points.

    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::setGeometry, M_GEOMOFOAM, m_mandbulk, 1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::setGeometry, M_GEOM, m_mandbulk, 1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOMOFOAM2, m_mandbndy, 2));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOM2, m_mandbndy, 2));
    built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAM>(this, &IOOFOAM_Kernel::setFacesMap, M_UMAPIDS));

    // creating output ports
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &BaseManipulation::getGeometry, M_GEOMOFOAM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &BaseManipulation::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOMOFOAM2));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAM>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOM2));
    built = (built && createPortOut<std::unordered_map<long,long>, IOOFOAM>(this, &IOOFOAM_Kernel::getFacesMap, M_UMAPIDS));

   m_arePortsBuilt = built;
};



/*!
 * Set writePointsOnly parameter. This option is valid only in write mode.
 * If true write mesh points only in a preexistent OpenFoam mesh case, if false
   write mesh from scratch.
 * \param[in] flag activation flag.
 */
void
IOOFOAM::setWritePointsOnly(bool flag){
	//m_writepointsonly = flag;
    m_writepointsonly = true; //TODO set to flag once write() is coded
}

/*!
 * Get WritePointsOnly parameter.  See setWritePointsOnly method.
 * \return writePointsOnly parameter content
 */
bool
IOOFOAM::getWritePointsOnly(){
	return m_writepointsonly;
}

/*!
 * Set overwrite parameter. This option is valid only in write mode and WritePointsOnly parameter enabled.
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
    if(slotXML.hasOption("WritePointsOnly")){
		input = slotXML.get("WritePointsOnly");
		bool value = true;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::string::trim(input));
			ss >> value;
		}
		setOverwrite(value);
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

	IOOFOAM_Kernel::flushSectionXML(slotXML, name);
    slotXML.set("WritePointsOnly", std::to_string(m_overwrite));
    slotXML.set("Overwrite", std::to_string(m_overwrite));
};

/*!Execution command.
 * It reads the geometry in reading mode, otherwise it writes, according to the write mode.
 */
void
IOOFOAM::execute(){
	bool check = true;
	if(!m_write){
        check = read();
    }else{
        if (!m_writepointsonly){
            check = write();
        }else{
            check = writePointsOnly();
        }
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

    //instantiate my bulk geometry container
    m_geometry = MimmoSharedPointer<MimmoObject>(new MimmoObject(2));
    MimmoSharedPointer<MimmoObject> mesh = m_geometry;

	Foam::Time *foamRunTime = 0;
	Foam::fvMesh *foamMesh = 0;

	//read mesh from OpenFoam case directory
	foamUtilsNative::initializeCase(m_path.c_str(), &foamRunTime, &foamMesh, m_geometry->getProcessorCount());

	//prepare my bulk geometry container
    mesh->getPatch()->reserveVertices(std::size_t(foamMesh->nPoints()));
	mesh->getPatch()->reserveCells(std::size_t(foamMesh->nCells()));

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
	    mesh->addConnectedCell(conn, eltype, PID, iDC);
	}


    std::vector<long> vectorVertexIds = mesh->getVertexFromCellList(mesh->getCellsIds(true));

#if MIMMO_ENABLE_MPI

	// Build ghost cells only in case of multi-processors run
	if (getProcessorCount() > 1){

	    // Build Ghost Cells from boundary cells of patches of type processor boundary

	    // Recovering boundary nodes for each partition
	    std::vector<long> boundary_vertices = mesh->extractBoundaryVertexID();

	    // Communicate boundary vertices to all other partitions
	    std::unordered_map<int, std::vector<long>*> boundaryVertexExchangeSources;
	    std::unordered_map<int, bitpit::PiercedVector<bitpit::Vertex> > receivedBoundaryVertices;

	    // Fill map to communicate
	    for (int irank = 0; irank < getProcessorCount(); irank++){
	        if (irank != getRank()){
	            boundaryVertexExchangeSources[irank] = &boundary_vertices;
	        }
	    }

	    {
	        // Send boundary vertices
	        std::unique_ptr<bitpit::DataCommunicator> dataCommunicator;
	        dataCommunicator = std::unique_ptr<bitpit::DataCommunicator>(new bitpit::DataCommunicator(getCommunicator()));

	        for (auto & entry : boundaryVertexExchangeSources){

	            int target_rank = entry.first;
	            std::vector<long> & vertexList = *entry.second;

	            // Prepare buffer for source cells and vertices of the source cells
	            size_t vertexBufferSize = 0;
	            for (long vertexId : vertexList){
	                vertexBufferSize += mesh->getPatch()->getVertex(vertexId).getBinarySize();
	            }
	            dataCommunicator->setSend(target_rank, sizeof(size_t) + vertexBufferSize);
	            bitpit::SendBuffer &buffer = dataCommunicator->getSendBuffer(target_rank);

	            // Fill the buffer with the source cells and vertices
	            buffer << vertexList.size();
	            for (long vertexId : vertexList){
	                buffer << mesh->getPatch()->getVertex(vertexId);
	            }
	            dataCommunicator->startSend(target_rank);

	        }

	        // Discover & start all the receives
	        dataCommunicator->discoverRecvs();
	        dataCommunicator->startAllRecvs();

	        // Receive the target boundary vertices
	        int nCompletedRecvs = 0;
	        while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
	            int source_rank = dataCommunicator->waitAnyRecv();
	            bitpit::RecvBuffer &buffer = dataCommunicator->getRecvBuffer(source_rank);

	            // Receive boundary vertices from mesh partitions
	            std::size_t vertexCount;
	            buffer >> vertexCount;
	            receivedBoundaryVertices[source_rank].reserve(vertexCount);
	            for (std::size_t ivertex = 0; ivertex < vertexCount; ivertex++){
	                bitpit::Vertex boundaryVertex;
	                buffer >> boundaryVertex;
	                receivedBoundaryVertices[source_rank].insert(boundaryVertex.getId(), boundaryVertex);
	            }

	            ++nCompletedRecvs;
	        }
	        // Wait for the sends to finish
	        dataCommunicator->waitAllSends();

	    }

	    // Update kdtree
	    if (mesh->getKdTreeSyncStatus() != SyncStatus::SYNC){
	        mesh->buildKdTree();
	    }

	    // Found list of vertices coincident
	    // with received boundary vertices for each partition

	    // Initialize ghost cells to communicate to each processor
	    std::unordered_map<int, long> ghostVertexExchangeSourcesCount;
	    std::unordered_map<int, std::unordered_set<long>> ghostCellExchangeSources;
	    std::unordered_map<int, std::unordered_set<long>> vertexCellExchangeSources;

	    // Loop over received bounday vertices for rank
	    for (int irank = 0; irank < getProcessorCount(); irank ++){

	        if (irank == getRank()) continue;

	        long received_vertices_count = receivedBoundaryVertices[irank].size();
	        std::vector<long> coincident_vertices;
	        coincident_vertices.reserve(received_vertices_count);

	        // Check received vertices with local kdtree
	        for (bitpit::Vertex vertex : receivedBoundaryVertices[irank]){

	            long coincident_vertex_id = bitpit::Vertex::NULL_ID;
	            mesh->getKdTree()->exist(&vertex, coincident_vertex_id);
	            if (coincident_vertex_id != bitpit::Vertex::NULL_ID){
	                coincident_vertices.push_back(coincident_vertex_id);
	            }
	            coincident_vertices.shrink_to_fit();

	        }

	        // Extract local source cells those are ghost cells for i-th rank
	        std::vector<long> ghosts = mesh->getCellFromVertexList(coincident_vertices, false);

	        // Fill communication structures
	        ghostVertexExchangeSourcesCount[irank] = 0;
	        for (long cellId : ghosts){
	            ghostCellExchangeSources[irank].insert(cellId);
	            for (long vertexId : mesh->getPatch()->getCell(cellId).getVertexIds()){
	                vertexCellExchangeSources[irank].insert(vertexId);
	                ghostVertexExchangeSourcesCount[irank]++;
	            }
	        }

	    }


	    {
	        // Send number of vertices
	        std::unique_ptr<bitpit::DataCommunicator> dataCommunicator;
	        dataCommunicator = std::unique_ptr<bitpit::DataCommunicator>(new bitpit::DataCommunicator(getCommunicator()));

	        for (auto & entry : ghostCellExchangeSources){

	            int target_rank = entry.first;
	            std::unordered_set<long> & cellList = entry.second;
	            std::unordered_set<long> & vertexList = vertexCellExchangeSources[target_rank];

	            // Prepare buffer for source cells and vertices of the source cells
	            long vertexCount = 0;
	            for (long vertexId : vertexList){
	                vertexCount++;
	            }
	            dataCommunicator->setSend(target_rank, vertexCount*sizeof(long));
	            bitpit::SendBuffer &buffer = dataCommunicator->getSendBuffer(target_rank);

	            // Fill the buffer with the source cells and vertices
	            buffer << vertexCount;
	            dataCommunicator->startSend(target_rank);
	        }

	        // Discover & start all the receives
	        dataCommunicator->discoverRecvs();
	        dataCommunicator->startAllRecvs();

	        // Receive the target count vertices

	        int nCompletedRecvs = 0;
	        while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
	            int source_rank = dataCommunicator->waitAnyRecv();
	            bitpit::RecvBuffer &buffer = dataCommunicator->getRecvBuffer(source_rank);

	            long vertexCount;
	            buffer >> vertexCount;

	            ghostVertexExchangeSourcesCount[source_rank] = vertexCount;

	            ++nCompletedRecvs;
	        }
	        // Wait for the sends to finish
	        dataCommunicator->waitAllSends();

	    } // end communication of receive vertex count


	    // Force re-build kdTree of vertices with new over-sized reserved vertex container
	    mesh->buildKdTree();

	    // Instantiate a to_insert structure for unique received vertices for each sending rank
	    // Initialize a kdTree tree for the received ghost vertices
	    std::unordered_map< long , std::vector<bitpit::Vertex> > unique_ghosts;
	    for (auto entry : ghostVertexExchangeSourcesCount){
	        unique_ghosts[entry.first].reserve(entry.second);
	    }
	    bitpit::KdTree<3, bitpit::Vertex, std::pair<int, long> > ghost_tree;

	    // Renumber connectivity of cells with new vertex ids for each sending rnak
	    std::unordered_map< long, std::unordered_map<long, long> > oldToNewVertexId;

	    // Some vertex are coincident with other ghost vertices already to be insert
	    // Link these vertices to renumber after the insertion of the linked ones of the related source rank
	    std::unordered_map< long, std::unordered_map<long, std::pair<int, long> > > oldToRenumberedGhostVertexId;

	    {
	        // Send vertices
	        std::unique_ptr<bitpit::DataCommunicator> dataCommunicator;
	        dataCommunicator = std::unique_ptr<bitpit::DataCommunicator>(new bitpit::DataCommunicator(getCommunicator()));

	        for (auto & entry : ghostCellExchangeSources){

	            int target_rank = entry.first;
	            std::unordered_set<long> & vertexList = vertexCellExchangeSources[target_rank];

	            // Prepare buffer for vertices of the source cells
	            size_t vertexBufferSize = 0;
	            for (long vertexId : vertexList){
	                vertexBufferSize += mesh->getPatch()->getVertex(vertexId).getBinarySize();
	            }
	            dataCommunicator->setSend(target_rank, sizeof(size_t) + vertexBufferSize);
	            bitpit::SendBuffer &buffer = dataCommunicator->getSendBuffer(target_rank);

	            // Fill the buffer with the source vertices
	            buffer << vertexList.size();
	            for (long vertexId : vertexList){
	                buffer << mesh->getPatch()->getVertex(vertexId);
	            }
	            dataCommunicator->startSend(target_rank);
	        }

	        // Discover & start all the receives
	        dataCommunicator->discoverRecvs();
	        dataCommunicator->startAllRecvs();

	        // Receive the target ghost vertices

	        int nCompletedRecvs = 0;
	        while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
	            int source_rank = dataCommunicator->waitAnyRecv();
	            bitpit::RecvBuffer &buffer = dataCommunicator->getRecvBuffer(source_rank);

	            // Receive ghost vertices
	            std::size_t vertexCount;
	            buffer >> vertexCount;
	            for (std::size_t ivertex = 0; ivertex < vertexCount; ivertex++){
	                bitpit::Vertex ghostVertex;
	                buffer >> ghostVertex;
	                long oldGhostVertexId = ghostVertex.getId();
	                // ADd vertex if not already in the tree
	                long ghostVertexId = oldGhostVertexId;
	                std::pair<int, long> ghost_label(source_rank, ghostVertexId);
	                std::pair<int, long> old_ghost_label(ghost_label);

	                // Check if already an insert ghost vertex
	                if (ghost_tree.exist(&ghostVertex, ghost_label) == -1){

	                    // Check if coincident qith a local vertex
	                    if ((mesh->getKdTree()->exist(&ghostVertex, ghostVertexId)) == -1){

	                        // To insert ghost vertex
	                        unique_ghosts[source_rank].push_back(ghostVertex);

	                        // Insert in the tree because the received vertices can be duplicated between partitions
	                        ghost_tree.insert(&unique_ghosts[source_rank].back(), ghost_label);

	                    }

	                    // Insert item in map oldid->newid for ghost vertices coincident with a vertex
	                    // already in the local tree
	                    if (ghostVertexId != oldGhostVertexId){
	                        oldToNewVertexId[source_rank].insert({{oldGhostVertexId, ghostVertexId}});
	                    }

	                }

	                // Insert item in map oldid->newid for ghost vertices coincident with a vertex
	                // already in the tree of ghosts [received by a previous partition])
	                if (ghost_label != old_ghost_label){
	                    oldToRenumberedGhostVertexId[source_rank].insert({{old_ghost_label.second, ghost_label}});
	                }

	            }
	            ++nCompletedRecvs;
	        }

	        // Wait for the sends to finish
	        dataCommunicator->waitAllSends();
	    }

	    // Insert unique vertices in mesh
	    for (auto & entry : unique_ghosts){
	        int source_rank = entry.first;
	        for (auto & ghostVertex : entry.second){
	            long oldGhostVertexId = ghostVertex.getId();
	            long ghostVertexId = bitpit::Vertex::NULL_ID;
	            ghostVertex.setId(bitpit::Vertex::NULL_ID);
	            ghostVertexId = mesh->addVertex(ghostVertex.getCoords(), bitpit::Vertex::NULL_ID);
	            // Insert new vertex id in map old->new to update connectivity
	            if (ghostVertexId != oldGhostVertexId){
	                oldToNewVertexId[source_rank].insert({{oldGhostVertexId, ghostVertexId}});
	            }
	        }
	    }

	    // Insert in the map the id of the ghost vertices dependent
	    // from the insertion of the ghost vertices received from other partitions
	    for (auto entry : oldToRenumberedGhostVertexId){
	        int source_rank = entry.first;
	        for (auto tuple : entry.second){
	            long old_id = tuple.first;
	            int linked_source_rank = tuple.second.first;
	            long linked_old_id = tuple.second.second;
	            oldToNewVertexId[source_rank].insert({{old_id, oldToNewVertexId[linked_source_rank][linked_old_id]}});
	        }
	    }

	    {
	        // Send cells
	        std::unique_ptr<bitpit::DataCommunicator> dataCommunicator;
	        dataCommunicator = std::unique_ptr<bitpit::DataCommunicator>(new bitpit::DataCommunicator(getCommunicator()));

	        for (auto & entry : ghostCellExchangeSources){

	            int target_rank = entry.first;
	            std::unordered_set<long> & cellList = entry.second;

	            // Prepare buffer for source cells
	            size_t cellBufferSize = 0;
	            for (long cellId : cellList){
	                cellBufferSize += mesh->getPatch()->getCell(cellId).getBinarySize();
	            }
	            dataCommunicator->setSend(target_rank, sizeof(size_t) + cellBufferSize);
	            bitpit::SendBuffer &buffer = dataCommunicator->getSendBuffer(target_rank);

	            // Fill the buffer with the source cells
	            buffer << cellList.size();
	            for (long cellId : cellList){
	                buffer << mesh->getPatch()->getCell(cellId);
	            }
	            dataCommunicator->startSend(target_rank);

	        }

	        // Discover & start all the receives
	        dataCommunicator->discoverRecvs();
	        dataCommunicator->startAllRecvs();

	        // Receive the target ghost cells
	        int nCompletedRecvs = 0;
	        while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
	            int source_rank = dataCommunicator->waitAnyRecv();
	            bitpit::RecvBuffer &buffer = dataCommunicator->getRecvBuffer(source_rank);

	            // Receive ghost cells add to mesh partition
	            std::size_t cellCount;
	            buffer >> cellCount;
	            for (std::size_t icell = 0; icell < cellCount; icell++){
	                bitpit::Cell ghostCell;
	                buffer >> ghostCell;
	                ghostCell.renumberVertices(oldToNewVertexId[source_rank]);
	                ghostCell.setId(bitpit::Cell::NULL_ID);
	                long ghostId = mesh->addCell(ghostCell, int(source_rank));

	            }
	            ++nCompletedRecvs;
	        }
	        // Wait for the sends to finish
	        dataCommunicator->waitAllSends();

	    }

	} // end if multi-processors

    // Update mesh
    mesh->update();

#endif

    mesh->updateAdjacencies();
    mesh->updateInterfaces();
    mesh->update();

    // TODO bitpit sort is bugged
    // Sort cells and vertices
//    mesh->getPatch()->sortVertices();
//    mesh->getPatch()->sortCells();

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
            bitpit::ConstProxyVector<long> vconn = bitInterfaces[bitFaceList[j]].getVertexIds();
            std::size_t vconnsize = vconn.size();
            std::vector<long> vListBIT(vconn.begin(), vconn.end());
            std::sort(vListBIT.begin(), vListBIT.end());

            if(std::equal(vListBIT.begin(), vListBIT.end(), vListOF.begin()) ){
                iBIT = bitFaceList[j];
            }else{
                ++j;
            }
        }
        m_OFbitpitmapfaces.insert(std::make_pair(long(iOF),iBIT) );
    }


	// From MimmoObject utilities, create the raw boundary mesh, storing it in m_boundary internal member.
	// PID will be every where 0. Need to be compiled to align with patch division of foamBoundaryMesh.
	// Every cell of the boundary mesh will have the same id of the border interfaces of the bulk.
    m_boundary = mesh->extractBoundaryMesh();

	//Once OFoam faces and bitpit Interfaces link is set,
	//Extract boundary patch info from foamBoundary and
	// get boundary patch division automatically pidding the class boundary mesh.

	const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
	long startIndex;
	long endIndex;

	forAll(foamBMesh, iBoundary){
		long PID = long(iBoundary+1);
		startIndex = foamBMesh[iBoundary].patch().start();
		endIndex = startIndex + long(foamBMesh[iBoundary].patch().size());
        for(long ind=startIndex; ind<endIndex; ++ind){
            m_boundary->setPIDCell(m_OFbitpitmapfaces[ind], PID);
		}
        m_boundary->setPIDName(PID,foamBMesh[iBoundary].name().c_str());
	}
	m_boundary->resyncPID();

	// Update bounfary mesh
	m_boundary->update();

	// Destroy interfaces to save memory
	m_geometry->destroyInterfaces();

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

	dvecarr3E points;

#if MIMMO_ENABLE_MPI
	// If MPI enable and multi-process run use only vertices of internal cells
	if (getProcessorCount() > 1){
        // Recover vertex ids of internal cells only
	    std::vector<long> vectorVertexIds = getGeometry()->getVertexFromCellList(getGeometry()->getCellsIds(true));
        // In order to maintain the same order of the original mesh
	    // put the find vertex ids in a set
	    points.reserve(vectorVertexIds.size());
	    // Run over vertices of internal cells and put coordinates of the vertex in the set
	    std::set<long> orderedVertexIds;
	    for (long vertexId : vectorVertexIds){
	        orderedVertexIds.insert(vertexId);
	    }
        for (long vertexId : orderedVertexIds){
            points.emplace_back(getGeometry()->getPatch()->getVertex(vertexId).getCoords());
        }
	}
	else
#endif
	{
	    points = getGeometry()->getVerticesCoords();
	}

    bool check = foamUtilsNative::writePointsOnCase(m_path.c_str(), points, m_overwrite);

    // In foamUtilsNative::writePointsOnCase the number of points of the initialized case are checked
    // vs the number of points of the input mesh.
    // It is the only check performed, the coherence between the order of the points is charged to the user
    // No coherence check on other structures and info of the mesh are done
    if(!check){
        (*m_log)<<"ERROR: "<<m_name<<" OpenFOAM mesh points coherence check failed during writing."<<std::endl;
        return false;
    }

    return true;

}

/*
 * ========================================================================================================
 */

/*!
 * Default constructor of IOOFOAMScalarField.
 \param[in] writemode if true, put the class in write mode, otherwise in read mode.
 */
IOOFOAMScalarField::IOOFOAMScalarField(bool writemode):IOOFOAM_Kernel(writemode){
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

    std::string fallback_name2 = "0";
	std::string input2 = rootXML.get("IOMode", fallback_name2);
	input2 = bitpit::utils::string::trim(input2);
    {
        std::stringstream ss(input2);
        ss>>m_write;
    }

	if(input == "mimmo.IOOFOAMScalarField"){
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
	bool built = true;

    //depending if the field to be read/written is bulk or boundary, you need at least bulk or boundary connected.
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::setGeometry, M_GEOMOFOAM, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::setGeometry, M_GEOM, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOMOFOAM2, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOM2, true,1));
    built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::setFacesMap, M_UMAPIDS));

    // creating output ports
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &BaseManipulation::getGeometry, M_GEOMOFOAM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &BaseManipulation::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOM2));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOMOFOAM2));
    built = (built && createPortOut<std::unordered_map<long,long>, IOOFOAMScalarField>(this, &IOOFOAM_Kernel::getFacesMap, M_UMAPIDS));
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
    if(!m_write){
        check = read();
    }else{
        check = write();
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
	}

	m_field.setGeometry(getGeometry());
	//Data on cells
	m_field.setDataLocation(2);

    //TODO ADD CELL TO POINT INTERPOLATION
    //TODO ADD COMMUNICATION ON GHOSTS

	// Read boundary field
    dmpvector1D boundaryFieldOnFace;

    //One field stored.
    if ( getBoundaryGeometry() != nullptr ){
        const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
        long startIndex;
        std::unordered_set<long> pids = getBoundaryGeometry()->getPIDTypeList();
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

#if MIMMO_ENABLE_MPI
    // Communicate ghost data
    m_boundaryField.communicateData();
#endif

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

    (*m_log)<<"WARNING: "<<m_name<<" Writing OF field is not available yet. Did nothing."<<std::endl;

    return true;
}



/*
 * ========================================================================================================
 */
 /*!
  * Default constructor of IOOFOAMVectorField.
  \param[in] writemode if true, put the class in write mode, otherwise in read mode.
 */
 IOOFOAMVectorField::IOOFOAMVectorField(bool writemode):IOOFOAM_Kernel(writemode){
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

    std::string fallback_name2 = "0";
	std::string input2 = rootXML.get("IOMode", fallback_name2);
	input2 = bitpit::utils::string::trim(input2);
    {
        std::stringstream ss(input2);
        ss>>m_write;
    }

 	if(input == "mimmo.IOOFOAMVectorField"){
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
    bool built = true;

    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::setGeometry, M_GEOMOFOAM, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::setGeometry, M_GEOM, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOMOFOAM2, true,1));
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::setBoundaryGeometry, M_GEOM2, true,1));
    built = (built && createPortIn<std::unordered_map<long,long>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::setFacesMap, M_UMAPIDS));

    // creating output ports
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &BaseManipulation::getGeometry, M_GEOMOFOAM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &BaseManipulation::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOMOFOAM2));
    built = (built && createPortOut<MimmoSharedPointer<MimmoObject>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::getBoundaryGeometry, M_GEOM2));
    built = (built && createPortOut<std::unordered_map<long,long>, IOOFOAMVectorField>(this, &IOOFOAM_Kernel::getFacesMap, M_UMAPIDS));
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
 	if(!m_write){
        check = read();
    }else{
        check = write();
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
	}

    m_field.setGeometry(getGeometry());
    m_field.setDataLocation(2);

    //TODO ADD CELL TO POINT INTERPOLATION
    //TODO ADD COMMUNICATION ON GHOSTS

    // Read boundary field
    dmpvecarr3E boundaryFieldOnFace;

    //One field stored.
    if ( getBoundaryGeometry() != nullptr ){

        const Foam::fvBoundaryMesh &foamBMesh = foamMesh->boundary();
        long startIndex;

        std::unordered_set<long> pids = getBoundaryGeometry()->getPIDTypeList();
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

#if MIMMO_ENABLE_MPI
    // Communicate ghost data
    m_boundaryField.communicateData();
#endif

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

	(*m_log)<<"WARNING: "<<m_name<<" Writing OF field is not available yet. Did nothing."<<std::endl;

	return true;
}

}
