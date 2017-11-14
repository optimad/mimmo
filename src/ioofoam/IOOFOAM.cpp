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

    setDefaults();

    if(input == "mimmo.IOOFOAM" && maybeIOOFMode){
        m_type = maybeIOOFMode->_to_integral();
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

    if(!slotXML.hasOption("IOMode")){
        warningXML(m_log, m_name);
        return;
    }

    input = slotXML.get("IOMode");
    if(input.empty()){
        warningXML(m_log, m_name);
        return;
    }
    better_enums::optional<IOOFMode> maybeIOMode;
    std::stringstream ss(bitpit::utils::string::trim(input));
    maybeIOMode = IOOFMode::_from_string_nothrow(input.c_str());

    if(!maybeIOMode) {
        warningXML(m_log, m_name);
        return;
    }
    if( m_type == maybeIOMode->_to_integral()){
        warningXML(m_log, m_name);
        return;
    }

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
IOOFOAM::mapEleVConnectivity(const livector1D & FoamConn, const bitpit::ElementType & eltype){
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
                temp.clear();
                Foam::label iFace = cells[iC][locC];
                temp.push_back((long)faces[iFace].size());
                forAll(faces[iFace], locF){
                    temp.push_back((long)faces[iFace][locF]);
                }
                conn.insert(conn.end(), temp.begin(), temp.end());
                
                if(iFace >= sizeNeighbours) continue;
                if(iC == faceOwner[iFace]){
                    adjacency[int(locC)][0] = faceNeighbour[iFace];
                }else{
                    adjacency[int(locC)][0] = faceOwner[iFace];
                }
            }
        }
        bitpit::PatchKernel::CellIterator it = mesh->addCell(eltype, true, conn, iDC);
        it->setPID(int(PID));
        it->setAdjacencies(adjacency);
    }
    //build as bitpit the Interfaces -> adjacency is automatically recoverd from openFoam info.
    mesh->buildInterfaces();
    
    //create OpenFoam faces - bitpit Interfaces local map to detect boundaries
    std::unordered_map<long,long> mapBitOF_faces;

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
            mapBitOF_faces[ofFaceList[i]] = bitFaceList[i];
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
            m_boundary->setPIDCell(mapBitOF_faces[ind], PID);
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


}
