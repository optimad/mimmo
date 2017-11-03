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
 */
IOOFOAM::IOOFOAM(){
    m_name          = "mimmo.IOOFOAM";
    setDefaults();
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOOFOAM::IOOFOAM(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.IOOFOAM";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.IOOFOAM"){
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
IOOFOAM::IOOFOAM(const IOOFOAM & other):BaseManipulation(other){
    m_read = other.m_read;
    m_readPath = other.m_readPath;

    m_write = other.m_write;
    m_writePath = other.m_writePath;
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
   std::swap(m_read, x.m_read);
   std::swap(m_readPath, x.m_readPath);

   std::swap(m_write, x.m_write);
   std::swap(m_writePath, x.m_writePath);
   std::swap(m_volmesh, x.m_volmesh);

   BaseManipulation::swap(x);
};

/*!
 * Default values for IOOFOAM.
 */
void
IOOFOAM::setDefaults(){

    m_read       = false;
    m_readPath   = ".";
    m_write      = false;
    m_writePath  = "./";
    
    m_OFE_supp["hex"]   = bitpit::ElementType::HEXAHEDRON;
    m_OFE_supp["tet"]   = bitpit::ElementType::TETRA;
    m_OFE_supp["prism"] = bitpit::ElementType::WEDGE;
    m_OFE_supp["pyr"]   = bitpit::ElementType::PYRAMID;
}

/*! It builds the input/output ports of the object
 */
void
IOOFOAM::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoObject*, IOOFOAM>(this, &IOOFOAM::setGeometry, M_GEOM));
//     built = (built && createPortIn<MimmoObject*, IOOFOAM>(this, &IOOFOAM::setSurfaceBoundary, M_GEOM2));
//     built = (built && createPortIn<dmpvector1D, IOOFOAM>(this, &IOOFOAM::setField, M_SCALARFIELD));

    built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &IOOFOAM::getGeometry, M_GEOM));
//     built = (built && createPortOut<MimmoObject*, IOOFOAM>(this, &IOOFOAM::getSurfaceBoundary, M_GEOM2));
//     built = (built && createPortOut<dmpvector1D, IOOFOAM>(this, &IOOFOAM::getField, M_SCALARFIELD));
    m_arePortsBuilt = built;
};

/*!
 * It sets to true to read the mesh during execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
IOOFOAM::setRead(bool read){
    m_read = read;
}

/*!
 * It sets the name of directory to read the OpenFOAM mesh.
 * \param[in] dir mesh input directory.
 */
void
IOOFOAM::setReadDir(std::string dir){
    m_readPath = dir;
}

/*!
 * It sets to true to write the mesh during execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
IOOFOAM::setWrite(bool write){
    m_write = write;
}

/*!
 * It sets the name of directory to write the OpenFOAM mesh.
 * \param[in] dir mesh output directory
 */
void
IOOFOAM::setWriteDir(std::string dir){
    m_writePath = dir;
}


/*!
 * Set current geometry to a Volume Mesh MimmoObject.
 * Method meant for writing class mode only.
 * Any prevoius mesh internally allocated will be destroyed.
 * \param[in] geom pointer to external volume mesh MimmoObject.
 */
void
IOOFOAM::setGeometry(MimmoObject * geom){
    if(geom == NULL || m_read)   return;
    if(geom->isEmpty())  return;
    if(geom->getType() != 2) return;

    BaseManipulation::setGeometry(geom);
}


/*!
 * Return current pointer to geometry.If read mode is active return local read volumetric mesh, else
 * otherwise return pointed externally geometry
 * \return pointer to linked volume mesh
 */
MimmoObject*
IOOFOAM::getGeometry(){
    if(m_read) return m_volmesh.get();
    else return BaseManipulation::getGeometry();
}

/*!
 * Clone actual internal mesh of the class, in an independent data structure.
 * \return unique pointer of the cloned mesh.
 */
std::unique_ptr<MimmoObject>
IOOFOAM::cloneInternalMesh(){
    return std::move(m_volmesh->clone());
}

/*!
 * It reads the OpenFOAM mesh from input file and store in a volume mesh INTERNAL MimmoObject.
 * \return false if errors occured during the reading.
 */
bool
IOOFOAM::read(){

    Foam::Time *foamRunTime = 0;
    Foam::fvMesh *foamMesh = 0;

    foamUtilsNative::initializeCase(m_readPath.c_str(), &foamRunTime, &foamMesh);

   
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

    m_volmesh = std::move(std::unique_ptr<MimmoObject>(new MimmoObject(2, mesh)));
    
    
    //boundary patch info.

    return true;

}

/*!
 * It writes the OpenFOAM mesh to an output file form an externally linked volume mesh MimmoObject.
 * \return false if errors occured during the writing.
 */
bool
IOOFOAM::write(){
    //do nothing for now
    return true;
}


/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
IOOFOAM::execute(){
    bool check = true;
    if (m_read) check = read();
    if (!check){
//         if (m_stopat == SHRT_MAX){
//             (*m_log) << m_name << " error: file not found : "<< m_rfilenameV << std::endl;
//             (*m_log) << " " << std::endl;
//             throw std::runtime_error (m_name + ": file not found : " + m_rfilenameV);
//         }
//         (*m_log) << m_name << " error: file not found : "<< m_rfilenameS[m_stopat] << std::endl;
//         (*m_log) << " " << std::endl;
        throw std::runtime_error (m_name + ": an error occured while reading from files");
    }
    if (m_write) check = write();
    if (!check){
//         if (m_stopat == 2){
//             (*m_log) << m_name << " error: write not done : surface and volume geometry not linked " << std::endl;
//             (*m_log) << " " << std::endl;
//             throw std::runtime_error (m_name + ": write not done : surface and volume geometry not linked ");
//         }
        throw std::runtime_error (m_name + ": an error occured while writing on files");
        
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
    
    if(slotXML.hasOption("ReadFlag")){
        input = slotXML.get("ReadFlag");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setRead(value);
    };

    if(slotXML.hasOption("ReadDir")){
        input = slotXML.get("ReadDir");
        input = bitpit::utils::string::trim(input);
        if(input.empty())   input = "./";
        setReadDir(input);
    };

    if(slotXML.hasOption("WriteFlag")){
        input = slotXML.get("WriteFlag");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        setWrite(value);
    };

    if(slotXML.hasOption("WriteDir")){
        input = slotXML.get("WriteDir");
        input = bitpit::utils::string::trim(input);
        if(input.empty())   input = "./";
        setWriteDir(input);
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

    std::string output;

    output = std::to_string(m_read);
    slotXML.set("ReadFlag", std::to_string(m_read));
    slotXML.set("ReadDir", m_readPath);
    slotXML.set("WriteFlag", std::to_string(m_write));
    slotXML.set("WriteDir", m_writePath);
};

/*!
 * Reorder OpenFoam vertex-cell connectivity of an elementary shape in 
 * a suitable one for corrispondent bitpit elementary shape. The reordering
 * preserve the local enumeration of faces which is the same in both OpenFOAM and bitpit.
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
    
}
