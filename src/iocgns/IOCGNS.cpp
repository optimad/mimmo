/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#include "IOCGNS.hpp"
#include <cgnslib.h>
#include <unordered_map>

using namespace std;
using namespace bitpit;

namespace mimmo{


/*!Default constructor of IOCGNS.
 */
IOCGNS::IOCGNS(bool read){
    setDefaults();
    m_read  = read;
    m_write = !read;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOCGNS::IOCGNS(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.IOCGNS";
    setDefaults();

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::trim(input);
    if(input == "mimmo.IOCGNS"){
        absorbSectionXML(rootXML);
    }else{
        std::cout<<"Warning in custom xml mimmo::IOCGNS constructor. No valid xml data found"<<std::endl;
    };
}

/*!Default destructor of IOCGNS.
 */
IOCGNS::~IOCGNS(){};

/*!Copy constructor of IOCGNS.
 */
IOCGNS::IOCGNS(const IOCGNS & other):BaseManipulation(){
    *this = other;
};

/*!Assignement operator of IOCGNS.
 */
IOCGNS & IOCGNS::operator=(const IOCGNS & other){

    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_read = other.m_read;
    m_rfilename = other.m_rfilename;
    m_write = other.m_write;
    m_wfilename = other.m_wfilename;
    m_rdir = other.m_rdir;
    m_wdir = other.m_wdir;
    m_surfmesh_not = other.m_surfmesh_not;
    m_storedInfo = new InfoCGNS((*other.m_storedInfo));
    m_storedBC = new BCCGNS((*other.m_storedBC));
    return *this;
};

/*!Default values for IOCGNS.
 */
void
IOCGNS::setDefaults(){

    m_name      = "mimmo.IOCGNS";
    m_read      = false;
    m_rfilename = "mimmoVolCGNS";
    m_write     = false;
    m_wfilename = "mimmoVolCGNS";
    m_rdir      = "./";
    m_wdir      = "./";
    m_surfmesh_not = NULL;
    m_storedInfo = new InfoCGNS;
    m_storedBC  = new BCCGNS;

    //Fill converters
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::TRIANGLE] = CG_ElementType_t::CG_TRI_3;
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::QUAD] = CG_ElementType_t::CG_QUAD_4;
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::TETRA] = CG_ElementType_t::CG_TETRA_4;
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::PYRAMID] = CG_ElementType_t::CG_PYRA_5;
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::WEDGE] = CG_ElementType_t::CG_PENTA_6;
    m_storedInfo->mcg_typeconverter[bitpit::ElementInfo::Type::HEXAHEDRON] = CG_ElementType_t::CG_HEXA_8;

    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_TRI_3] = "Elem_tri";
    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_QUAD_4] = "Elem_quad";
    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_TETRA_4] = "Elem_tetra";
    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_PYRA_5] = "Elem_pyra";
    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_PENTA_6] = "Elem_prism";
    m_storedInfo->mcg_typetostring[CG_ElementType_t::CG_HEXA_8] = "Elem_hexa";

}

/*! It builds the input/output ports of the object
 */
void
IOCGNS::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, IOCGNS>(this, &IOCGNS::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<MimmoObject*, IOCGNS>(this, &IOCGNS::setSurfaceBoundary, PortType::M_GEOM2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortIn<BCCGNS*, IOCGNS>(this, &IOCGNS::setBoundaryConditions, PortType::M_BCCGNS, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BCCGNS_));

    built = (built && createPortOut<MimmoObject*, IOCGNS>(this, &IOCGNS::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<MimmoObject*, IOCGNS>(this, &IOCGNS::getSurfaceBoundary, PortType::M_GEOM2, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<BCCGNS*, IOCGNS>(this, &IOCGNS::getBoundaryConditions, PortType::M_BCCGNS, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BCCGNS_));
    m_arePortsBuilt = built;
};

/*!
 * Return all the surface bounding the current volume mesh.
 * During reading mode returns info contained in cgns, marking with PID all the possible boundary patches, 
 * while in writing mode refers to the actual object pointed by the User.
 */
MimmoObject*
IOCGNS::getSurfaceBoundary(){
    if(m_read)	return m_surfmesh.get();
    //    if(m_surfmesh != NULL) return m_surfmesh.get();
    else	return m_surfmesh_not;
}

/*!
 * Return current pointer to geometry.If read mode return local read volumetric mesh, else
 * if in write mode return pointed externally geometry
 */
MimmoObject*
IOCGNS::getGeometry(){
    if(m_read) return m_volmesh.get();
    //    if(m_volmesh != NULL) return m_volmesh.get();
    else return BaseManipulation::getGeometry();
}

/*!
 * Return the Info object (class BCCGNS) for boundary conditions of surface patches
 * as PIDs associated to Surface boundary mesh.
 * Meaningful only in reading mode.
 * \return Pointer to BCCGNS object with boundary conditions information read from file.
 */
BCCGNS*
IOCGNS::getBoundaryConditions(){
    return m_storedBC;
};

/*!
 * Return true if current class is in reading mode
 */
bool
IOCGNS::isReadingMode(){
    return m_read;
}

/*!
 * Return true if current class is in writing mode
 */
bool
IOCGNS::isWritingMode(){
    return m_write;
}

/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
IOCGNS::setRead(bool read){
    m_read = read;
    m_write = !read;
}

/*!It sets the name of directory to read the geometry.
 * If writing mode is on, set the path to read the template cgns file
 * \param[in] dir Name of directory.
 */
void
IOCGNS::setReadDir(string dir){
    m_rdir = dir;
}

/*!It sets the name of file to read the geometry.
 * If writing mode is on, set the filename where read the template cgns file
 * \param[in] filename Name of input file.
 */
void
IOCGNS::setReadFilename(string filename){
    m_rfilename = filename;
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
IOCGNS::setWrite(bool write){
    setRead(!write);
}

/*!It sets the name of directory to write the geometry.
 * \param[in] dir Name of directory.
 */
void
IOCGNS::setWriteDir(string dir){
    m_wdir = dir;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filename Name of output file.
 */
void
IOCGNS::setWriteFilename(string filename){
    m_wfilename = filename;
}

/*!
 * Set current geometry to an external volume mesh.
 */
void
IOCGNS::setGeometry(MimmoObject * geo){
    if(geo == NULL || m_read)	return;
    if(geo->isEmpty())	return;
    if(geo->getType() != 2)	return;

    BaseManipulation::setGeometry(geo);
    //m_volmesh.reset(nullptr);
}

/*!
 * Set boundary surface relative to the volume mesh.Option active only in writing mode.
 * Preexistent boundary surfaces read from file, when class is set in read mode,  will be erased.
 */
void			
IOCGNS::setSurfaceBoundary(MimmoObject* geosurf){
    if(geosurf == NULL || m_read)	return;
    if(geosurf->isEmpty())	return;
    if(geosurf->getType() != 1)	return;

    m_surfmesh_not = geosurf;
    //m_surfmesh.reset(nullptr);
};

/*!
 * Set the Info object (class BCCGNS) for boundary conditions of surface patches
 * as PIDs associated to Surface boundary mesh.
 * Meaningful only in writing mode.
 * \param[in] Pointer to BCCGNS object with boundary conditions information read from file.
 */
void
IOCGNS::setBoundaryConditions(BCCGNS* bccgns){
    delete m_storedBC;
    m_storedBC = NULL;
    m_storedBC = new BCCGNS(*bccgns);
};


/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
IOCGNS::execute(){
    bool check = true;
    if (m_read) {
        check = read();
    }
    if (!check){
        std::cout << "mimmo : ERROR : file corrupted or not found : "<< m_rfilename << std::endl;
        std::cout << " " << std::endl;
        exit(10);
    }
    if (m_write) {
        check = write();
    }
    if (!check){
        std::cout << "mimmo : ERROR : write not done : geometry not linked " << std::endl;
        std::cout << " " << std::endl;
        exit(11);
    }
}

/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists, doesn't hold a cgns unique base, zone, unstructured volume mesh.
 */
bool
IOCGNS::read(){

    std::string file = m_rdir+"/"+m_rfilename+".cgns";
    std::string error_string = "read CGNS grid: " + file;

    long nVertices;
    //     long nCells, nBoundVertices;
    int nCoords; //, nConnSections;
    std::array< std::vector<double>,3 > coords;
    std::vector<CG_ElementType_t> orderedConns;
    std::unordered_map<int, ivector1D> conns;
    std::unordered_map<int, ivector1D > bcLists;

    //Open cgns file
    int indexfile;
    if(cg_open(file.c_str(), CG_MODE_READ, &indexfile) != CG_OK){
        return false;
    }

    //Read number of bases
    int nbases;
    if(cg_nbases(indexfile, &nbases) != CG_OK){
        return false;
    }
    if(nbases > 1){
        std::cout << "More than one base found in grid file -> only the first zone will be read" << std::endl;
    }

    //Read name of basis and physical dimension
    char basename[33];
    int physdim, celldim;
    if(cg_base_read(indexfile, 1, basename, &celldim, &physdim) != CG_OK){
        return false;
    }
    if(celldim != 3 || physdim !=3){
        //Only volume mesh supported
        //TODO NO SURFACE MESH (ONLY SURFACE WITHOUT VOLUME)?
        return false;
    }

    //Read number of zone in basis, check if only one.
    int nzones;
    if(cg_nzones(indexfile, 1, &nzones)!=CG_OK){
        return false;
    }
    if(nzones > 1){
        //Only single volume zone supported for now
        return false;
    }

    //Read type mesh of zone
    CG_ZoneType_t zoneType;
    int index_dim;
    if(cg_zone_type(indexfile,1,1, &zoneType)!= CG_OK){
        return false;
    }
    if(cg_index_dim(indexfile,1,1, &index_dim)!= CG_OK){
        return false;
    }
    if(zoneType != CG_ZoneType_t::CG_Unstructured || index_dim != 1){
        //Only unstructured mesh supported for now (index_dim == 1 for unstructured)
        return false;
    }

    //Read size of zone (n nodes, n cells, n boundary nodes)
    std::vector<cgsize_t> sizeG(3);
    char zonename[33];
    if(cg_zone_read(indexfile,1,1, zonename, sizeG.data()) != CG_OK ){
        return false;
    }


    nVertices 		= sizeG[0];
    //     nCells 			= sizeG[1];
    //     nBoundVertices	= sizeG[2];

    //Read Vertices.
    if(cg_ncoords(indexfile,1,1, &nCoords)!= CG_OK){
        return false;
    }
    for(int i = 1; i <= nCoords; ++i){
        cgsize_t startIndex = 1;
        cgsize_t finishIndex = sizeG[0];
        CG_DataType_t datatype;
        char name[33];
        coords[i-1].resize(nVertices);
        if(cg_coord_info(indexfile,1,1,i, &datatype, name)!=CG_OK){
            return false;
        }
        if(cg_coord_read(indexfile,1,1,name, datatype, &startIndex, &finishIndex, coords[i-1].data() )!=CG_OK){
            return false;
        }
    }

    //Read connectivities.
    //They are read starting from 1, fortran style. When matching up w/ coords
    //vector positions remember to diminish conn value of 1.
    //Read number of sections of connectivity structure
    int nSections;
    if(cg_nsections(indexfile,1,1, &nSections)!= CG_OK){
        return false;
    }


    for(int sec = 1; sec <= nSections; ++sec){

        //Read section sec
        char elementname[33];
        CG_ElementType_t type;
        cgsize_t eBeg, eEnd, enBdry;
        int parent_flag;
        std::vector<cgsize_t>	connlocal;
        cgsize_t size;

        //Read elements name, type and range
        if(cg_section_read(indexfile,1,1, sec, elementname, &type, &eBeg, &eEnd,&enBdry, & parent_flag)!= CG_OK){
            return false;
        }

        //Read size of connectivity data
        if(cg_ElementDataSize(indexfile,1,1, sec,&size)!= CG_OK){
            return false;
        }
        //cgsize_t nelements = eEnd-eBeg + 1;

        connlocal.resize((size_t) size);
        int *ptr = NULL;
        if(cg_elements_read(indexfile,1,1,sec, connlocal.data(),ptr) !=CG_OK){
            return false;
        }

        conns[type].resize(connlocal.size());
        int count=0;
        for(auto &val: connlocal){
            conns[type][count] = (int)val;
            ++count;
        }

        orderedConns.push_back(type);
    }

    //explore superficial boundary conditions definition, for boundary surface extraction.
    int nBcs;
    if(cg_nbocos(indexfile,1,1, &nBcs)!= CG_OK){
        return false;
    }

    //up to now i'm not able to encoding all the possible variants of bc node structure.
    // i rely on the element list . Need to study better the cgns docs.
    // TODO ????
    for(int bc=1; bc<=nBcs; ++bc){

        char name[33];
        CG_BCType_t bocotype;
        CG_PointSetType_t ptset_type;
        std::vector<cgsize_t> nBCElements(2);
        int normalIndex;
        cgsize_t normalListSize;
        CG_DataType_t normalDataType;
        int ndataset;
        CG_GridLocation_t location;

        if(cg_boco_gridlocation_read(indexfile,1,1,bc, &location) != CG_OK){
            return false;
        }

        if(cg_boco_info(indexfile,1,1,bc, name, &bocotype, &ptset_type, nBCElements.data(),
                &normalIndex, &normalListSize, &normalDataType, &ndataset) != CG_OK){
            return false;
        }

        if(ptset_type == CG_PointSetType_t::CG_ElementList || ptset_type == CG_PointSetType_t::CG_ElementRange){

            cgsize_t dim;
            if(ptset_type == CG_PointSetType_t::CG_ElementList){
                dim = nBCElements[0];
            }else{
                dim = nBCElements[1] - nBCElements[0] +1;
            }

            std::vector<cgsize_t> localbc((size_t) dim);
            int *ptr = NULL;

            if(cg_boco_read(indexfile,1,1,bc, localbc.data(), ptr )!= CG_OK){
                return false;
            }

            bcLists[bc-1].resize(localbc.size());
            int count=0;
            for(auto &val: localbc){
                // WHY INT AND NOT LONG INT (IT IS A BITPIT PID) ??
                bcLists[bc-1][count] = (int)val;
                ++count;
            }

        }

        m_storedBC->mcg_pidtobc[bc] = bocotype;
        m_storedBC->mcg_pidtoname[bc] = name;

    }

    m_storedBC->mcg_pidtobc[0] = CG_BCTypeNull;
    m_storedBC->mcg_pidtoname[0] = "undefined";

    //Finish reading CGNS file
    cg_close(indexfile);

    //Reverse info in your grids.
    std::unique_ptr<MimmoObject> patchVol(new MimmoObject(2));
    std::unique_ptr<MimmoObject> patchBnd(new MimmoObject(1));
    long id;
    darray3E temp;

    //	Stock vertices in volume grid
    for(int i=0; i<nVertices; ++i){

        for(int j=0; j<3; ++j)	temp[j] = coords[j][i];
        id = i+1; //not C indexing, labeling coherent with connectivity.
        patchVol->addVertex(temp,id);
    }

    //debug only;
    //TODO ????

    id = 1;
    //Unpack 3D elements connectivities and store it in volume grid. Label same species elements w PID.
    for(auto & val: orderedConns){

        livector1D lConn;
        bitpit::ElementInfo::Type btype;
        short PID = 0;
        int size = conns[(int)val].size();

        switch(val){
        case CGNS_ENUMV(MIXED):

		                        unpack3DElementsMixedConns(patchVol.get(), patchBnd.get(), conns[(int)val], id);
        break;

        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):

        btype = bitpit::ElementInfo::Type::TETRA;
        lConn.resize(4);
        for(int i=0; i<size; i+=4){
            for(int j=0; j<4; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchVol->addConnectedCell(lConn, btype, id );
            patchVol->setPIDCell(id,PID);
            id++;
        }
        break;

        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_14):

        btype = bitpit::ElementInfo::Type::PYRAMID;
        lConn.resize(5);
        for(int i=0; i<size; i+=5){
            for(int j=0; j<5; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchVol->addConnectedCell(lConn, btype, id );
            patchVol->setPIDCell(id,PID);
            id++;
        }
        break;

        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):

        btype = bitpit::ElementInfo::Type::WEDGE;
        lConn.resize(6);
        for(int i=0; i<size; i+=6){
            for(int j=0; j<6; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchVol->addConnectedCell(lConn, btype, id );
            patchVol->setPIDCell(id,PID);
            id++;
        }

        break;

        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):

        btype = bitpit::ElementInfo::Type::HEXAHEDRON;
        lConn.resize(8);
        for(int i=0; i<size; i+=8){
            for(int j=0; j<8; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchVol->addConnectedCell(lConn, btype, id );
            patchVol->setPIDCell(id,PID);
            id++;
        }
        break;

        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):

        btype = bitpit::ElementInfo::Type::TRIANGLE;
        lConn.resize(3);

        for(int i=0; i<size; i+=3){
            for(int j=0; j<3; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchBnd->addConnectedCell(lConn, btype, id );
            patchBnd->setPIDCell(id,PID);
            id++;
        }
        break;

        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):

        btype = bitpit::ElementInfo::Type::QUAD;
        lConn.resize(4);
        for(int i=0; i<size; i+=4){
            for(int j=0; j<4; ++j){
                lConn[j] = conns[(int)val][i+j];
            }
            patchBnd->addConnectedCell(lConn, btype, id );
            patchBnd->setPIDCell(id,PID);
            id++;
        }
        break;

        default:
            //do nothing
            break;
        }
    }

    //divide patchBnd in subpatch if any bc is present.
    if(bcLists.size() > 0)	{
        for(auto & sel: bcLists){
            if(sel.second.size() > 0){
                for(auto &ind : sel.second){
                    id = ind;
                    short PID = sel.first + 1;
                    patchBnd->setPIDCell(id,PID);
                }
            }
        }//end for

    }


    //adding vstand alone vertices to boundary patches and release all structures
    {
        livector2D tempConn = patchBnd->getConnectivity();
        std::set<long> ordIndex;

        for(auto & val: tempConn){
            ordIndex.insert(val.begin(), val.end());
        }

        for(auto & val: ordIndex){
            temp = patchVol->getVertexCoords(val);
            patchBnd->addVertex(temp,val);
        }
    }

    //release the meshes
    m_volmesh = std::move(patchVol);
    m_surfmesh = std::move(patchBnd);

    return true;
}

/*!It writes the mesh geometry on output .cgns file.
 * If boundary surface, PID subdivided, is available, write all patches as CGNS Wall boundary condition .
 *\return False if volume geometry at least is not linked .
 */
bool
IOCGNS::write(){

    /* Clear old Info if stored. */
    m_storedInfo->mcg_number.clear();

    std::string file = m_wdir+"/"+m_wfilename+".cgns";
    std::string error_string = "write CGNS grid: " + file;

    long nVertices, nCells, nBoundVertices;
    //int nCoords;
    std::array< std::vector<double>,3 > coords;

    /* Fill this structures (surface boundary connectivity has
     * to be filled with same vertex indices of volume mesh) */
    MimmoObject * vol = getGeometry();
    MimmoObject * bnd = getSurfaceBoundary();

    if( vol == NULL && bnd == NULL ) return false;
    if( vol->isEmpty() && bnd->isEmpty() ) return false;

    /* Verify if a surface boundary mesh exists. */
    //bool flagBnd = (bnd != NULL);

    /* Check number of volume cells and global vertices. */
    nVertices = vol->getNVertex();
    nCells = vol->getNCells();
    nBoundVertices=0;

    //fill coordinates;
    coords[0].resize(nVertices);
    coords[1].resize(nVertices);
    coords[2].resize(nVertices);
    {
        dvecarr3E vert = vol->getVertexCoords();
        int count=0;
        for(auto &val:vert){
            for(int j=0; j<3; ++j)	coords[j][count] = val[j];
            ++count;
        }
    }

    //Write Base and Zone Info
    int indexfile;
    if(cg_open(file.c_str(), CG_MODE_WRITE, &indexfile) != CG_OK){
        exit(11);
    }

    int baseindex = 1;
    char basename[33] = "Base";
    int physdim=3, celldim=3;
    if(cg_base_write(indexfile,basename, celldim, physdim, &baseindex) != CG_OK){
        exit(11);
    }

    int zoneindex=1;
    CG_ZoneType_t zoneType =CG_ZoneType_t::CG_Unstructured ;
    //int index_dim = 1;
    char zonename[33] = "Zone  1";
    std::vector<cgsize_t> sizeG(3);
    sizeG[0] = nVertices;
    sizeG[1] = nCells;
    sizeG[2] = nBoundVertices;

    if(cg_zone_write(indexfile,baseindex, zonename, sizeG.data(), zoneType, &zoneindex) != CG_OK ){
        exit(11);
    }

    //writing vertices.
    svector1D names(3, "CoordinateX");
    names[1] = "CoordinateY";
    names[2] = "CoordinateZ";

    CG_DataType_t datatype= CG_DataType_t::CG_RealDouble;

    for(int i=1; i<=3; ++i){
        if(cg_coord_write(indexfile,baseindex,zoneindex, datatype, names[i-1].data(), coords[i-1].data(), &i )!=CG_OK){
            exit(11);
        }
    }

    //Recover CGNS Info from volume and surface mesh
    recoverCGNSInfo();

    /* Write volume elements */
    int sec = 0;
    cgsize_t eBeg = 1, eEnd;
    int nbdry = 0;
    for(auto &connMap : m_storedInfo->mcg_typetoconn){
        sec++;
        if(connMap.second.size() == 0) continue;

        CG_ElementType_t type = (connMap.first);

        eEnd = eBeg + m_storedInfo->mcg_number[type] - 1;

        if(cg_section_write(indexfile,baseindex,zoneindex,m_storedInfo->mcg_typetostring[type].data(), type,eBeg,eEnd, nbdry, connMap.second.data(), &sec) !=CG_OK){
            exit(11);
        }
        eBeg = eEnd+1;
    }

    /* Write surface elements */
    for(auto &connMap : m_storedInfo->mcg_bndtypetoconn){
        sec++;
        if(connMap.second.size() == 0) continue;

        CG_ElementType_t type = (connMap.first);

        eEnd = eBeg + m_storedInfo->mcg_bndnumber[type] - 1;

        if(cg_section_write(indexfile,baseindex,zoneindex,m_storedInfo->mcg_typetostring[type].data(), type,eBeg,eEnd, nbdry, connMap.second.data(), &sec) !=CG_OK){
            exit(11);
        }
        eBeg = eEnd+1;
    }

    /* Write boundary conditions if thay are present.
     */
    int bcid;
    for(auto & bc: m_storedBC->mcg_bndpidtoindex){

        /* Default tag for boundary conditions = Wall.
         * Default name for boundary conditions = wall_PID.
         */
        CG_BCType_t bocotype = CG_BCType_t::CG_BCWall;
        CG_PointSetType_t ptset_type = CG_ElementList;


        cgsize_t nelements = bc.second.size();
        int bcPID = bc.first;
        if (m_storedBC->mcg_pidtobc.count(bcPID)) bocotype = m_storedBC->mcg_pidtobc[bcPID];

        std::string bcname = "wall_"+std::to_string(bcPID);
        if (m_storedBC->mcg_pidtoname.count(bcPID)) bcname = m_storedBC->mcg_pidtoname[bcPID];

        if(cg_boco_write(indexfile, baseindex, zoneindex, bcname.data(),
                bocotype, ptset_type, nelements, bc.second.data(), &bcid )!= CG_OK){
            exit(11);
        }

    }

    /*
     * Default location = faceCenter ---> deactivated
    for (int ibc = 1; ibc <= bcid; ibc++){
        cg_goto(indexfile,baseindex,"Zone_t",1,"ZoneBC_t",1,"BC_t",ibc,"end");
        cg_gridlocation_write(CG_CellCenter);

    }
     */

    /* Finish writing CGNS file */
    cg_close(indexfile);
    return true;
}

/*!
 * Extract mixed connectivity of 3D elements volume mesh//surface boundary mesh in two separated objects, given an array 
 * of mixed connectivity of CGNS reader
 * \param[in,out]	patchVol pointer to Volume MimmoObject handler
 * \param[in,out]	patchSurf pointer to Surface MimmoObject handler
 * \param[in,out]	conn	list of vertex index connectivity of CGNS_ENUMV(MIXED)type
 * \param[in,out]	startID starting ID to labeling mixed cells. Returns incremented of found mixed 3D cells	
 */
void
IOCGNS::unpack3DElementsMixedConns(MimmoObject * patchVol, MimmoObject* patchSurf, ivector1D & conn, long & startId){

    livector1D lConn;
    bitpit::ElementInfo::Type btype;
    short PID=0;
    long id = startId;
    ivector1D::iterator it=conn.begin(), itE=conn.end();

    while(it !=itE){

        CGNS_ENUMT(ElementType_t) et = static_cast<CGNS_ENUMT(ElementType_t)>(*it);
        ++it;

        switch(et){
        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):

        btype = bitpit::ElementInfo::Type::TETRA;
        lConn.resize(4);
        for(int i=0; i<4; i++){
            lConn[i] = *it;
            ++it;
        }

        patchVol->addConnectedCell(lConn, btype, id );
        patchVol->setPIDCell(id,PID);
        id++;
        break;

        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_14):

        btype = bitpit::ElementInfo::Type::PYRAMID;
        lConn.resize(5);
        for(int i=0; i<5; i++){
            lConn[i] = *it;
            ++it;
        }

        patchVol->addConnectedCell(lConn, btype, id );
        patchVol->setPIDCell(id,PID);
        id++;
        break;

        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):

        btype = bitpit::ElementInfo::Type::WEDGE;
        lConn.resize(6);
        for(int i=0; i<6; i++){
            lConn[i] = *it;
            ++it;
        }

        patchVol->addConnectedCell(lConn, btype, id );
        patchVol->setPIDCell(id,PID);
        id++;
        break;

        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):

        btype = bitpit::ElementInfo::Type::HEXAHEDRON;
        lConn.resize(8);
        for(int i=0; i<8; i++){
            lConn[i] = *it;
            ++it;
        }

        patchVol->addConnectedCell(lConn, btype, id );
        patchVol->setPIDCell(id,PID);
        id++;
        break;

        case CGNS_ENUMV(NODE):
			                                                        ++it;
        break;

        case CGNS_ENUMV(BAR_2):
        case CGNS_ENUMV(BAR_3):

        for(int j=0; j<2; ++j)	++it;
        break;

        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):

        btype = bitpit::ElementInfo::Type::TRIANGLE;
        lConn.resize(3);
        for(int i=0; i<3; i++){
            lConn[i] = *it;
            ++it;
        }

        patchSurf->addConnectedCell(lConn, btype, id );
        patchSurf->setPIDCell(id,PID);
        id++;
        break;

        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):

        btype = bitpit::ElementInfo::Type::QUAD;
        lConn.resize(4);
        for(int i=0; i<4; i++){
            lConn[i] = *it;
            ++it;
        }

        patchSurf->addConnectedCell(lConn, btype, id );
        patchSurf->setPIDCell(id,PID);
        id++;
        break;

        default:
            std::cout<< "mimmo : ERROR : "<< m_name << " found unrecognized CGNS element while reading. Impossible to absorb further mixed elements. "<<std::endl;
            return;
            break;
        }
    } //end while

};	

/*! It recovers CGNS Info from linked Mimmo Objects.
 * It fills all the structures needed to export the volume and (eventually) surface
 * boundary meshes in CGNS format.
 *
 */
void
IOCGNS::recoverCGNSInfo(){

    MimmoObject * vol = getGeometry();
    MimmoObject * bnd = getSurfaceBoundary();

    /* Verifiy if a surface mesh (boundary mesh) is set. */
    bool flagBnd = (bnd != NULL);

    if( vol == NULL && bnd == NULL ) return;
    if( vol->isEmpty() && bnd->isEmpty() ) return;

    int nVolElements = 0;

    /* Fill volume map info. */
    CG_ElementType_t cgtype;
    long int ID;
    for (auto & cell : vol->getCells()){
        cgtype = m_storedInfo->mcg_typeconverter[cell.getType()];
        ID = cell.getId();
        m_storedInfo->mcg_number[cgtype]++;
        m_storedInfo->mcg_typetoid[cgtype].push_back(ID);
        long int * conn = cell.getConnect();
        for (int iV=0; iV<cell.getVertexCount(); iV++){
            m_storedInfo->mcg_typetoconn[cgtype].push_back(conn[iV]);
        }
    }

    /*Iterate map volume info and fill local CGNS index info
     * (volume elements will be appended before boundary elements)
     */
    std::map<CG_ElementType_t, std::vector<long int> >::iterator it, itend;
    std::vector<long int>::iterator itV, itVend;
    itend = m_storedInfo->mcg_typetoid.end();
    int cgnsidx = 0;
    for ( it = m_storedInfo->mcg_typetoid.begin(); it != itend; ++it ){
        itVend = it->second.end();
        for ( itV = it->second.begin(); itV != itVend; ++itV ){
            cgnsidx++;
            ID = (*itV);
            m_storedInfo->mcg_idtoindex[ID] = cgnsidx;
            m_storedInfo->mcg_indextoid.push_back(ID);
        }
    }
    nVolElements = cgnsidx;


    if (flagBnd){
        /* Fill surface mesh (boundary) map info. */
        bitpit::PatchKernel* bndpatch = bnd->getPatch();
        CG_ElementType_t cgtype;
        long int ID;
        for (auto & cell : bnd->getCells()){
            cgtype = m_storedInfo->mcg_typeconverter[cell.getType()];
            ID = cell.getId();
            m_storedInfo->mcg_bndnumber[cgtype]++;
            m_storedInfo->mcg_bndtypetoid[cgtype].push_back(ID);
            long int * conn = cell.getConnect();
            for (int iV=0; iV<cell.getVertexCount(); iV++){
                m_storedInfo->mcg_bndtypetoconn[cgtype].push_back(conn[iV]);
            }
        }

        /*Iterate map surface (boundary) info and fill local CGNS index info
         * (surface elements will be appended after volume elements)
         */
        std::map<CG_ElementType_t, std::vector<long int> >::iterator it, itend;
        std::vector<long int>::iterator itV, itVend;
        itend = m_storedInfo->mcg_bndtypetoid.end();
        int bndcgnsidx = 0;
        for ( it = m_storedInfo->mcg_bndtypetoid.begin(); it != itend; ++it ){
            itVend = it->second.end();
            for ( itV = it->second.begin(); itV != itVend; ++itV ){
                bndcgnsidx++;
                ID = (*itV);
                m_storedInfo->mcg_bndidtoindex[ID] = nVolElements + bndcgnsidx;
                m_storedInfo->mcg_bndindextoid.push_back(ID);
                m_storedBC->mcg_bndpidtoindex[bndpatch->getCell(ID).getPID()].push_back(bndcgnsidx + nVolElements);
            }
        }
    }

}




/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class reads the following parameters:
 *
 *  --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>ReadFlag</B>: activate reading mode boolean 1-reading mode, 0-writing mode
 * - <B>ReadDir</B>: reading directory path
 * - <B>ReadFilename</B>: name of file for reading
 * - <B>WriteDir</B>: writing directory path
 * - <B>WriteFilename</B>: name of file for writing
 *
 * \param[in]   slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void IOCGNS::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss>>value;
        }
        setPriority(value);
    };

    if(slotXML.hasOption("ReadFlag")){
        input = slotXML.get("ReadFlag");
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::trim(input));
            ss >> value;
        }
        setRead(value);
    };

    if(slotXML.hasOption("ReadDir")){
        input = slotXML.get("ReadDir");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setReadDir(input);
    };

    if(slotXML.hasOption("ReadFilename")){
        input = slotXML.get("ReadFilename");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "mimmoGeometry";
        setReadFilename(input);
    };

    if(slotXML.hasOption("WriteDir")){
        input = slotXML.get("WriteDir");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "./";
        setWriteDir(input);
    };

    if(slotXML.hasOption("WriteFilename")){
        input = slotXML.get("WriteFilename");
        input = bitpit::utils::trim(input);
        if(input.empty())   input = "mimmoGeometry";
        setWriteFilename(input);
    };

};

/*!
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class writes the following parameters(if different from default):
 *
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.IOCGNS"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>ReadFlag</B>: activate reading mode boolean 1-reading mode, 0-writing mode
 * - <B>ReadDir</B>: reading directory path
 * - <B>ReadFilename</B>: name of file for reading
 * - <B>WriteDir</B>: writing directory path
 * - <B>WriteFilename</B>: name of file for writing
 *
 * \param[in]   slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void IOCGNS::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));


    std::string output;

    output = std::to_string(m_read);
    slotXML.set("ReadFlag", output);

    slotXML.set("ReadDir", m_rdir);

    slotXML.set("ReadFilename", m_rfilename);

    slotXML.set("WriteDir", m_wdir);

    slotXML.set("WriteFilename", m_wfilename);

};


};
