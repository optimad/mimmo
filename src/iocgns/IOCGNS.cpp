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

#include "IOCGNS.hpp"
#include <cgnslib.h>

namespace mimmo{

/*!
 * Default constructor of IOCGNS.
 * \param[in] mode Define the working mode.
 */
IOCGNS::IOCGNS(IOCGNS::IOCGNS_Mode mode){
    setDefaults();
    m_mode = mode;
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
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.IOCGNS"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of IOCGNS.
 */
IOCGNS::~IOCGNS(){};

/*!
 * Copy constructor of IOCGNS.
 */
IOCGNS::IOCGNS(const IOCGNS & other):BaseManipulation(other){
    m_mode = other.m_mode;
    m_dir = other.m_dir;
    m_filename = other.m_filename;
    m_surfmesh_not = other.m_surfmesh_not;
    m_writeOnFile = other.m_writeOnFile;
    m_wtype = other.m_wtype;
    m_multizone = other.m_multizone;
    m_elementsSectionName = other.m_elementsSectionName;

    m_storedBC = std::move(std::unique_ptr<BCCGNS>(new BCCGNS(*(other.m_storedBC.get()))));
};

/*!
 * Assignement operator of IOCGNS.
 */
IOCGNS & IOCGNS::operator=(IOCGNS other){
    swap(other);
    return *this;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void IOCGNS::swap(IOCGNS & x) noexcept
{
    std::swap(m_mode, x.m_mode);
    std::swap(m_dir, x.m_dir);
    std::swap(m_filename, x.m_filename);
    std::swap(m_surfmesh_not, x.m_surfmesh_not);
    std::swap(m_volmesh, x.m_volmesh);
    std::swap(m_surfmesh, x.m_surfmesh);
    std::swap(m_storedBC, x.m_storedBC);

    std::swap(m_writeOnFile, x.m_writeOnFile);
    std::swap(m_wtype, x.m_wtype);
    std::swap(m_multizone, x.m_multizone);
    std::swap(m_elementsSectionName, x.m_elementsSectionName);

    BaseManipulation::swap(x);
}

/*!
 * Default values for IOCGNS.
 */
void
IOCGNS::setDefaults(){

    m_name      = "mimmo.IOCGNS";
    m_mode      = IOCGNS_Mode::READ;
    m_filename  = "mimmoVolCGNS";
    m_dir       = "./";
    m_surfmesh_not = NULL;
//    m_storedInfo = std::move(std::unique_ptr<InfoCGNS>(new InfoCGNS()));
    m_storedBC  = std::move(std::unique_ptr<BCCGNS>(new BCCGNS()));

    m_writeOnFile = false;
    m_wtype = IOCGNS_WriteType::ADF;
    m_multizone = false;

    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(TETRA_4))] = "Elem_tetra";
    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(PYRA_5))]  = "Elem_pyra";
    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(PENTA_6))] = "Elem_prism";
    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(HEXA_8))]  = "Elem_hexa";
    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(TRI_3))]   = "Elem_tri";
    m_elementsSectionName[static_cast<int>(CGNS_ENUMV(QUAD_4))]  = "Elem_quad";

}

/*!
 * It builds the input/output ports of the object
 */
void
IOCGNS::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoObject*, IOCGNS>(this, &IOCGNS::setGeometry, M_GEOM));
    built = (built && createPortIn<MimmoObject*, IOCGNS>(this, &IOCGNS::setSurfaceBoundary, M_GEOM2));
    built = (built && createPortIn<BCCGNS*, IOCGNS>(this, &IOCGNS::setBoundaryConditions, M_BCCGNS));

    built = (built && createPortOut<MimmoObject*, IOCGNS>(this, &IOCGNS::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoObject*, IOCGNS>(this, &IOCGNS::getSurfaceBoundary, M_GEOM2));
    built = (built && createPortOut<BCCGNS*, IOCGNS>(this, &IOCGNS::getBoundaryConditions, M_BCCGNS));
    m_arePortsBuilt = built;
};

/*!
 * Return all the surface bounding the current volume mesh.
 * During reading mode returns info contained in cgns, marking with PID all the possible boundary patches,
 * while in writing mode refers to the actual object pointed by the User.
 * \return pointer to surface boundary mesh
 */
MimmoObject*
IOCGNS::getSurfaceBoundary(){

    MimmoObject* res = nullptr;

    switch(m_mode){
        case IOCGNS_Mode::READ:
        case IOCGNS_Mode::RESTORE:
            res = m_surfmesh.get();
        break;
        case IOCGNS_Mode::WRITE:
        case IOCGNS_Mode::DUMP:
            res = m_surfmesh_not;
        break;
    }

    return res;
}

/*!
 * Return current pointer to geometry.If read mode return local read volumetric mesh, else
 * if in write mode return pointed externally geometry
 * \return pointer to volume mesh
 */
MimmoObject*
IOCGNS::getGeometry(){

    MimmoObject* res = nullptr;

    switch(m_mode){
        case IOCGNS_Mode::READ:
        case IOCGNS_Mode::RESTORE:
            res = m_volmesh.get();
        break;
        case IOCGNS_Mode::WRITE:
        case IOCGNS_Mode::DUMP:
            res = BaseManipulation::getGeometry();
        break;
    }

    return res;
}

/*!
 * Return the Info object (class BCCGNS) for boundary conditions of surface patches
 * as PIDs associated to Surface boundary mesh.
 * Meaningful only in reading mode.
 * \return pointer to BCCGNS object with boundary conditions information read from file.
 */
BCCGNS*
IOCGNS::getBoundaryConditions(){
    return m_storedBC.get();
};

/*!
 * Return the working mode of the class. See setMode method doc for detailed explanation.
 * \return working mode
 */
IOCGNS::IOCGNS_Mode
IOCGNS::whichMode(){
    return m_mode;
}

/*!
 *Overloading of whichMode, returning int flag;
 * \return working mode
 */
int
IOCGNS::whichModeInt(){
    return static_cast<int>(m_mode);
}

/*!
  Check if the class is writing native cgns with HDF5 format (1) or
  ADF format(2) or ADF2 format(3). NONE format is marked as 0.
  The method is meaningful only in class mode IOCGNS_Mode::WRITE.
  \return format type.
 */
IOCGNS::IOCGNS_WriteType    IOCGNS::whatWritingFormat(){
    return m_wtype;
}
/*!
  Check if the class is set to write cgns multizone or not.
  The method is meaningful only in class mode IOCGNS_Mode::WRITE.
  \return boolean true-writing multizone, false-writing single zone.
 */
bool    IOCGNS::isWritingMultiZone(){
    return m_multizone;
}


/*!It sets the  working directory path for IO operation.
   File to be read or file to be written will be located here.
   \param[in] dir path to directory
 */
void
IOCGNS::setDir(const std::string &dir){
    m_dir = dir;
}

/*!It sets the  filename without tag (.***) for IO operation.
   Given a target name "test", for each mode ew will have:
   - READ        : read from Dir the cgns file test.cgns (MPI with 0-rank only).
   - RESTORE     : read from Dir the dump file test.xxx.dump.
   - WRITE       : write to Dir the cgns file test.cgns, containing the inner mesh (MPI with 0-rank only).
   - DUMP        : write to Dir the dump file test.xxx.dump, containing the inner mesh.

   In case setWriteOnFileMeshInfo is set to true (or WriteInfo = 1 in xml) the mesh info file
   test_MeshInfo.dat will be written on Dir, for READ and WRITE mode only .

 * \param[in] filename Name of target file.
 */
void
IOCGNS::setFilename(const std::string & filename){
    m_filename = filename;
}


/*!It sets the  working mode of the class.
   Given a target filename "test" and a target dir "pathdir/"
   - READ        : read the mesh from abs path file pathdir/test.cgns
   - RESTORE     : restore the mesh from abs path file pathdir/test.xxx.dump.
   - WRITE       : write the mesh to abs path file pathdir/test.cgns.
   - DUMP        : write the mesh to dump abs path pathdir/test.xxx.dump.

 * \param[in] mode working mode of the class
 */
void
IOCGNS::setMode(IOCGNS::IOCGNS_Mode mode){
    m_mode = mode;
}


/*!Write info of the mesh on file, such as zone names, bc names, etc... .
 * The writing is active for READ and WRITE working mode. The save directory path is
 * specified through setDir method. File is writter by
 * \param[in] write boolean, if true write the info file.
 */
void
IOCGNS::setWriteOnFileMeshInfo(bool write){
    m_writeOnFile = write;
}

/*!
 * Set current geometry to an external volume mesh.
 * \param[in] geo Pointer to input volume mesh.
 */
void
IOCGNS::setGeometry(MimmoObject * geo){
    if(geo == nullptr )    return;
    if(geo->getType() != 2)    return;
    switch(m_mode){
        case IOCGNS_Mode::WRITE:
        case IOCGNS_Mode::DUMP:
            BaseManipulation::setGeometry(geo);
        break;
        default:
            //do nothing;
        break;
    }
}

/*!
 * Set boundary surface relative to the volume mesh.Option active only in writing mode.
 * Pre-existent boundary surfaces read from file, when class is set in read mode,  will be erased.
 * \param[in] geosurf Pointer to input volume mesh.
 */
void
IOCGNS::setSurfaceBoundary(MimmoObject* geosurf){
    if(geosurf == nullptr)    return;
    if(geosurf->getType() != 1)      return;
    switch(m_mode){
        case IOCGNS_Mode::WRITE:
        case IOCGNS_Mode::DUMP:
            m_surfmesh_not = geosurf;
        break;
        default:
            //do nothing;
        break;
    }
};

/*!
 * Set the Info object (class BCCGNS) for boundary conditions of surface patches
 * as PIDs associated to Surface boundary mesh.
 * Meaningful only in writing mode.
 * \param[in] bccgns Pointer to BCCGNS object with boundary conditions information read from file.
 */
void
IOCGNS::setBoundaryConditions(BCCGNS* bccgns){
    if(bccgns != NULL){
        std::unique_ptr<BCCGNS> temp(new BCCGNS(*bccgns));
        m_storedBC = std::move(temp);
    }
};


/*!
    Set write type cgns native format.
    The method is meaningful only in class mode IOCGNS_Mode::WRITE.
    See whatWritingFormat method doc for options available.
    \param[in] type of writing format
*/
void    IOCGNS::setWritingFormat(IOCGNS::IOCGNS_WriteType type){
    m_wtype = type;
}
/*!
    Set writing MultiZone cgns native format.
    The method is meaningful only in class mode IOCGNS_Mode::WRITE.
    WARNING: the class writes in single zone format for now.
    \param[in] multizone true write multi zone, false write single zone.
*/

void    IOCGNS::setWritingMultiZone(bool multizone){
    //TODO uncomment
    //m_multizone = multizone;
    BITPIT_UNUSED(multizone);
    m_multizone = false;
}


/*!Execution command.
 * It reads the geometry if the condition m_read is true.
 * It writes the geometry if the condition m_write is true.
 */
void
IOCGNS::execute(){

    std::string target = m_dir+"/"+m_filename+".cgns";

    switch(m_mode){
        case IOCGNS_Mode::READ:
            if(!read(target)){
                (*m_log) << "Error IOCGNS Reading mode: file corrupted or not found : "<< target << std::endl;
                (*m_log) << " " << std::endl;
                throw std::runtime_error ("file corrupted or not found : " + target);
            }
            writeInfoFile();
        break;
        case IOCGNS_Mode::DUMP:
            {
                int archiveVersion = 1;
                std::string header(m_name);
                std::string filedump = m_dir+"/"+m_filename;
        #if MIMMO_ENABLE_MPI
            	bitpit::OBinaryArchive binaryWriter(filedump, "dump", archiveVersion, header, m_rank);
        #else
            	bitpit::OBinaryArchive binaryWriter(filedump, "dump", archiveVersion, header);
        #endif
                if(!dump(binaryWriter.getStream())){
                    (*m_log) << "Error IOCGNS Dumping mode: impossible to write on : "<< filedump << std::endl;
                    (*m_log) << " " << std::endl;
                    throw std::runtime_error ("impossible to write dumpfile : " + filedump);
                }
                binaryWriter.close();
            }
        break;
        case IOCGNS_Mode::RESTORE:
        {
            std::string filedump = m_dir+"/"+m_filename;
            #if MIMMO_ENABLE_MPI
                	bitpit::IBinaryArchive binaryReader(filedump,"dump", m_rank);
            #else
                	bitpit::IBinaryArchive binaryReader(filedump,"dump");
            #endif

                m_surfmesh_not = nullptr;
                m_geometry = nullptr;
                m_volmesh = std::unique_ptr<MimmoObject>(new MimmoObject(2));
                m_surfmesh = std::unique_ptr<MimmoObject>(new MimmoObject(1));

                if(!restore(binaryReader.getStream())){
                    (*m_log) << "Error IOCGNS Restoring mode: impossible to restore from dump : "<< filedump << std::endl;
                    (*m_log) << " " << std::endl;
                    throw std::runtime_error ("impossible to read dumpfile : " + filedump);
                }
                binaryReader.close();
        }
        break;
        case IOCGNS_Mode::WRITE:
            if(!write(target)){
                (*m_log) << "Error IOCGNS Writing mode: impossible to write on : "<< target << std::endl;
                (*m_log) << " " << std::endl;
                throw std::runtime_error ("impossible to write : " + target);
            }
            writeInfoFile();

        break;
    }

}

/*!
    Write info file
 */
void IOCGNS::writeInfoFile(){

    std::string writeInfoFilename = m_dir+"/"+m_filename+"_MeshInfo.dat";
    std::string target = m_dir+"/"+m_filename+".cgns";

#if MIMMO_ENABLE_MPI
    if(m_rank == 0){
#endif
    if(m_writeOnFile){

        std::ofstream out;
        out.open(writeInfoFilename);
        if(out.is_open()){
            out<<"Info on :" << target<< " CGNS unstructured volume mesh"<<std::endl;
            out<<std::endl;
            out<<"Number of Vertices                    : "<<getGeometry()->getPatch()->getVertexCount()<<std::endl;
            out<<"Number of Cells                       : "<<getGeometry()->getPatch()->getCellCount()<<std::endl;
            out<<"Number of Patched Boundary Face Cells : "<<getSurfaceBoundary()->getPatch()->getCellCount()<<std::endl;
            out<<std::endl;
            out<<std::endl;

            if(m_storedBC){
                out<<"Zones defined on the mesh :"<<std::endl;
                out<<std::endl;
                for(auto & pl : m_storedBC->mcg_zonepidnames){
                    out<<"PID : "<<pl.first<<" Name : "<<pl.second<<std::endl;
                }
                out<<std::endl;
                out<<std::endl;

                out<<"Boundary Patches defined on the mesh :"<<std::endl;
                out<<std::endl;
                for(auto & pl : m_storedBC->mcg_bcpidnames){
                    out<<"PID : "<<pl.first<<" Name : "<<pl.second<<std::endl;
                }
            }

            out.close();
        }else{
            (*m_log) << "Warning IOCGNS : cannot open "<<writeInfoFilename<<" to flush mesh Info" << std::endl;
        }
    }
#if MIMMO_ENABLE_MPI
    } //endif m_rank==0
#endif
}

/*!It reads the mesh geometry from an input file.
   \param[in] file abs path to read cgns.
   \return False if problems occur during reading stage.
 */
bool
IOCGNS::read(const std::string & file){

    m_volmesh = std::unique_ptr<MimmoObject>(new MimmoObject(2));
    m_surfmesh = std::unique_ptr<MimmoObject>(new MimmoObject(1));
    std::unique_ptr<MimmoObject> patchVol(new MimmoObject(2));

#if MIMMO_ENABLE_MPI
    if(m_rank == 0){
#endif
    std::string error_string = "read CGNS grid: " + file;

    std::vector<long> nVertices;
    std::vector<long> nCells;
    //     long nCells, nBoundVertices;
    std::vector<int> nCoords; // how many coordinates 3 - x y z
    std::vector< std::array< std::vector<double>,3 > > zonecoords; //node coordinates for each zone
    std::vector<int> nSections; // haw many grouped sections of homogeneous cell elements
    std::vector< std::vector<CGNS_ENUMT(ElementType_t)> > orderedConns; //orderedConnectivities - retain element type at index i, for each zone.
    std::vector<std::unordered_map<int, ivector1D> > conns; //connectivity to be absorbed for each zone
    std::vector<int> nBcs; // number of bc for each zone.
    std::vector<std::unordered_map<int, ivector1D >> bcLists; //list of elements composing bcs, for each zone.
    std::vector<std::unordered_map<int, std::string >> bcNames; // list of names of bcs, for each zone
    std::vector<std::unordered_map<int, CGNS_ENUMT(BCType_t) >> bcTypes; // list of types of bcs, for each zone
    std::vector<std::unordered_map<int, bool > > bcOnElements;

    //FIRST STEP ABSORB INFO from FILE//

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
        (*m_log) << m_name << " more than one base found in grid file -> only the first zone will be read" << std::endl;
    }

    //Read name of basis and physical dimension
    char basename[33];
    int physdim, celldim;
    if(cg_base_read(indexfile, 1, basename, &celldim, &physdim) != CG_OK){
        return false;
    }
    if(celldim != 3 || physdim !=3){
        //Only volume mesh supported
        return false;
    }

    //Read number of zone in basis, check if only one.
    int nzones;
    if(cg_nzones(indexfile, 1, &nzones)!=CG_OK){
        return false;
    }

    //Read type mesh for each zones

    CGNS_ENUMT(ZoneType_t) zoneType[nzones];
    int index_dim[nzones];
    std::vector<std::vector<cgsize_t>> sizeG(nzones, std::vector<cgsize_t>(3));
    std::vector<std::string> zonenames(nzones);
    nVertices.resize(nzones);
    nCells.resize(nzones);
    nCoords.resize(nzones);
    zonecoords.resize(nzones);
    nSections.resize(nzones);
    conns.resize(nzones);
    orderedConns.resize(nzones);
    nBcs.resize(nzones);
    bcLists.resize(nzones);
    bcNames.resize(nzones);
    bcTypes.resize(nzones);
    bcOnElements.resize(nzones);

    for(int indexzone=0; indexzone<nzones; ++indexzone){
        if(cg_zone_type(indexfile,1,indexzone+1, &zoneType[indexzone])!= CG_OK){
            return false;
        }
        if(cg_index_dim(indexfile,1,indexzone+1, &index_dim[indexzone])!= CG_OK){
            return false;
        }
        if(zoneType[indexzone] != CGNS_ENUMT(ZoneType_t)::CGNS_ENUMV(Unstructured) || index_dim[indexzone] != 1){
            //Only unstructured mesh supported for now (index_dim == 1 for unstructured)
            return false;
        }
        //Read size of zone (n nodes, n cells, n boundary nodes)
        char fzz[33];
        if(cg_zone_read(indexfile,1,indexzone+1, fzz, sizeG[indexzone].data()) != CG_OK ){
            return false;
        }
        zonenames[indexzone] = std::string(fzz);

        nVertices[indexzone]= sizeG[indexzone][0];
        nCells[indexzone]= sizeG[indexzone][1];
        //Read Vertices.
        if(cg_ncoords(indexfile,1,indexzone+1, &nCoords[indexzone])!= CG_OK){
            return false;
        }


        for(int i = 0; i < nCoords[indexzone]; ++i){
            cgsize_t startIndex = 1;
            cgsize_t finishIndex = sizeG[indexzone][0];
            CGNS_ENUMT(DataType_t) datatype;
            char name[33];
            zonecoords[indexzone][i].resize(nVertices[indexzone]);
            if(cg_coord_info(indexfile,1,indexzone+1,i+1, &datatype, name)!=CG_OK){
                return false;
            }
            if(cg_coord_read(indexfile,1,indexzone+1,name, datatype, &startIndex, &finishIndex, zonecoords[indexzone][i].data() )!=CG_OK){
                return false;
            }
        }

        //Read connectivities.
        //They are read starting from 1, fortran style. When matching up w/ coords
        //vector positions remember to diminish conn value of 1.
        //Read number of sections of connectivity structure

        if(cg_nsections(indexfile,1,indexzone+1, &nSections[indexzone])!= CG_OK){
            return false;
        }

        for(int sec = 0; sec < nSections[indexzone]; ++sec){

            //Read section sec
            char elementname[33];
            CGNS_ENUMT(ElementType_t) type;
            cgsize_t eBeg, eEnd;
            int enBdry;
            int parent_flag;
            std::vector<cgsize_t>    connlocal;
            cgsize_t size;

            //Read elements name, type and range
            if(cg_section_read(indexfile,1,indexzone+1, sec+1, elementname, &type, &eBeg, &eEnd,&enBdry, & parent_flag)!= CG_OK){
                return false;
            }

            //Read size of connectivity data
            if(cg_ElementDataSize(indexfile,1,indexzone+1, sec+1,&size)!= CG_OK){
                return false;
            }

            connlocal.resize((size_t) size);
            if(cg_elements_read(indexfile,1,indexzone+1,sec+1, connlocal.data(),nullptr) !=CG_OK){
                return false;
            }

            conns[indexzone][sec].resize(connlocal.size());
            int count=0;
            for(const auto &val: connlocal){
                conns[indexzone][sec][count] = (int)val-1; //from fortran to c indexing.
                ++count;
            }
            orderedConns[indexzone].push_back(type);
        }


        //explore superficial boundary conditions definition, for boundary surface extraction.
        if(cg_nbocos(indexfile,1,indexzone+1, &nBcs[indexzone])!= CG_OK){
            return false;
        }

        //Everything marked as bc for me it's a wall.
        for(int indexbc=0; indexbc<nBcs[indexzone]; ++indexbc){

            char name[33];
            CGNS_ENUMT(BCType_t) bocotype;
            CGNS_ENUMT(PointSetType_t) ptset_type;
            std::vector<cgsize_t> nBCElements(2);
            int normalIndex;
            cgsize_t normalListSize;
            CGNS_ENUMT(DataType_t) normalDataType;
            int ndataset;
            //CGNS_ENUMT(GridLocation_t) location;

            //reading information of the boundary. What i need here is
            // the name (for sewing betwenn zones after) and the ptset_type-nBCElements.
            if(cg_boco_info(indexfile,1,indexzone+1,indexbc+1, name, &bocotype, &ptset_type, nBCElements.data(),
                    &normalIndex, &normalListSize, &normalDataType, &ndataset) != CG_OK){
                return false;
            }
            GridLocation_t bclocation;
            if(cg_boco_gridlocation_read(indexfile,1,indexzone+1,indexbc+1, &bclocation) != CG_OK){
                return false;
            }

            bcNames[indexzone][indexbc] = std::string(name);
            bcTypes[indexzone][indexbc] = bocotype;

            bool validSetType = true;
            std::vector<cgsize_t> localbclist;

            switch(ptset_type){
                case CGNS_ENUMV(PointList):
                case CGNS_ENUMV(ElementList):
                {
                    localbclist.resize((size_t) nBCElements[0]);
                    if(cg_boco_read(indexfile,1,indexzone+1,indexbc+1, localbclist.data(), nullptr )!= CG_OK){
                        return false;
                    }
                }
                break;
                case CGNS_ENUMV(PointRange):
                case CGNS_ENUMV(ElementRange):
                {
                    cgsize_t dim = nBCElements[1] - nBCElements[0] + 1;
                    localbclist.reserve((size_t) dim);
                    for (cgsize_t idx = nBCElements[0]; idx <= nBCElements[1]; idx++){
                        localbclist.push_back(idx);
                    }
                }
                break;
                default:
                    validSetType = false;
                break;
            }
            if(!validSetType){
                (*m_log)<<"IOCGNS reader cannot support BC PointSetType_t different from PointList, PointRange, ElementList and ElementRange.Aborting"<<std::endl;
                return false;
            }

            bool flag = ( ptset_type == CGNS_ENUMV(ElementList) );
            flag = flag ||  (ptset_type == CGNS_ENUMV(ElementRange) );
            flag = flag ||  ( (ptset_type == CGNS_ENUMV(PointList) || ptset_type == CGNS_ENUMV(PointRange) )
                               && bclocation == CGNS_ENUMV(FaceCenter) ) ;
            bcOnElements[indexzone].insert(std::make_pair(indexbc,flag));
            bcLists[indexzone][indexbc].resize(localbclist.size());
            int count=0;
            for(const auto &val: localbclist){
                bcLists[indexzone][indexbc][count] = (int)val-1;//from fortran to c style
                ++count;
            }
        }
    }

    //Finish reading CGNS file
    cg_close(indexfile);

    //PUT THIS INFO IN MimmoObject
    long PIDZoneOffset = 0;
    long PIDBCOffset = 1;
    long idVertexOffset = 0;
    long idCellOffset = 0;

    long totV(0), totC(0);
    for(int i=0; i<nzones; ++i){
        totV += nVertices[i];
        totC += nCells[i];
    }
    m_volmesh->getPatch()->reserveVertices(totV);
    m_volmesh->getPatch()->reserveCells(totC);
    std::unordered_map<long, std::map<int, long> > mapCellFacePid; //used for point list. idcell - face - PID
    std::unordered_map<long, std::unordered_map<long,std::vector<long> > > mapBndFaceDCPid; //used for element list. PID - id -2D element connectivity
    std::unordered_map<long, std::string > bndNames;

    for(int target_zone = 0; target_zone< nzones; ++target_zone){

        int nbc = nBcs[target_zone];

        for(int j=0; j<nbc; ++j){
             m_storedBC->mcg_pidtobc[PIDBCOffset +j] = bcTypes[target_zone][j];
             m_storedBC->mcg_bcpidnames[PIDBCOffset +j] = bcNames[target_zone][j];
             m_storedBC->mcg_zonetobndpid[target_zone].push_back(PIDBCOffset+j);
        }

        //Reverse info in your grids.
        patchVol->resetPatch();
        patchVol->getPatch()->reserveVertices(nVertices[target_zone]);
        patchVol->getPatch()->reserveCells(nCells[target_zone]);

        std::unordered_map<long, std::vector<long>> flaggedBCConns;

        long id;
        darray3E temp;

        //    Stock vertices in volume grid
        for(int i=0; i<nVertices[target_zone]; ++i){
            for(int j=0; j<3; ++j)    temp[j] = zonecoords[target_zone][j][i];
            id = idVertexOffset + i;
            patchVol->addVertex(temp,id);
        }

        id = 0;
        int idsection = -1;
        //Unpack 3D elements connectivities and store it in volume grid. Label same species elements w PID.
        for(const auto & val: orderedConns[target_zone]){
            idsection++;

            livector1D lConn;
            bitpit::ElementType btype;
            long PIDZone = target_zone + PIDZoneOffset;
            int size = conns[target_zone][idsection].size();

            switch(val){

            case CGNS_ENUMV(TETRA_4):
            case CGNS_ENUMV(TETRA_10):
            case CGNS_ENUMV(TETRA_16):
            case CGNS_ENUMV(TETRA_20):
            case CGNS_ENUMV(TETRA_22):
            case CGNS_ENUMV(TETRA_34):
            case CGNS_ENUMV(TETRA_35):

            btype = bitpit::ElementType::TETRA;
            lConn.resize(4);
            for(int i=0; i<size; i+=4){ //first 4 element are needed only
                for(int j=0; j<4; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                patchVol->addConnectedCell(lConn, btype, id+idCellOffset );
                patchVol->setPIDCell(id+idCellOffset,PIDZone);
                id++;
            }
            break;

            case CGNS_ENUMV(PYRA_5):
            case CGNS_ENUMV(PYRA_13):
            case CGNS_ENUMV(PYRA_14):
            case CGNS_ENUMV(PYRA_21):
            case CGNS_ENUMV(PYRA_29):
            case CGNS_ENUMV(PYRA_30):
            case CGNS_ENUMV(PYRA_50):
            case CGNS_ENUMV(PYRA_55):

            btype = bitpit::ElementType::PYRAMID;
            lConn.resize(5);
            for(int i=0; i<size; i+=5){
                for(int j=0; j<5; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                patchVol->addConnectedCell(lConn, btype, id+idCellOffset );
                patchVol->setPIDCell(id+idCellOffset,PIDZone);
                id++;
            }
            break;

            case CGNS_ENUMV(PENTA_6):
            case CGNS_ENUMV(PENTA_15):
            case CGNS_ENUMV(PENTA_18):
            case CGNS_ENUMV(PENTA_24):
            case CGNS_ENUMV(PENTA_38):
            case CGNS_ENUMV(PENTA_40):
            case CGNS_ENUMV(PENTA_33):
            case CGNS_ENUMV(PENTA_66):
            case CGNS_ENUMV(PENTA_75):

            btype = bitpit::ElementType::WEDGE;
            lConn.resize(6);
            for(int i=0; i<size; i+=6){
                for(int j=0; j<6; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                //remap in bitpit conn. TODO complete ref element mapper for connectivity.
                std::swap(lConn[1], lConn[2]);
                std::swap(lConn[4], lConn[5]);

                patchVol->addConnectedCell(lConn, btype, id+idCellOffset );
                patchVol->setPIDCell(id+idCellOffset,PIDZone);
                id++;
            }

            break;

            case CGNS_ENUMV(HEXA_8):
            case CGNS_ENUMV(HEXA_20):
            case CGNS_ENUMV(HEXA_27):
            case CGNS_ENUMV(HEXA_32):
            case CGNS_ENUMV(HEXA_56):
            case CGNS_ENUMV(HEXA_64):
            case CGNS_ENUMV(HEXA_44):
            case CGNS_ENUMV(HEXA_98):
            case CGNS_ENUMV(HEXA_125):

            btype = bitpit::ElementType::HEXAHEDRON;
            lConn.resize(8);
            for(int i=0; i<size; i+=8){
                for(int j=0; j<8; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                patchVol->addConnectedCell(lConn, btype, id+idCellOffset );
                patchVol->setPIDCell(id+idCellOffset,PIDZone);
                id++;
            }
            break;

            case CGNS_ENUMV(TRI_3):
            case CGNS_ENUMV(TRI_6):
            case CGNS_ENUMV(TRI_9):
            case CGNS_ENUMV(TRI_10):
            case CGNS_ENUMV(TRI_12):
            case CGNS_ENUMV(TRI_15):

            btype = bitpit::ElementType::TRIANGLE;
            lConn.resize(3);

            for(int i=0; i<size; i+=3){
                for(int j=0; j<3; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                flaggedBCConns[id] = lConn;
                id++;
            }
            break;

            case CGNS_ENUMV(QUAD_4):
            case CGNS_ENUMV(QUAD_8):
            case CGNS_ENUMV(QUAD_9):
            case CGNS_ENUMV(QUAD_12):
            case CGNS_ENUMV(QUAD_16):
            case CGNS_ENUMV(QUAD_25):

            btype = bitpit::ElementType::QUAD;
            lConn.resize(4);
            for(int i=0; i<size; i+=4){
                for(int j=0; j<4; ++j){
                    lConn[j] = idVertexOffset + conns[target_zone][idsection][i+j];
                }
                flaggedBCConns[id] = lConn;
                id++;
            }
            break;

            case CGNS_ENUMV(MIXED):
                unpackMixedConns(conns[target_zone][idsection], patchVol.get(), flaggedBCConns, idVertexOffset, idCellOffset, PIDZone, id);
            break;
            default:
                (*m_log)<<"Warning IOCGNS : found unsupported Element Type in CGNS connectivity"<<std::endl;
            break;
            }
        }

        // prepare data to be passed to surface mesh.
        // Beware: if not defined by points, surface elements are specified in
        // the previous step of connectivity and stored in flaggedBCConns.
        // Now you have to recognize volume faces that holds this boundary points/ elements.

        //you have bc on points, you need to pass from patchVol first.
        bitpit::PiercedVector<bitpit::Cell> & cells = patchVol->getCells();
        //bitpit::PiercedVector<bitpit::Vertex> & verts = patchVol->getVertices();
        std::unordered_map<long, std::set<int> > borderFaceCells;

        //loop on local bc;
        for(int j=0; j<nbc; ++j){

             m_storedBC->mcg_pidtolisttype[int(PIDBCOffset +j)] = int(bcOnElements[target_zone][j]);

             if(!bcOnElements[target_zone][j]){

                //compute borderFaceCells once.
                if(borderFaceCells.empty()) {
                    borderFaceCells = patchVol->extractBoundaryFaceCellID();
                }
                // reverse each singular list of bcs in a more suitable container
                std::set<long> pool;
                for(long idV: bcLists[target_zone][j]){
                    pool.insert(idV+idVertexOffset);
                }

                // run over boundary faces, get their connectivity. If each conn vertex is in one of the pools
                // mark it as a cell element of boundary
                for(const auto & cellPair : borderFaceCells){
                    bitpit::Cell & cell = cells.at(cellPair.first);

                    for(int iface : cellPair.second){
                        bitpit::ConstProxyVector<long> conn = cell.getFaceConnect(iface);
                        if( belongToPool(conn, pool)){
                            mapCellFacePid[cellPair.first].insert(std::make_pair(iface, PIDBCOffset +long(j)));
                        }
                    }
                }

            }else{

                for(long idC: bcLists[target_zone][j]){
                        mapBndFaceDCPid[PIDBCOffset +long(j)].insert(std::make_pair(idC, flaggedBCConns.at(idC)) );
                }
            }
        } //end of looping bcs.
        flaggedBCConns.clear();
        //clean up adjacencies. These portion is local, you need to append this structure to the real mesh manager after.
        if(patchVol->areAdjacenciesBuilt()) patchVol->resetAdjacencies();

        //reversing patchVol inside volmesh .
        for(auto it=patchVol->getPatch()->vertexBegin(); it != patchVol->getPatch()->vertexEnd(); ++it){
            m_volmesh->addVertex(*it, it->getId());
        }

        for(auto it=patchVol->getPatch()->cellBegin(); it != patchVol->getPatch()->cellEnd(); ++it){
            m_volmesh->addCell(*it, it->getId());
        }

        //store the name of volumetric pid (zone names)
        m_storedBC->mcg_zonepidnames[PIDZoneOffset+target_zone] = zonenames[target_zone];

        //update the offsets;
        PIDZoneOffset = 0;
        PIDBCOffset += nbc ;
        idVertexOffset = m_volmesh->getPatch()->getVertexCount();
        idCellOffset = m_volmesh->getPatch()->getCellCount();

        // release intermediate structure of cgns for target_zone.
        zonecoords[target_zone][0].clear();
        zonecoords[target_zone][1].clear();
        zonecoords[target_zone][2].clear();
        orderedConns[target_zone].clear();
        conns[target_zone].clear();
        bcLists[target_zone].clear();
        bcNames[target_zone].clear();
        bcTypes[target_zone].clear();

    }

    //delete coincident vertices in the mother volume.
    m_volmesh->getPatch()->deleteCoincidentVertices();

    //Squeeze the mother volume
    m_volmesh->getPatch()->squeeze();

    //now create the surface mesh, using the mapCellFacePid information.
    long totSV, totSC(0);
    for(auto & pp : mapCellFacePid){
        totSC += pp.second.size();
    }
    for(auto & pp : mapBndFaceDCPid){
        totSC += pp.second.size();
    }

    totSV = 4*totSC;

    m_surfmesh->getPatch()->reserveVertices(totSV);
    m_surfmesh->getPatch()->reserveCells(totSC);


    bitpit::PiercedVector<bitpit::Cell> & volCells = m_volmesh->getCells();
    bitpit::PiercedVector<bitpit::Vertex> & volVerts = m_volmesh->getVertices();
    bitpit::PiercedVector<bitpit::Vertex> & surfVerts = m_surfmesh->getVertices();

    long PIDSurf;
    long idCellCounter(0);
    int iface;
    bitpit::ElementType et;
    bitpit::ConstProxyVector<long> conn;

    // fill data from point list to surface
    for(auto & mapp : mapCellFacePid){
        bitpit::Cell & cell = volCells.at(mapp.first);

        for(auto & info: mapp.second){
            iface = info.first;
            PIDSurf = info.second;
            conn = cell.getFaceConnect(iface);
            et = cell.getFaceType(iface);

            //push vertices in m_surfmesh;
            for(long idV : conn){
                if(!surfVerts.exists(idV)){
                    m_surfmesh->addVertex(volVerts.at(idV), idV);
                }
            }
            m_surfmesh->addConnectedCell(std::vector<long>(conn.begin(), conn.end()), et, PIDSurf, idCellCounter);
            ++idCellCounter;
        }
    }

    // fill data from element list to surface
    for(auto & rootmap : mapBndFaceDCPid){
        PIDSurf = rootmap.first;
        for(auto & info: rootmap.second){
            bitpit::ElementType et = bitpit::ElementType::QUAD;
            if(info.second.size()<4){
                et = bitpit::ElementType::TRIANGLE;
            }
            //push vertices in m_surfmesh;
            for(long idV : info.second){
                if(!surfVerts.exists(idV)){
                    m_surfmesh->addVertex(volVerts.at(idV), idV);
                }
            }
            m_surfmesh->addConnectedCell(info.second, et, PIDSurf, idCellCounter);
            ++idCellCounter;
        }
    }


#if MIMMO_ENABLE_MPI
    } //endif mpi
    MPI_Barrier(m_communicator);
    //make sure all procs know the m_storedBC info absorbed while reading.
    communicateAllProcsStoredBC();
#endif

    //Squeeze the surface mesh
    m_surfmesh->getPatch()->squeeze();

    return true;
}

/*!It writes the mesh geometry on output .cgns file.
   If boundary surface, PID subdivided, is available, write all patches as CGNS Wall boundary condition .
  \param[in] file abs path file to write mesh.
  \return False if valid volume geometry or surface geometry is not found.
 */
bool
IOCGNS::write(const std::string & file){

switch(m_wtype){
    case IOCGNS_WriteType::HDF5:
        cg_set_file_type(CG_FILE_HDF5);
        break;
    case IOCGNS_WriteType::ADF:
        cg_set_file_type(CG_FILE_ADF);
        break;
    case IOCGNS_WriteType::ADF2:
        cg_set_file_type(CG_FILE_ADF2);
        break;
    default: //NONE
        (*m_log)<<"Error in IOCGNS::write: a NONE file type is set as writing type format"<<std::endl;
        break;
}

if(m_multizone){
    //TODO : write better the multizone part and add the zone-zone connectivity
    //information maybe faces? 1to1Connectivity.
}else{

#if MIMMO_ENABLE_MPI
    if(m_rank == 0){
#endif

    std::string error_string = "write CGNS grid: " + file;

    /* Fill this structures (surface boundary connectivity has
     * to be filled with same vertex indices of volume mesh) */
    MimmoObject * vol = getGeometry();
    MimmoObject * bnd = getSurfaceBoundary();


    if( vol == nullptr || bnd == nullptr ) return false;

    //just in case resynchronize the internal pids - to be sure
    vol->resyncPID();
    bnd->resyncPID();

    //Open index nad Write Unique Base Info
    int indexfile;
    if(cg_open(file.c_str(), CG_MODE_WRITE, &indexfile) != CG_OK){
        (*m_log) << "error: cgns error during write: opening file " << file << std::endl;
        throw std::runtime_error ("cgns error during write " + file);
    }

    int baseindex = 1;
    char basename[33] = "Base0001";
    int physdim=3, celldim=3;
    if(cg_base_write(indexfile,basename, celldim, physdim, &baseindex) != CG_OK){
        (*m_log) << "error: cgns error during write : base  " << file << std::endl;
        throw std::runtime_error ("cgns error during write " + file);
    }

    //std::unordered_map<long, std::string> & zonePidList = vol->getPIDTypeListWNames();

    std::string zonename = "Zone0001";

    livector1D cellIds, vertIds;
    std::unordered_map<long, int> globToLoc;
    std::map<long, std::vector<cgsize_t>> bndPools;

    int zoneindex=1;
    //take out volume cell ids in the target zone.
    cellIds = vol->getCells().getIds();
    vertIds = vol->getVertices().getIds();


    // fill the inverse point map. Start from 1 because CGNS is a fortran buddy
    int countvert = 1;
    for(long idV : vertIds){
        globToLoc[idV] = countvert;
        ++countvert;
    }

    int bndsize = m_storedBC->mcg_pidtobc.size();
    livector1D temp;
    if(bndsize>0){
        for(auto & val : m_storedBC->mcg_pidtobc){
            if(m_storedBC->mcg_pidtolisttype[val.first] > 0){
                temp = bnd->extractPIDCells(val.first);
                bndPools[val.first].resize(temp.size());
                int llcc = 0;
                for(long id : temp){
                    bndPools[val.first][llcc] = id;
                    ++llcc;
                }
            }else{
                temp = bnd->getVertexFromCellList(bnd->extractPIDCells(val.first));
                bndPools[val.first].resize(temp.size());
                int llcc = 0;
                for(long id : temp){
                    bndPools[val.first][llcc] = globToLoc[id];
                    ++llcc;
                }
            }
        }
        // cgns boundaries for each zone are always unique.
    }

    CGNS_ENUMT(ZoneType_t) zoneType =CGNS_ENUMT(ZoneType_t)::CGNS_ENUMV(Unstructured) ;
    std::vector<cgsize_t> sizeG(3);
    sizeG[0] = vertIds.size();
    sizeG[1] = cellIds.size();
    sizeG[2] = 0; //unsorted elements.

    if(cg_zone_write(indexfile,baseindex, zonename.data(), sizeG.data(), zoneType, &zoneindex) != CG_OK ){
        (*m_log) << "error: cgns error during write: zone " << file << std::endl;
        throw std::runtime_error ("cgns error during write " + file);
    }

    //writing vertices.
    svector1D names(3, "CoordinateX");
    names[1] = "CoordinateY";
    names[2] = "CoordinateZ";

    CGNS_ENUMT(DataType_t) datatype= CGNS_ENUMV(RealDouble);
    {
        //put vertices coordinates in a more suitable structure
        std::array<std::vector<double>,3 > coords;
        coords[0].reserve(vertIds.size());
        coords[1].reserve(vertIds.size());
        coords[2].reserve(vertIds.size());

        darray3E temp;
        for(long idV : vertIds){
            temp = vol->getVertexCoords(idV);
            coords[0].push_back(temp[0]);
            coords[1].push_back(temp[1]);
            coords[2].push_back(temp[2]);
        }

        for(int i=1; i<=3; ++i){
            int index;
            if(cg_coord_write(indexfile,baseindex,zoneindex, datatype, names[i-1].data(), coords[i-1].data(), &index )!=CG_OK){
                (*m_log) << "error: cgns error during write: node coordinates " << file << std::endl;
                throw std::runtime_error ("cgns error during write " + file);
            }
        }
    }//end scope vertices

    std::map<int, std::size_t> mmcgns_ncells;
    std::map<int, std::vector<std::size_t> > mmcgns = getZoneConn(cellIds, globToLoc, mmcgns_ncells);
    livector1D bndCellIds;
    for(auto & mapp : bndPools){
        if(m_storedBC->mcg_pidtolisttype[int(mapp.first)] > 0){
            bndCellIds.insert(bndCellIds.end(), mapp.second.begin(), mapp.second.end());
        }
    }
    std::map<int, std::size_t> bndcgns_ncells;
    std::unordered_map<long, int> surfCellGlobToLoc;
    std::map<int, std::vector<std::size_t> > bndcgns = getBCElementsConn(bndCellIds, globToLoc, bndcgns_ncells, surfCellGlobToLoc);

    /* Write volume elements */
    int sec, countsec(1);
    cgsize_t eBeg = 1, eEnd;
    std::string sectionname;
    std::vector<cgsize_t> cgtemp;
    std::size_t locsize;

    //write voume elements.
    for(auto &connMap : mmcgns){

        if(connMap.second.empty()) continue;

        locsize = connMap.second.size();
        eEnd = eBeg + mmcgns_ncells[connMap.first] - 1;
        cgtemp.clear();
        cgtemp.reserve(locsize);
        for(std::size_t val : connMap.second){
            cgtemp.push_back(cgsize_t(val));
        }

        sectionname = "UndefElements";
        if(m_elementsSectionName.count(connMap.first) > 0){
            sectionname = m_elementsSectionName[connMap.first];
        }

        if(cg_section_write(indexfile,baseindex,zoneindex,sectionname.data(), static_cast<CGNS_ENUMT(ElementType_t)>(connMap.first),
                            eBeg,eEnd,0, cgtemp.data(), &sec) !=CG_OK )
        {
            cg_error_print();
            (*m_log) << "error: cgns error during write: section " << file << std::endl;
            throw std::runtime_error ("cgns error during write " + file);
        }
        eBeg = eEnd+1;
        ++countsec;
    }

    mmcgns.clear();
    cgsize_t startSurfElementsOffset(eBeg);
    //write surface elements.
    for(auto &connMap : bndcgns){

        if(connMap.second.empty()) continue;

        locsize = connMap.second.size();
        eEnd = eBeg + bndcgns_ncells[connMap.first] - 1;
        cgtemp.clear();
        cgtemp.reserve(locsize);
        for(std::size_t val : connMap.second){
            cgtemp.push_back(cgsize_t(val));
        }

        sectionname = "UndefElements";
        if(m_elementsSectionName.count(connMap.first) > 0){
            sectionname = m_elementsSectionName[connMap.first];
        }

        if(cg_section_write(indexfile,baseindex,zoneindex,sectionname.data(), static_cast<CGNS_ENUMT(ElementType_t)>(connMap.first),
                            eBeg,eEnd,0, cgtemp.data(), &sec) !=CG_OK )
        {
            cg_error_print();
            (*m_log) << "error: cgns error during write: section " << file << std::endl;
            throw std::runtime_error ("cgns error during write " + file);
        }
        eBeg = eEnd+1;
        ++countsec;
    }

    bndcgns.clear();

    /* Write boundary conditions
     */

    int bccount(0), bcid;

    for(auto & pool: bndPools){

        int pid = pool.first;
        CGNS_ENUMT(BCType_t) bocotype = static_cast<CGNS_ENUMT(BCType_t)>(m_storedBC->mcg_pidtobc[pid]);
        std::string bcname = "Undefined_BC_"+ std::to_string(pid);
        if(m_storedBC->mcg_bcpidnames.count(pid) > 0) bcname = m_storedBC->mcg_bcpidnames[pid];
        CGNS_ENUMT(PointSetType_t) ptset_type;
        if(m_storedBC->mcg_pidtolisttype[pid] == 0){
            ptset_type = CGNS_ENUMV(PointList);
        }else{
            ptset_type = CGNS_ENUMV(ElementList);
            //remap the ids of elements into the pool.
            for(cgsize_t & idSC : pool.second){
                idSC = surfCellGlobToLoc[idSC] + int(startSurfElementsOffset);
            }
        }
        cgsize_t nelems = pool.second.size();

        if(cg_boco_write(indexfile, baseindex, zoneindex, bcname.data(),
                bocotype, ptset_type, nelems, pool.second.data(), &bcid )!= CG_OK)
        {
            (*m_log) << "error: cgns error during write: bc " << file << std::endl;
            throw std::runtime_error ("cgns error during write " + file);
        }
        ++bccount;
    }

    /* Finish writing CGNS file */
    cg_close(indexfile);

#if MIMMO_ENABLE_MPI
    } //ENDIF M-RANK mpi
    MPI_Barrier(m_communicator);
#endif

}//endif multizone;

    return true;
}



/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCGNS::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("IOCGNS_Mode")){
        input = slotXML.get("IOCGNS_Mode");
        int value = 0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss >> value;
        }
        value = std::min(3, std::max(0, value));
        setMode(static_cast<IOCGNS_Mode>(value));
    };

    if(slotXML.hasOption("Dir")){
        input = slotXML.get("Dir");
        input = bitpit::utils::string::trim(input);
        if(input.empty())   input = "./";
        setDir(input);
    };

    if(slotXML.hasOption("Filename")){
        input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        if(input.empty())   input = "mimmoGeometry";
        setFilename(input);
    };


    if(slotXML.hasOption("WriteInfo")){
        input = slotXML.get("WriteInfo");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >>value;
        };
        setWriteOnFileMeshInfo(value);
    };

    if(slotXML.hasOption("WriteFormat")){
        input = slotXML.get("WriteFormat");
        input = bitpit::utils::string::trim(input);
        int value = 2;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >>value;
        };

        value = std::min(3, std::max(value, 0));
        setWritingFormat(static_cast<IOCGNS_WriteType>(value));
    };

    if(slotXML.hasOption("WriteMultiZone")){
        input = slotXML.get("WriteMultiZone");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >>value;
        };
        setWritingMultiZone(value);
    };
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
IOCGNS::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("IOCGNS_Mode", std::to_string(whichModeInt()));
    slotXML.set("Dir", m_dir);
    slotXML.set("Filename", m_filename);
    slotXML.set("WriteInfo", std::to_string(int(m_writeOnFile)));
    slotXML.set("WriteFormat", std::to_string(static_cast<int>(whatWritingFormat())));
    slotXML.set("WriteMultiZone", std::to_string(int(isWritingMultiZone())));

};

/*!
    Check if all the vertices of an element connectivity belongs to a pool of vertices
    \param[in] elementconn elementConnectivity
    \param[in] pool        ordered list of mesh vertices
    \return true if all the conn vertices belongs to the pool, false otherwise.
*/
bool IOCGNS::belongToPool(const bitpit::ConstProxyVector<long> & elementconn, const std::set<long> &pool){

    bool check = (elementconn.size() > 0);
    for(const long idV : elementconn){
        check = check && (pool.count(idV) > 0);
    }

    return check;
}

/*!
    Given a portion of cells referring to volmesh , and the list of vertices involved in it
    renumbered from 0 to nmaxVert-1, return a cgns map connectivity having as key the type of elements
    and as argument the compact vector of the connectivity.

    \param[in] cellIds  list of cells
    \param[in] mapToLocVert map from volmesh globla enumeration to local one starting from 0 up to max Number of vertices involved.
    \param[out] ncells map of number of cells for each type of elements.

    \return the Zone connectivity of a certain portion of volume mesh.

 */
std::map<int, std::vector<std::size_t> >
IOCGNS::getZoneConn(const livector1D& cellIds,
                    const std::unordered_map<long, int> & mapToLocVert,
                    std::map<int, std::size_t> &ncells)
{
    std::map<int, std::vector<std::size_t> > mm;
    ncells.clear();
    bitpit::PiercedVector<bitpit::Cell> & cells = getGeometry()->getCells();

    bitpit::ElementType btype;
    CGNS_ENUMT(ElementType_t) cgnstype;
    bool unsupported = false;
    std::vector<std::size_t> connConverted;

    for(long idC : cellIds){
        bitpit::Cell & cell = cells.at(idC);
        btype = cell.getType();
        switch(btype){
            case bitpit::ElementType::TETRA:
                cgnstype = CGNS_ENUMV(TETRA_4);
            break;
            case bitpit::ElementType::PYRAMID:
                cgnstype = CGNS_ENUMV(PYRA_5);
            break;
            case bitpit::ElementType::WEDGE:
                cgnstype = CGNS_ENUMV(PENTA_6);
            break;
            case bitpit::ElementType::HEXAHEDRON:
                cgnstype = CGNS_ENUMV(HEXA_8);
            break;
            default:
                unsupported = true;
            break;
        }

        if(unsupported){
            (*m_log)<<"WARNING IOCGNS:write() : Skipping a polyhedral element in volume mesh. Unsupported writign at the moment."<<std::endl;
            continue;
        }
        long * conn = cell.getConnect();
        int iconn = cell.getConnectSize();
        connConverted.resize(iconn);
        for(int i=0; i< iconn; ++i){
            connConverted[i] = std::size_t(mapToLocVert.at(conn[i]));
        }
        if(cgnstype == CGNS_ENUMV(PENTA_6)){
            std::swap(connConverted[1], connConverted[2]);
            std::swap(connConverted[4], connConverted[5]);
        }
        int tt = static_cast<int>(cgnstype);
        mm[tt].insert(mm[tt].end(), connConverted.begin(), connConverted.end());
        if(ncells.count(tt) < 1){
            ncells[tt] = 0;
        }
        ncells[tt] += 1;
    }

    return mm;
};

/*!
    Given a portion of cells referring to volmesh , and the list of vertices involved in it
    renumbered from 0 to nmaxVert-1, return a cgns map connectivity having as key the type of elements
    and as argument the compact vector of the connectivity.

    \param[in] cellIds  list of cells
    \param[in] mapToLocVert map from volmesh globla enumeration to local one starting from 0 up to max Number of vertices involved.
    \param[out] ncells map of number of cells for each type of elements.
    \param[out] surfCellGlobToLoc  map form surface cell id to local compact consecutive id.
    \return the Zone connectivity of a certain portion of volume mesh.

 */
std::map<int, std::vector<std::size_t> >
IOCGNS::getBCElementsConn(const livector1D& cellIds,
                          const std::unordered_map<long, int> & mapToLocVert,
                          std::map<int, std::size_t> &ncells,
                          std::unordered_map<long, int> & surfCellGlobToLoc)
{
    std::map<int, std::vector<std::size_t> > mm;
    std::map<int, std::vector<long> > idInsertion;
    ncells.clear();
    bitpit::PiercedVector<bitpit::Cell> & cells = getSurfaceBoundary()->getCells();

    bitpit::ElementType btype;
    CGNS_ENUMT(ElementType_t) cgnstype;
    bool unsupported = false;
    std::vector<std::size_t> connConverted;

    for(long idC : cellIds){
        bitpit::Cell & cell = cells.at(idC);
        btype = cell.getType();
        switch(btype){
            case bitpit::ElementType::TRIANGLE:
                cgnstype = CGNS_ENUMV(TRI_3);
            break;
            case bitpit::ElementType::QUAD:
                cgnstype = CGNS_ENUMV(QUAD_4);
            break;
            default:
                unsupported = true;
            break;
        }

        if(unsupported){
            (*m_log)<<"WARNING IOCGNS:write() : Skipping a polyhedral element in volume mesh. Unsupported writign at the moment."<<std::endl;
            continue;
        }
        long * conn = cell.getConnect();
        int iconn = cell.getConnectSize();
        connConverted.resize(iconn);
        for(int i=0; i< iconn; ++i){
            connConverted[i] = std::size_t(mapToLocVert.at(conn[i]));
        }
        int tt = static_cast<int>(cgnstype);
        mm[tt].insert(mm[tt].end(), connConverted.begin(), connConverted.end());
        idInsertion[tt].push_back(idC);
    }

    int counter = 0;
    for(auto &mapp: idInsertion){
        for(long & val: mapp.second){
            surfCellGlobToLoc.insert(std::make_pair(val, counter) );
            ++counter;
        }
        ncells[mapp.first] = mapp.second.size();
    }
    return mm;
};

/*!
    Unpack conns array of mixed elements.
    \param[in] connsArray of mixed conns as read by cgns.
    \param[in,out] patchVol mimmoObject to store volume elements.
    \param[in,out] surfElem map to store surface elements for boundary purposes (key ordered according to the element type)
    \param[in] idVertexOffset id offset for remapping local vertices.
    \param[in] idCellOffset id offset for remapping local cells.
    \param[in] PIDZoneVolume PID to be assigned to volumtric cells
    \param[in,out] idwork id counter for new cells added.
*/

void
IOCGNS::unpackMixedConns( const ivector1D & connsArray,
                          MimmoObject * patchVol,
                          std::unordered_map<long, std::vector<long>> surfElem,
                          const long & idVertexOffset,
                          const long & idCellOffset,
                          const long & PIDZoneVolume,
                          long & idwork)
{
    livector1D lConn;
    bitpit::ElementType btype;
    auto it=connsArray.begin(), itE=connsArray.end();

    while(it != itE){

        CGNS_ENUMT(ElementType_t) et = static_cast<CGNS_ENUMT(ElementType_t)>(*it);
        ++it;

        switch(et){

        case CGNS_ENUMV(TETRA_4):
        case CGNS_ENUMV(TETRA_10):
        case CGNS_ENUMV(TETRA_16):
        case CGNS_ENUMV(TETRA_20):
        case CGNS_ENUMV(TETRA_22):
        case CGNS_ENUMV(TETRA_34):
        case CGNS_ENUMV(TETRA_35):

            btype = bitpit::ElementType::TETRA;
            lConn.resize(4);
            for(long & val: lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            patchVol->addConnectedCell(lConn, btype, idwork+idCellOffset );
            patchVol->setPIDCell(idwork + idCellOffset,PIDZoneVolume);
            idwork++;
        break;

        case CGNS_ENUMV(PYRA_5):
        case CGNS_ENUMV(PYRA_13):
        case CGNS_ENUMV(PYRA_14):
        case CGNS_ENUMV(PYRA_21):
        case CGNS_ENUMV(PYRA_29):
        case CGNS_ENUMV(PYRA_30):
        case CGNS_ENUMV(PYRA_50):
        case CGNS_ENUMV(PYRA_55):

            btype = bitpit::ElementType::PYRAMID;
            lConn.resize(5);
            for(long & val : lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            patchVol->addConnectedCell(lConn, btype, idwork+idCellOffset );
            patchVol->setPIDCell(idwork+idCellOffset,PIDZoneVolume);
            idwork++;
        break;

        case CGNS_ENUMV(PENTA_6):
        case CGNS_ENUMV(PENTA_15):
        case CGNS_ENUMV(PENTA_18):
        case CGNS_ENUMV(PENTA_24):
        case CGNS_ENUMV(PENTA_38):
        case CGNS_ENUMV(PENTA_40):
        case CGNS_ENUMV(PENTA_33):
        case CGNS_ENUMV(PENTA_66):
        case CGNS_ENUMV(PENTA_75):

            btype = bitpit::ElementType::WEDGE;
            lConn.resize(6);
            for(long & val : lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            //remap in bitpit conn. TODO complete ref element mapper for connectivity.
            std::swap(lConn[1], lConn[2]);
            std::swap(lConn[4], lConn[5]);

            patchVol->addConnectedCell(lConn, btype, idwork+idCellOffset );
            patchVol->setPIDCell(idwork+idCellOffset,PIDZoneVolume);
            idwork++;

        break;

        case CGNS_ENUMV(HEXA_8):
        case CGNS_ENUMV(HEXA_20):
        case CGNS_ENUMV(HEXA_27):
        case CGNS_ENUMV(HEXA_32):
        case CGNS_ENUMV(HEXA_56):
        case CGNS_ENUMV(HEXA_64):
        case CGNS_ENUMV(HEXA_44):
        case CGNS_ENUMV(HEXA_98):
        case CGNS_ENUMV(HEXA_125):

            btype = bitpit::ElementType::HEXAHEDRON;
            lConn.resize(8);
            for(long & val : lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            patchVol->addConnectedCell(lConn, btype, idwork+idCellOffset );
            patchVol->setPIDCell(idwork+idCellOffset,PIDZoneVolume);
            idwork++;
        break;

        case CGNS_ENUMV(TRI_3):
        case CGNS_ENUMV(TRI_6):
        case CGNS_ENUMV(TRI_9):
        case CGNS_ENUMV(TRI_10):
        case CGNS_ENUMV(TRI_12):
        case CGNS_ENUMV(TRI_15):

            lConn.resize(3);
            for(long & val : lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            surfElem[idwork] = lConn;
            idwork++;
        break;

        case CGNS_ENUMV(QUAD_4):
        case CGNS_ENUMV(QUAD_8):
        case CGNS_ENUMV(QUAD_9):
        case CGNS_ENUMV(QUAD_12):
        case CGNS_ENUMV(QUAD_16):
        case CGNS_ENUMV(QUAD_25):

            lConn.resize(4);
            for(long & val : lConn){
                val = idVertexOffset + *it;
                ++it;
            }
            surfElem[idwork] = lConn;
            idwork++;
        break;

        case CGNS_ENUMV(NODE):
            ++it;
            idwork++;
        break;

        case CGNS_ENUMV(BAR_2):
            for(int j=0; j<2; ++j){
                ++it;
            }
            idwork++;
        break;
        case CGNS_ENUMV(BAR_3):
            for(int j=0; j<3; ++j){
                ++it;
            }
            idwork++;
        break;
        case CGNS_ENUMV(BAR_4):
            for(int j=0; j<4; ++j){
                ++it;
            }
            idwork++;
        break;
        case CGNS_ENUMV(BAR_5):
            for(int j=0; j<5; ++j){
                ++it;
            }
            idwork++;
        break;

        default:
            (*m_log)<< "error: "<< m_name << " found unrecognized CGNS element while reading. Impossible to absorb further mixed elements. "<<std::endl;
            throw std::runtime_error (m_name + " : found unrecognized CGNS element while reading. Impossible to absorb further mixed elements. ");
            return;
            break;
        }
    } //end while
}



/*!
    Dump class contents on a stream
    \param[in] stream output stream
*/
bool IOCGNS::dump(std::ostream &stream){
    try{

        getGeometry()->dump(stream);
        getSurfaceBoundary()->dump(stream);
        m_storedBC->dump(stream);

        // m_storedInfo->dump(stream); ITS USELESS FOR NOW!

    }catch(std::exception & ee){
        return false;
    }
    return true;
}

/*!
    Restore class contents from a stream
    \param[in] stream input stream
*/
bool IOCGNS::restore(std::istream &stream){

    try{

        m_volmesh->restore(stream);
        m_surfmesh->restore(stream);
        m_storedBC->restore(stream);
        // m_storedInfo->restore(stream); //it's useless for now.

    }catch(std::exception & ee){
        return false;
    }

    m_volmesh->buildPatchInfo();
    m_surfmesh->buildPatchInfo();

#if MIMMO_ENABLE_MPI
    m_volmesh->updatePointGhostExchangeInfo();
    m_surfmesh->updatePointGhostExchangeInfo();

#endif
    return true;
}

#if MIMMO_ENABLE_MPI

/*!
 * Makes rank 0 communicate m_storedBC info to all other procs.
 */
void IOCGNS::communicateAllProcsStoredBC(){

    if(m_rank == 0){

        //create char output data buffer and reverse data into it.
        bitpit::OBinaryStream dataBuffer_pidtolisttype;
        bitpit::OBinaryStream dataBuffer_pidtobc;
        bitpit::OBinaryStream dataBuffer_zonetobndpid;
        bitpit::OBinaryStream dataBuffer_bcpidnames;
        bitpit::OBinaryStream dataBuffer_zonepidnames;

        dataBuffer_pidtobc << m_storedBC->mcg_pidtobc;
        dataBuffer_zonetobndpid << m_storedBC->mcg_zonetobndpid;
        dataBuffer_pidtolisttype << m_storedBC->mcg_pidtolisttype;
        dataBuffer_bcpidnames << m_storedBC->mcg_bcpidnames;
        dataBuffer_zonepidnames << m_storedBC->mcg_zonepidnames;

        long dbs1 = dataBuffer_pidtobc.getSize();
        long dbs2 = dataBuffer_zonetobndpid.getSize();
        long dbs3 = dataBuffer_pidtolisttype.getSize();
        long dbs4 = dataBuffer_bcpidnames.getSize();
        long dbs5 = dataBuffer_zonepidnames.getSize();

        //Send data to all other procs
        for (int sendRank=1; sendRank<m_nprocs; sendRank++){
           MPI_Send(&dbs1, 1, MPI_LONG, sendRank, 100, m_communicator);
           MPI_Send(dataBuffer_pidtobc.data(), dataBuffer_pidtobc.getSize(), MPI_CHAR, sendRank, 110, m_communicator);
           MPI_Send(&dbs2, 1, MPI_LONG, sendRank, 200, m_communicator);
           MPI_Send(dataBuffer_zonetobndpid.data(), dataBuffer_zonetobndpid.getSize(), MPI_CHAR, sendRank, 210, m_communicator);
           MPI_Send(&dbs3, 1, MPI_LONG, sendRank, 300, m_communicator);
           MPI_Send(dataBuffer_pidtolisttype.data(), dataBuffer_pidtolisttype.getSize(), MPI_CHAR, sendRank, 310, m_communicator);
           MPI_Send(&dbs4, 1, MPI_LONG, sendRank, 400, m_communicator);
           MPI_Send(dataBuffer_bcpidnames.data(), dataBuffer_bcpidnames.getSize(), MPI_CHAR, sendRank, 410, m_communicator);
           MPI_Send(&dbs5, 1, MPI_LONG, sendRank, 500, m_communicator);
           MPI_Send(dataBuffer_zonepidnames.data(), dataBuffer_zonepidnames.getSize(), MPI_CHAR, sendRank, 510, m_communicator);

        }
        //hey 0, your job is done.
    }else{

        m_storedBC  = std::move(std::unique_ptr<BCCGNS>(new BCCGNS()));

        long dbs1,dbs2, dbs3,dbs4, dbs5;

        MPI_Recv(&dbs1, 1, MPI_LONG, 0, 100, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer_pidtobc(dbs1);
        MPI_Recv(dataBuffer_pidtobc.data(), dataBuffer_pidtobc.getSize(), MPI_CHAR, 0, 110, m_communicator, MPI_STATUS_IGNORE);
        dataBuffer_pidtobc >> m_storedBC->mcg_pidtobc;

        MPI_Recv(&dbs2, 1, MPI_LONG, 0, 200, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer_zonetobndpid(dbs2);
        MPI_Recv(dataBuffer_zonetobndpid.data(), dataBuffer_zonetobndpid.getSize(), MPI_CHAR, 0, 210, m_communicator, MPI_STATUS_IGNORE);
        dataBuffer_zonetobndpid >> m_storedBC->mcg_zonetobndpid;

        MPI_Recv(&dbs3, 1, MPI_LONG, 0, 300, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer_pidtolisttype(dbs3);
        MPI_Recv(dataBuffer_pidtolisttype.data(), dataBuffer_pidtolisttype.getSize(), MPI_CHAR, 0, 310, m_communicator, MPI_STATUS_IGNORE);
        dataBuffer_pidtolisttype >> m_storedBC->mcg_pidtobc;

        MPI_Recv(&dbs4, 1, MPI_LONG, 0, 400, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer_bcpidnames(dbs4);
        MPI_Recv(dataBuffer_bcpidnames.data(), dataBuffer_bcpidnames.getSize(), MPI_CHAR, 0, 410, m_communicator, MPI_STATUS_IGNORE);
        dataBuffer_bcpidnames >> m_storedBC->mcg_bcpidnames;

        MPI_Recv(&dbs5, 1, MPI_LONG, 0, 500, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer_zonepidnames(dbs5);
        MPI_Recv(dataBuffer_zonepidnames.data(), dataBuffer_zonepidnames.getSize(), MPI_CHAR, 0, 510, m_communicator, MPI_STATUS_IGNORE);
        dataBuffer_zonepidnames >> m_storedBC->mcg_zonepidnames;

    }

}

#endif

/*!
    Dump class contents on a stream
    \param[in] out output stream
*/
void BCCGNS::dump(std::ostream & out){

    int pidtolisttype_size = mcg_pidtolisttype.size();
    int pidtobc_size = mcg_pidtobc.size();
    int zonetobndpid_size = mcg_zonetobndpid.size();
    int bcpidnames_size = mcg_bcpidnames.size();
    int zonepidnames_size = mcg_zonepidnames.size();

    {
        bitpit::utils::binary::write(out, pidtobc_size);
        auto it = mcg_pidtobc.begin();
        for(int i=0; i<pidtobc_size; ++i){
            bitpit::utils::binary::write(out, it->first);
            bitpit::utils::binary::write(out, it->second);
            ++it;
        }
    }

    {
        bitpit::utils::binary::write(out, zonetobndpid_size);
        auto it = mcg_zonetobndpid.begin();
        for(int i=0; i<zonetobndpid_size; ++i){
            bitpit::utils::binary::write(out, it->first);
            bitpit::utils::binary::write(out, it->second);
            ++it;
        }
    }
    {
        bitpit::utils::binary::write(out, pidtolisttype_size);
        auto it = mcg_pidtolisttype.begin();
        for(int i=0; i<pidtolisttype_size; ++i){
            bitpit::utils::binary::write(out, it->first);
            bitpit::utils::binary::write(out, it->second);
            ++it;
        }
    }
    {
        bitpit::utils::binary::write(out, bcpidnames_size);
        auto it = mcg_bcpidnames.begin();
        for(int i=0; i<bcpidnames_size; ++i){
            bitpit::utils::binary::write(out, it->first);
            bitpit::utils::binary::write(out, it->second);
            ++it;
        }
    }
    {
        bitpit::utils::binary::write(out, zonepidnames_size);
        auto it = mcg_zonepidnames.begin();
        for(int i=0; i<zonepidnames_size; ++i){
            bitpit::utils::binary::write(out, it->first);
            bitpit::utils::binary::write(out, it->second);
            ++it;
        }
    }

}

/*!
    Restore class contents from a stream
    \param[in] in input stream
*/
void BCCGNS::restore(std::istream & in){
    int pidtobc_size, zonetobndpid_size, pidtolisttype_size, bcpidnames_size, zonepidnames_size;
    clear();

    {
        bitpit::utils::binary::read(in, pidtobc_size);
        int key;
        M_CG_BCType_t val;

        for(int i=0; i<pidtobc_size; ++i){
            bitpit::utils::binary::read(in, key);
            bitpit::utils::binary::read(in, val);
            mcg_pidtobc.insert(std::make_pair(key,val));
        }
    }


    {
        bitpit::utils::binary::read(in, zonetobndpid_size);
        int key;
        std::vector<int> val;

        for(int i=0; i<zonetobndpid_size; ++i){
            bitpit::utils::binary::read(in, key);
            bitpit::utils::binary::read(in, val);
            mcg_zonetobndpid.insert(std::make_pair(key,val));
        }
    }

    {
        bitpit::utils::binary::read(in, pidtolisttype_size);
        int key;
        int val;

        for(int i=0; i<pidtolisttype_size; ++i){
            bitpit::utils::binary::read(in, key);
            bitpit::utils::binary::read(in, val);
            mcg_pidtolisttype.insert(std::make_pair(key,val));
        }
    }

    {
        bitpit::utils::binary::read(in, bcpidnames_size);
        int key;
        std::string val;

        for(int i=0; i<bcpidnames_size; ++i){
            bitpit::utils::binary::read(in, key);
            bitpit::utils::binary::read(in, val);
            mcg_bcpidnames.insert(std::make_pair(key,val));
        }
    }

    {
        bitpit::utils::binary::read(in, zonepidnames_size);
        int key;
        std::string val;

        for(int i=0; i<zonepidnames_size; ++i){
            bitpit::utils::binary::read(in, key);
            bitpit::utils::binary::read(in, val);
            mcg_zonepidnames.insert(std::make_pair(key,val));
        }
    }


}




};
