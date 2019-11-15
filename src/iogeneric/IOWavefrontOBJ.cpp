/*
    Get a License Header here- Temporary POC for Click-Ins
*/
#include "IOWavefrontOBJ.hpp"
#include <algorithm>

namespace mimmo{


/*!
 Constructor
*/
WavefrontObjData::WavefrontObjData(){
    materialfile="";
    materials.setName("Materials");
    smoothids.setName("SmooothingGroupIds");
    materials.setDataLocation(MPVLocation::CELL);
    smoothids.setDataLocation(MPVLocation::CELL);

}

/*!
    Swap operator. Swap content with a twin class.
    \param[in] other WavefrontObjData object.
*/
void
WavefrontObjData::swap(WavefrontObjData & x) noexcept{
    materials.swap(x.materials);
    smoothids.swap(x.smoothids);
    std::swap(materialsList, x.materialsList);
    std::swap(smoothidsList, x.smoothidsList);
    std::swap(materialfile, x.materialfile);

}


/*!
    Synchronize the lists with the inner data currently present in the class.
*/
void
WavefrontObjData::syncListsOnData(){
    //materials
    materialsList.clear();
    long id = 0;
    for(auto it = materials.begin(); it!= materials.end(); ++it){
        if (materialsList.count(*it) ==0){
            materialsList.insert({{*it, id}});
            ++id;
        }
    }

    //smoothids
    smoothidsList.clear();
    for(auto it = smoothids.begin(); it!= smoothids.end(); ++it){
        smoothidsList.insert(*it);
    }
}

/*!
    Dump class contents to an output stream in binary format
    \param[in] out binary output stream
*/
void
WavefrontObjData::dump(std::ostream & out){
    int location  = static_cast<int>(MPVLocation::CELL);
    std::string name;
    std::size_t totSize, locSize;

    //dump materials
    name = materials.getName();
    totSize = materials.size();

    bitpit::utils::binary::write(out, location);
    bitpit::utils::binary::write(out, name);
    bitpit::utils::binary::write(out, totSize);

    for(auto it=materials.begin(); it!=materials.end(); ++it){
        bitpit::utils::binary::write(out, it.getId());
        bitpit::utils::binary::write(out, *it);
    }

    //dump smoothids
    name = smoothids.getName();
    totSize = smoothids.size();

    bitpit::utils::binary::write(out, location);
    bitpit::utils::binary::write(out, name);
    bitpit::utils::binary::write(out, totSize);

    for(auto it=smoothids.begin(); it!=smoothids.end(); ++it){
        bitpit::utils::binary::write(out, it.getId());
        bitpit::utils::binary::write(out, *it);
    }

    // dump material file
    bitpit::utils::binary::write(out, materialfile);
}
/*!
    Restore class contents from an input stream in binary format
    \param[in] in binary input stream
*/
void
WavefrontObjData::restore(std::istream & in){

    int location;
    std::string name;
    std::size_t totSize, locSize;
    long id;

    materials.clear();
    //restore materials
    {
        bitpit::utils::binary::read(in, location);
        bitpit::utils::binary::read(in, name);
        bitpit::utils::binary::read(in, totSize);

        materials.setName(name);
        materials.setDataLocation(static_cast<MPVLocation>(location));
        materials.reserve(totSize);

        std::string data;
        for(int i=0; i<totSize; ++i){
            bitpit::utils::binary::read(in, id);
            bitpit::utils::binary::read(in, data);
            materials.insert(id, data);
        }
    }

    smoothids.clear();
    //restore smoothids
    {
        bitpit::utils::binary::read(in, location);
        bitpit::utils::binary::read(in, name);
        bitpit::utils::binary::read(in, totSize);

        smoothids.setName(name);
        smoothids.setDataLocation(static_cast<MPVLocation>(location));
        smoothids.reserve(totSize);

        long data;
        for(int i=0; i<totSize; ++i){
            bitpit::utils::binary::read(in, id);
            bitpit::utils::binary::read(in, data);
            smoothids.insert(id, data);
        }
    }

    //RESTORE material file
    bitpit::utils::binary::read(in, materialfile);

    // resync lists
    syncListsOnData();

}


/*!
    Constructor
    \param[in] mode working mode of the class - see IOWavefrontOBJ::IOMode enum
*/
IOWavefrontOBJ::IOWavefrontOBJ(IOWavefrontOBJ::IOMode mode){
    m_name = "mimmo.IOWavefrontOBJ";
    m_mode = mode;
    m_dir = ".";
    m_filename = "iowavefrontobj";
    m_resume = false;
    m_tol = 1.0E-8;
    m_txttol = 1.0E-8;
    m_skiptexture = false;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
IOWavefrontOBJ::IOWavefrontOBJ(const bitpit::Config::Section & rootXML){
    m_name = "mimmo.IOWavefrontOBJ";
    m_mode =  IOMode::READ;
    m_dir = ".";
    m_filename = "iowavefrontobj";
    m_resume = false;
    m_skiptexture = false;
    m_tol = 1.0E-8;
    m_txttol = 1.0E-8;

    std::string fallback_name = "ClassNONE";
    std::string fallback_mode = "0";
    std::string input = rootXML.get("ClassName", fallback_name);
    std::string mode = rootXML.get("IOMode", fallback_mode);

    input = bitpit::utils::string::trim(input);
    mode  = bitpit::utils::string::trim(mode);

    if(input == "mimmo.IOWavefrontOBJ"){
        std::stringstream ss(mode);
        int locmode;
        ss >> locmode;
        locmode = std::max(0, std::min(3, locmode));
        m_mode = static_cast<IOMode>(locmode);
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
    Destructor
*/
IOWavefrontOBJ::~IOWavefrontOBJ(){};

/*!
    Swap utility
    \param[in] x other object to swap from
*/
void IOWavefrontOBJ::swap(IOWavefrontOBJ & x) noexcept{
    std::swap(m_mode, x.m_mode);
    std::swap(m_dir, x.m_dir);
    std::swap(m_filename, x.m_filename);
    std::swap(m_resume, x.m_resume);
    std::swap(m_tol, x.m_tol);
    std::swap(m_tol, x.m_txttol);
    std::swap(m_skiptexture, x.m_skiptexture);

    std::swap(m_intPatch, x.m_intPatch);
    std::swap(m_intTexture, x.m_intTexture);
    std::swap(m_intData, x.m_intData);
    std::swap(m_extTexture, x.m_extTexture);
    std::swap(m_extData, x.m_extData);
    BaseManipulation::swap(x);
}

/*!
    Building class ports
*/
void IOWavefrontOBJ::buildPorts(){
    bool built = true;
    bool mandatoryWrite;
    switch(m_mode){
        case IOMode::WRITE:
        case IOMode::DUMP:
            mandatoryWrite = true;
            break;
        default: //all other cases, read, restore
            mandatoryWrite = false;
            break;
    }

    built = (built && createPortIn<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setGeometry, M_GEOM, mandatoryWrite));
    built = (built && createPortIn<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setTexture, X_TEXTURE));
    built = (built && createPortIn<WavefrontObjData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setData, X_WDATA));
    built = (built && createPortIn<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setMaterials, M_STRINGFIELD, false, 2));
    built = (built && createPortIn<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setSmoothIds, M_LONGFIELD, false, 2));

    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getGeometry, M_GEOM));
    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getTexture, X_TEXTURE));
    built = (built && createPortOut<WavefrontObjData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getData, X_WDATA));
    built = (built && createPortOut<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getMaterials, M_STRINGFIELD));
    built = (built && createPortOut<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getSmoothIds, M_LONGFIELD));
    m_arePortsBuilt = built;

}

/*!
    \return working mode of the class, see IOWavefrontOBJ::IOMode enum.
*/
IOWavefrontOBJ::IOMode IOWavefrontOBJ::whichMode(){
    return m_mode;
};

/*!
    \return working mode of the class as integer, see IOWavefrontOBJ::IOMode enum.
*/
int IOWavefrontOBJ::whichModeInt(){
    return static_cast<int>(whichMode());
};

/*!
    Return pointer to MimmoObject geometry.
    - In Read mode : return pointer to the what resulted from .obj file reading
    - in Write mode: return pointer to what is set with setGeometry method.
    \return MimmoObject mesh
*/
MimmoObject*   IOWavefrontOBJ::getGeometry(){
    if(whichModeInt()> 1) return m_geometry;
    else                  return m_intPatch.get();
};

/*!
    Return pointer to MimmoObject storing the texture. Texture is treated as
    an external surface geometry who shares the same cell ids then polygonal geoemetry
    which refers. The cell type is identical between mesh
    - In Read/Restore mode : return pointer to Texture read from .obj file
    - in Write/Dump mode: return pointer to what is set with setTexture method.
    \return texture data
*/
MimmoObject*   IOWavefrontOBJ::getTexture(){
    if(m_skiptexture)   return nullptr;
    if(whichModeInt()> 1) return m_extTexture;
    else                  return m_intTexture.get();

};

/*!
    Return pointer to WavefrontObjData (materials, group, smoothids) attached to polygonal mesh.
    Data are referred always to cell.
    - In Read/Restore mode : return pointer to data read from .obj file
    - in Write/Dump mode: return pointer to what is set with setData method.
    \return data attached to the mesh
*/
WavefrontObjData*   IOWavefrontOBJ::getData(){
    if(whichModeInt()> 1) return m_extData;
    else                  return m_intData.get();

};

/*!
    Return pointer to WavefrontObjData materials attached to polygonal mesh.
    Data are referred always to cell.
    - In Read/Restore mode : return pointer to data read from .obj file
    - in Write/Dump mode: return pointer to what is set with setData method.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<std::string>*   IOWavefrontOBJ::getMaterials(){
    if(whichModeInt()> 1) return &(m_extData->materials);
    else                  return &(m_intData->materials);
};

/*!
    Return pointer to WavefrontObjData smoothids attached to polygonal mesh.
    Data are referred always to cell.
    - In Read/Restore mode : return pointer to data read from .obj file
    - in Write/Dump mode: return pointer to what is set with setData method.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<long>*   IOWavefrontOBJ::getSmoothIds(){
    if(whichModeInt()> 1) return &(m_extData->smoothids);
    else                  return &(m_intData->smoothids);
};

/*!
    Return the map with PID and name of the object which the mesh is subdivided in.
    In read/restore mode, this subdivision is extrapolated from file obj, in write/dump mode this refers to
    subdivision hold by the geometry set with setGeometry method.
    \return pid/name mesh internal subdivision.
*/
std::unordered_map<long, std::string> IOWavefrontOBJ::getSubParts(){
    if(!getGeometry()) return std::unordered_map<long, std::string>();
    return getGeometry()->getPIDTypeListWNames();
}

/*!
    \return true if the class is forcibly skipping reading or writing texture
    associated to the mesh.
*/
bool    IOWavefrontOBJ::isSkippingTexture(){
    return m_skiptexture;
}
/*!
    Set the geometry meant to be written.  Does nothing in read/restore mode.
    \param[in] geo MimmoObject surface mesh of type 1
*/
void    IOWavefrontOBJ::setGeometry(MimmoObject * geo){
    if(!geo || whichModeInt()<2) return;
    if(geo->getType() != 1) return;
    m_geometry = geo;
}

/*!
    Set texture mesh related to surface polygonal mesh. They should share
    the same cell-ids.
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] texture mesh.
*/
void    IOWavefrontOBJ::setTexture(MimmoObject* texture){
    if(!texture || whichModeInt()<2) return;
    m_extTexture = texture;
};

/*!
    Set data attached to surface polygonal mesh, i.e:
    - materials associated to cells
    - cell group labels associated to cells
    - smoothing group ids associated to cells
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] data.
*/
void    IOWavefrontOBJ::setData(WavefrontObjData* data){
    if(!data || whichModeInt()<2) return;
    m_extData = data;
};

/*!
    Set materials data attached to surface polygonal mesh, i.e the materials associated to cells
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] data.
*/
void    IOWavefrontOBJ::setMaterials(MimmoPiercedVector<std::string>* data){
    if(!data || whichModeInt()<2) return;
    // I can use extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData == nullptr){
        *(m_log)<<m_name<<" : external data structure not yet linked before set materials calling."<<std::endl;
        throw std::runtime_error(m_name + " : external data structure not yet linked before set materials calling.");
    }
    m_extData->materials = *data;
};

/*!
    Set materials data attached to surface polygonal mesh, i.e the smoothing group ids associated to cells
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] data.
*/
void    IOWavefrontOBJ::setSmoothIds(MimmoPiercedVector<long>* data){
    if(!data || whichModeInt()<2) return;
    // I can use extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData == nullptr){
        *(m_log)<<m_name<<" : external data structure not yet linked before set smoothids calling."<<std::endl;
        throw std::runtime_error(m_name + " : external data structure not yet linked before set smoothids calling.");
    }
    m_extData->smoothids = *data;
};

/*!
    Set the directory where the read file is searched
    or the written file has to be places.
    \param[in] dir path to reference directory for I/O purposes
*/
void    IOWavefrontOBJ::setDir(const std::string & pathdir){
    m_dir= pathdir;
}

/*!
    Set the name of the file to be read or written.
    No format tag .obj has to be specified.
    \param[in] name of the file
*/
void    IOWavefrontOBJ::setFilename(const std::string & name){
    m_filename= name;
}

/*!
    If set to true print a resume file of the data referenced by the mesh.
    \param[in] print true/false
*/
void    IOWavefrontOBJ::printResumeFile(bool print){
    m_resume= print;
}

/*!
    Skip texture reading/writing no matter the mode of the class of texture
    object linked through setTexture.
    \param[in] skip if true forcibly skip IO of texture.
*/
void    IOWavefrontOBJ::skipTexture(bool skip){
    m_skiptexture= skip;
}


/*!
    set the geometric tolerance for duplicated vertices collapsing
    \param[in] tolerance
*/
void    IOWavefrontOBJ::setGeomTolerance(double tolerance){
    m_tol= std::max(std::numeric_limits<double>::min(), tolerance);
}

/*!
    Set the Texture geometric tolerance for duplicated vertices collapsing
    \param[in] tolerance
*/
void    IOWavefrontOBJ::setTextureTolerance(double tolerance){
    m_txttol= std::max(std::numeric_limits<double>::min(), tolerance);
}

/*!
    Absorb class parameters from a section xml
    \param[in] slotXML xml section
    \param[in] name unused string
*/
void    IOWavefrontOBJ::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name ){
    BITPIT_UNUSED(name);
    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("Dir")){
        input = slotXML.get("Dir");
        input = bitpit::utils::string::trim(input);
        if(input.empty())    input = "./";
        setDir(input);
    };


    if(slotXML.hasOption("Filename")){
        input = slotXML.get("Filename");
        input = bitpit::utils::string::trim(input);
        if(input.empty())    input = "iowavefrontobj";
        setFilename(input);
    };

    if(slotXML.hasOption("PrintResumeFile")){
        input = slotXML.get("PrintResumeFile");
        bool value = false;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        printResumeFile(value);
    };

    if(slotXML.hasOption("SkipTexture")){
        input = slotXML.get("SkipTexture");
        bool value = false;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        skipTexture(value);
    };

    if(slotXML.hasOption("GeomTolerance")){
        input = slotXML.get("GeomTolerance");
        double value = 1.0E-8;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setGeomTolerance(value);
    };

    if(slotXML.hasOption("TextureTolerance")){
        input = slotXML.get("TextureTolerance");
        double value = 1.0E-8;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setTextureTolerance(value);
    };

}

/*!
    Flush class parameters to a section xml
    \param[in] slotXML xml section to be filled
    \param[in] name unused string
*/
void    IOWavefrontOBJ::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("IOMode", std::to_string(static_cast<int>(m_mode)));
    slotXML.set("Dir", m_dir);
    slotXML.set("Filename", m_filename);
    slotXML.set("PrintResumeFile", std::to_string(m_resume));
    slotXML.set("GeomTolerance", std::to_string(m_tol));
    slotXML.set("TextureTolerance", std::to_string(m_txttol));
    slotXML.set("SkipTexture", std::to_string(m_skiptexture));
}


/*!
    Class workflow execution
*/
void   IOWavefrontOBJ::execute(){
    switch(m_mode) {
        case IOMode::READ:
            //instantiation of a brand new MimmoObject to absorb grid
            m_intPatch = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
            m_intTexture = nullptr;
            if(!m_skiptexture){
                m_intTexture = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            }
            read(m_dir+"/"+m_filename+".obj");
            break;
        case IOMode::RESTORE:
        {
            std::string filedump = m_dir+"/"+m_filename;
// #if MIMMO_ENABLE_MPI
//             bitpit::IBinaryArchive binaryReader(filedump,"dump", m_rank);
// #else
            bitpit::IBinaryArchive binaryReader(filedump,"dump");
//#endif
            //instantiation of a brand new MimmoObject to absorb grid
            m_intPatch = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
            m_intTexture = nullptr;
            if(!m_skiptexture){
                m_intTexture = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            }
            restore(binaryReader.getStream());
            binaryReader.close();
        }
            break;
        case IOMode::WRITE:
            write(m_dir+"/"+m_filename+".obj");
            break;

        case IOMode::DUMP:
        {
            int archiveVersion = 1;
            std::string header(m_name);
            std::string filedump = m_dir+"/"+m_filename;
// #if MIMMO_ENABLE_MPI
//             bitpit::OBinaryArchive binaryWriter(filedump, "dump", archiveVersion, header, m_rank);
// #else
            bitpit::OBinaryArchive binaryWriter(filedump, "dump", archiveVersion, header);
//#endif
            dump(binaryWriter.getStream());
            binaryWriter.close();
        }
            break;
        default:
            //never been reached;
            break;
    }


    if(m_resume){
        writeResumeFile();
    }
}

/*!
    Read mesh and data attached from file obj. In mpi version, read through the 0-rank only
*/
void IOWavefrontOBJ::read(const std::string & filename){

    std::ifstream in(filename);
    if(in.is_open()){

        //search and store the material file fullpath.
        m_intData->materialfile = searchMaterialFile(in);
        //search and store the stream position of the sub-objects
        std::vector<std::streampos> objectPositions;
        std::vector<std::string> objectNames;
        long nVertTot, nCellTot;
        searchObjectPosition(in, objectPositions, objectNames, nVertTot, nCellTot);

        //reserve stuff
        m_intPatch->getVertices().reserve(nVertTot);
        m_intPatch->getCells().reserve(nCellTot);
        if(!m_skiptexture){
            m_intTexture->getVertices().reserve(nVertTot);
            m_intTexture->getCells().reserve(nCellTot);
        }
        m_intData->materials.reserve(nCellTot);
        m_intData->smoothids.reserve(nCellTot);

        //compile safely data options
        m_intData->materials.setGeometry(m_intPatch.get());
        m_intData->smoothids.setGeometry(m_intPatch.get());
        m_intData->materials.setDataLocation(MPVLocation::CELL);
        m_intData->smoothids.setDataLocation(MPVLocation::CELL);
        m_intData->materials.setName("Materials");
        m_intData->smoothids.setName("SmooothingGroupIds");


        long pidObject(0);
        long vOffset(1), vTxtOffset(1), cOffset(1);
        //read file object by object
        for(std::streampos & pos : objectPositions){
            (*m_log)<<m_name<<" : reading object "<<objectNames[pidObject]<<"...";
            readObjectData(in, pos, pidObject,vOffset, vTxtOffset, cOffset);
            m_intPatch->setPIDName(pidObject, objectNames[pidObject]);
            ++pidObject;
            (*m_log)<<"done"<<std::endl;

        }

        in.close();

    }else{
        *(m_log)<<m_name<<" : impossible to read from obj file "<<filename<<std::endl;
        throw std::runtime_error("IOWavefrontOBJ::read(), impossible reading obj file");
    }



    //shrink to fit all the reserve stuffs;
    m_intPatch->getVertices().shrinkToFit();
    m_intPatch->getCells().shrinkToFit();
    if(!m_skiptexture){
        m_intTexture->getVertices().shrinkToFit();
        m_intTexture->getCells().shrinkToFit();
    }
    m_intData->materials.shrinkToFit();
    m_intData->smoothids.shrinkToFit();

    // sync the lists of intData
    m_intData->syncListsOnData();

    //collapse duplicated vertices if any
    m_intPatch->getPatch()->setTol(m_tol);
    m_intPatch->getPatch()->deleteCoincidentVertices();
    m_intPatch->getPatch()->deleteOrphanVertices();
    //resyncPID
    m_intPatch->resyncPID();

    if(!m_skiptexture){
        m_intTexture->getPatch()->setTol(m_txttol);
        if(m_intTexture->getNVertices()> 0) {
            m_intTexture->getPatch()->deleteCoincidentVertices();
            m_intTexture->getPatch()->deleteOrphanVertices();
        }
        //resyncPID
        m_intTexture->resyncPID();
    }
}

/*!
    Write mesh and data attached to file obj
*/
void IOWavefrontOBJ::write(const std::string & filename){
    if(!m_geometry){
        (*m_log)<<"WARNING: no geometry linked in "<<m_name<<". Nothing to write on file "<<filename<<std::endl;
        return;
    }
    if(!m_geometry->areAdjacenciesBuilt())  m_geometry->buildAdjacencies();

    std::ofstream out(filename);
    if(out.is_open()){

        out<<"# Optimad's mimmo : "<<m_name<<" OBJ file"<<std::endl;
        out<<"# www.github.com/optimad/mimmo"<<std::endl;
        std::string materialfile= "UnknownFile.mtl";
        if(m_extData) materialfile = m_extData->materialfile;
        out<<"mtllib "<< materialfile<<std::endl;

        auto mapParts = getSubParts();
        long vOffset(1), vTxtOffset(1), cOffset(1);
        long refPid;
        //write object by object
        for(auto & entryPart : mapParts){
            //write header object;
            out<<"o "<<entryPart.second<<std::endl;

            refPid = entryPart.first;
            //select list of vertices and cells referring to this pid.
            livector1D cellList = m_geometry->extractPIDCells(refPid, true);
            livector1D vertList = m_geometry->getVertexFromCellList(cellList);
            livector1D txtVertList;
            if(m_extTexture && !m_skiptexture){
                txtVertList = m_extTexture->getVertexFromCellList(cellList);
            }
            (*m_log)<<m_name<<" : writing object "<<entryPart.second<<"...";
            writeObjectData(out,vertList,txtVertList, cellList, vOffset, vTxtOffset, cOffset);
            (*m_log)<<"done"<<std::endl;

        }

        out.close();

    }else{
        *(m_log)<<m_name<<" : impossible to write obj file "<<filename<<std::endl;
        throw std::runtime_error("IOWavefrontOBJ::write(), impossible writing obj file");
    }
}

/*!
    Read mesh and data attached from class own dump file
*/
void IOWavefrontOBJ::restore(std::istream & in){
    // m_intPatch, m_intTexture and m_intData are supposed to be
    //initialized

    bitpit::utils::binary::read(in, m_tol);
    bitpit::utils::binary::read(in, m_txttol);
    bitpit::utils::binary::read(in, m_skiptexture);

    bool geoMark, txtMark, dataMark;

    bitpit::utils::binary::read(in, geoMark);
    if(geoMark){
        m_intPatch->restore(in);
    }else{
        m_intPatch = nullptr;
    }

    bitpit::utils::binary::read(in, txtMark);
    if(txtMark){
        m_intTexture->restore(in);
    }else{
        m_intTexture = nullptr;
    }

    bitpit::utils::binary::read(in, dataMark);
    if(dataMark){
        m_intData->restore(in);
        //attach polygonal geometry to fields
        m_intData->materials.setGeometry(m_intPatch.get());
        m_intData->smoothids.setGeometry(m_intPatch.get());

    }else{
        m_intData = nullptr;
    }


}

/*!
    write mesh and data attached to class own dump file
*/
void IOWavefrontOBJ::dump(std::ostream & out){

    bitpit::utils::binary::write(out, m_tol);
    bitpit::utils::binary::write(out, m_txttol);
    bitpit::utils::binary::write(out, m_skiptexture);

    bool geoMark = (m_geometry != nullptr);
    bool txtMark = (m_extTexture != nullptr && !m_skiptexture);
    bool dataMark = (m_extData != nullptr);

    bitpit::utils::binary::write(out, geoMark);
    if(geoMark) m_geometry->dump(out);

    bitpit::utils::binary::write(out, txtMark);
    if(txtMark) m_extTexture->dump(out);

    bitpit::utils::binary::write(out, dataMark);
    if(dataMark) m_extData->dump(out);

}

/*!
    Print the resume file
*/
void    IOWavefrontOBJ::writeResumeFile(){
    std::ofstream out;
    std::string path = m_outputPlot+"/"+m_filename+"_RESUME.dat";
    out.open(path);
    if(out.is_open()){
        out<< "#IO log for mimmo: "<<m_name <<" class execution"<<std::endl;
        out<< "#"<<std::endl;
        out<< "#Reference I/O file: "<<m_filename<<std::endl;
        out<< "#"<<std::endl;
        out<< "#"<<std::endl;
        out<< "#Polygonal mesh info:"<<std::endl;
        if(getGeometry()){
            out<< "#    N vertices:     "<< getGeometry()->getNVertices()<<std::endl;
            out<< "#    N cells:        "<< getGeometry()->getNCells()<<std::endl;
            auto pidlist = getGeometry()->getPIDTypeListWNames();
            out<< "#    N objects(PID): "<< pidlist.size()<<std::endl;
            out<< "#    Object list:    "<<std::endl;
            for(auto & entry : pidlist){
                out<< "#        "<<entry.first<<" - "<< entry.second<<std::endl;
            }
        }
        out<< "#"<<std::endl;
        out<< "#"<<std::endl;
        out<< "#Texture mesh info:"<<std::endl;
        if(getTexture()){
            out<< "#    N vertices:     "<< getTexture()->getNVertices()<<std::endl;
            out<< "#    N cells:        "<< getTexture()->getNCells()<<std::endl;
        }
        out<< "#"<<std::endl;
        out<< "#"<<std::endl;
        out<< "#Data mesh info:"<<std::endl;
        if(getData()){
            getData()->syncListsOnData();
            std::map<long, std::string> locmap;
            for(auto & entry : getData()->materialsList){
                locmap.insert({{entry.second, entry.first}});
            }

            out<< "#    N Materials:     "<< locmap.size()<<std::endl;
            out<< "#    Material List:   "<<std::endl;
            for(auto & entry : locmap){
                out<< "#        "<<entry.first<<" - "<<entry.second<<std::endl;
            }
            out<< "#"<<std::endl;

            out<< "#    N SmoothGroups:    "<< getData()->smoothidsList.size()<<std::endl;
        }

        out.close();
    }else{
        (*m_log)<<"WARNING in "<<m_name<<" : not able to write Resume File. Aborting..."<<std::endl;
    }

}

/*!
    Search for material file name marked with entry key "mtllib".
    If none, return an empty string
    \return materials file name
*/
std::string IOWavefrontOBJ::searchMaterialFile(std::ifstream & in)
{
    if(!in.good()){
        in.clear();
        in.seekg(0);
    };

    std::string line,key, result("");
    bool breakloop = false;
    while(in && !breakloop){
        std::getline(in, line);
        line = bitpit::utils::string::trim(line);
        std::size_t fsplit = line.find(' ');
        key = line.substr(0, fsplit);
        if(key == "mtllib"){
            result = line.substr(fsplit);
            result = bitpit::utils::string::trim(result);
            breakloop = true;
        }
    }
    //restore stream to its beggining;
    in.clear();
    in.seekg(0);
    return result;
}

/*!
    While reading, scan the obj file to get the number of subparts, their stream positions
    and names eventually. Restore the stream to file begin after searching.
    \param[in] in reading obj file stream
    \param[in,out] mapPos map pid/streampos to be filled
    \param[in,out] mapNames map pid/names to be filled
    \param[out] nVertTot total vertices in the obj file
    \param[out] nCellTot total cells in the obj file
*/
void IOWavefrontOBJ::searchObjectPosition(std::ifstream & in,
                          std::vector<std::streampos> & mapPos,
                          std::vector<std::string>& mapNames,
                          long &nVertTot, long &nCellTot)
{
    nVertTot = 0;
    nCellTot = 0;
    if(!in.good()){
        in.clear();
        in.seekg(0);
    };

    std::string line,nn,key;
    while(in){
        std::getline(in, line);
        line = bitpit::utils::string::trim(line);
        std::size_t fsplit = line.find(' ');
        key = line.substr(0, fsplit);
        if(key == "o"){
            mapPos.push_back(in.tellg());
            nn = line.substr(fsplit);
            mapNames.push_back(bitpit::utils::string::trim(nn));
        }
        if(key == "v")  ++nVertTot;
        if(key == "f")  ++nCellTot;

    }
    mapPos.shrink_to_fit();
    mapNames.shrink_to_fit();
    //restore stream to its beggining;
    in.clear();
    in.seekg(0);
}

/*!
    Get the stream to the begin of a Object section and start absorb its data, i.e.
    - v: vertices coordinates
    - vt: texture coordinates (if any)
    - vn: normals referred to vertices (if any)
    - f : facet connectivity
    - usemtl : all the cells after the keyword will be assigned to the material specified by mtllib file (if any)
    - g: all the cells after the keyword will be assigned to cell group specified by g (if any)
    - s: all the cells after the keyword will be assigned to smooth group specified by s (if any)

    - vp, vertex in free form statement, is ignored. A warning will be raised on logger.
    - l, line connectivity, is ignored. A warning will be raised on logger.
    - p, point connectivity, is ignored. A warning will be raised on logger.

    all other keywords are ignored.

    Data are directly stored in m_intPatch (for main mesh), m_intTexture(for texture) and
    m_intData(for attached WavefrontObjData datatypes)
    \param[in] in reading stream
    \param[in] begObjectStream position of the stream to start absorbing data
    \param[in] PID part identifier assigned to the current object.
    \param[in, out] vOffset i/o current id of the next-to-be inserted mesh vertex.
    \param[in, out] vTxtOffset i/o current id of the next-to-be inserted texture vertex.
    \param[in, out] cOffset i/o current id of the next-to-be inserted mesh cell (same for mesh and texture)
*/
void IOWavefrontOBJ::readObjectData(std::ifstream & in, const std::streampos &begObjectStream, const long &PID,
                     long &vOffset, long &vTxtOffset, long &cOffset)
{
    if (!in.good()){
        in.clear();
    }
    in.seekg(begObjectStream);
    std::string line, key(" ");

    std::string activeMaterial = "NONE", stringsmooth;
    long activeSmooth = 0;
    std::vector<long> locConn, txtConn;

    while(in.good() && key!="o"){
        std::getline(in, line);
        line = bitpit::utils::string::trim(line);
        std::stringstream ss(line);
        ss >> key;
        darray3E temp;
        std::string label;
        std::vector<std::string> tempstring;
        switch( convertKeyEntryToInt(key) ){
            case 0: //vertex v
                temp.fill(0.0);
                ss>>temp[0]>>temp[1]>>temp[2];
                m_intPatch->addVertex(temp, vOffset);
                ++vOffset;
                break;
            case 1: //texture vt
                if(!m_skiptexture){
                    temp.fill(0.0);
                    int c = 0;
                    while(ss.good() && c<3){
                        ss>>temp[c];
                        ++c;
                    }
                    m_intTexture->addVertex(temp, vTxtOffset);
                    ++vTxtOffset;
                }
                break;
            case 2: //normals vn, I don't really need it
                break;
            case 3: //facet cell fn
                {
                    tempstring.clear();
                    while(ss.good()){
                        ss >> label;
                        tempstring.push_back(label);
                    }

                    if(tempstring.size() < 3) break;
                    locConn.resize(tempstring.size());
                    txtConn.resize(tempstring.size());

                    int counter(0);
                    bool textureOn = false;
                    for(std::string & str: tempstring){
                        str = str.substr(0,str.find("//"));
                        std::replace(str.begin(), str.end(), '/', ' ');
                        std::stringstream ss(str);
                        ss>>locConn[counter];

                        //handling negative vertex indexing rule of obj file
                        if(locConn[counter] < 0) locConn[counter] += vOffset;

                        if(ss.good() && !m_skiptexture){
                            textureOn=true;
                            ss>>txtConn[counter];
                            if(txtConn[counter] < 0) txtConn[counter] += vTxtOffset;
                        }
                        //increment counter;
                        ++counter;
                    }

                    //adding cell desuming cell type from locConn.
                    switch(long(locConn.size())){
                        case 3: //triangles
                            m_intPatch->addConnectedCell(locConn, bitpit::ElementType::TRIANGLE, PID, cOffset, -1);
                            if(textureOn){
                                m_intTexture->addConnectedCell(txtConn, bitpit::ElementType::TRIANGLE, PID, cOffset, -1);
                            }
                        break;
                        case 4: //quads
                            m_intPatch->addConnectedCell(locConn, bitpit::ElementType::QUAD, PID, cOffset, -1);
                            if(textureOn){
                                m_intTexture->addConnectedCell(txtConn, bitpit::ElementType::QUAD, PID, cOffset, -1);
                            }
                        break;
                        default: //polygons
                            {
                                std::vector<long> tt(1,locConn.size());
                                tt.insert(tt.end(), locConn.begin(), locConn.end());
                                m_intPatch->addConnectedCell(tt, bitpit::ElementType::POLYGON, PID, cOffset, -1);
                                if(textureOn){
                                    tt.resize(1);
                                    tt.insert(tt.end(), txtConn.begin(), txtConn.end());
                                    m_intTexture->addConnectedCell(tt, bitpit::ElementType::POLYGON, PID, cOffset, -1);
                                }
                            }
                        break;
                    }
                    //adjusting data;
                    if(activeMaterial != "NONE"){
                        m_intData->materials.insert(cOffset, activeMaterial);
                    }
                    if(activeSmooth != 0){
                        m_intData->smoothids.insert(cOffset, activeSmooth);
                    }
                    ++cOffset;
                }
                break;
            case 4: //usemtl material name
                activeMaterial = " ";
                ss >> activeMaterial;
                if(activeMaterial.empty()) activeMaterial="NONE";
                break;
            case 5: //g cellgroup string labels
                break;
            case 6: //s smoothgroup id
                stringsmooth = "";
                ss >>stringsmooth;
                stringsmooth = bitpit::utils::string::trim(stringsmooth);
                if(stringsmooth.empty() || stringsmooth == "off"){
                    activeSmooth = 0;
                }else{
                    activeSmooth = std::stol(stringsmooth);
                }
                break;
            case 7: //vertex in freeform st vp
            case 8: //line conn l
            case 9: //p conn l
            default: //unsupported flag
                *(m_log)<<"WARNING "<<m_name<<" : unsupported flag declaration "<<key<<" while reading obj file. Ignoring..."<<std::endl;
                break;
        }
    }

    return;
}

/*!
    Write data of an object to OBJ output stream. Write, only:
    - v: vertices coordinates
    - vt: texture coordinates (if any)
    - vn: normals referred to vertices
    - f : facet connectivity
    - usemtl : all the cells after the keyword will be assigned to the material specified by mtllib file (if any)
    - g: all the cells after the keyword will be assigned to cell group specified by g (if any)
    - s: all the cells after the keyword will be assigned to smooth group specified by s (if any)

    Data are taken from m_geometry (for main mesh), m_extTexture(for texture) and
    m_extData(for attached WavefrontObjData datatypes)
    \param[in] out writin stream
    \param[in] vertList list of mesh vertex ids belonging to current object
    \param[in] txtVertList list of texture vertex ids belonging to current object. If list is empty texture is ignored.
    \param[in] cellList list of mesh(and texture) cell ids belonging to current object
    \param[in, out] vOffset i/o current id of the next-to-be inserted mesh vertex.
    \param[in, out] vTxtOffset i/o current id of the next-to-be inserted texture vertex.
    \param[in, out] cOffset i/o current id of the next-to-be inserted mesh cell (same for mesh and texture)
*/
void IOWavefrontOBJ::writeObjectData(std::ofstream & out, const livector1D &vertList,
                                     const livector1D & txtVertList, const livector1D & cellList,
                                    long &vOffset, long &vTxtOffset, long &cOffset)
{
    // create map for current mesh vertex id, and order of insertion counter
    bool shutTexture = txtVertList.empty();
    std::unordered_map<long, long> mapVID_index, mapVTxtID_index;
    darray3E temp;

    //calculate vertex normals
    MimmoPiercedVector<std::array<double,3> > normals;
    normals.reserve(vertList.size());
    bitpit::SurfaceKernel * skpoly = static_cast<bitpit::SurfaceKernel *>(m_geometry->getPatch());
    for(long idC : cellList){
        bitpit::Cell & cell = skpoly->getCell(idC);
        bitpit::ConstProxyVector<long> vIds = cell.getVertexIds();
        int i=0;
        for(long vid: vIds){
            if(normals.exists(vid)) continue;
            temp = skpoly->evalVertexNormal(idC, i);
            normals.insert(vid, temp);
            ++i;
        }
    }

    //write vertices and create the insertion map mapVID_index
    int count(0);
    for(long id: vertList){
        temp = m_geometry->getVertexCoords(id);
        out<<"v "<<std::fixed<<std::setprecision(6)<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<std::endl;
        mapVID_index[id] = vOffset+count;
        ++count;
    }
    vOffset += count;

    if(!shutTexture){
        //texture vertices with its own insertion map
        count = 0;
        for(long id: txtVertList){
            temp = m_extTexture->getVertexCoords(id);
            out<<"vt "<<std::fixed<<std::setprecision(6)<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<std::endl;
            mapVTxtID_index[id] = vTxtOffset+count;
            ++count;
        }
        vTxtOffset += count;
    }

    // vertex normals
    for(long id: vertList){
        out<<"vn "<<std::fixed<<std::setprecision(6)<<normals[id][0]<<" "<<normals[id][1]<<" "<<normals[id][2]<<std::endl;
    }
    normals.clear();

    //connectivity time. Use insertion maps and regroup by materials in m_extData.
    //pay attention also to smoothing ids.
    std::unordered_map<std::string, std::vector<long>> matgCells = regroupCellsByMaterials(cellList);
    MimmoPiercedVector<long> * sid;
    bool deleteSid;

    if(m_extData){
        sid = &(m_extData->smoothids);
        deleteSid = false;
    }else{
        sid = new MimmoPiercedVector<long>(m_geometry, MPVLocation::CELL);
        deleteSid = true;
    }
    sid->completeMissingData(0);

    std::cout<<" material g cells "<<matgCells.size()<<std::endl;

    for(auto & matEntry : matgCells){
        // write material
        if(matEntry.first != "NONE"){
            out<<"usemtl "<<matEntry.first<<std::endl;
        }

        long activeSID = 0;
        std::cout<<matEntry.second.size()<<std::endl;
        for(long idCell : matEntry.second){

            // write smoothing ids
            if(sid->at(idCell) != activeSID){
                if(activeSID == 0){
                    out<<"s off"<<std::endl;
                }else{
                    out<<"s "<<std::to_string(activeSID)<<std::endl;
                }
                activeSID = sid->at(idCell);
            }

            // prepare to write facets
            bitpit::Cell & meshCell = m_geometry->getCells()[idCell];
            std::size_t countV = meshCell.getVertexCount();

            std::vector<std::string> mconn(countV, ""), tconn(countV, "");

            bitpit::ConstProxyVector<long> mvIds= meshCell.getVertexIds();
            for(int i=0; i<countV; ++i ){
                mconn[i] = std::to_string(mapVID_index[mvIds[i]]);
            }
            if(!shutTexture){
                bitpit::ConstProxyVector<long> tvIds= m_extTexture->getCells()[idCell].getVertexIds();
                for(int i=0; i<countV; ++i ){
                    tconn[i] = std::to_string(mapVTxtID_index[tvIds[i]]);
                }
            }

            //write facet line v1/vt1/vn1 v2/vt2,vn2 ...;
            out<<"f ";
            for(int i=0; i<countV; ++i){
                out<<mconn[i]<<"/"<<tconn[i]<<"/"<<mconn[i]<<" ";
            }
            out<<std::endl;
        }// ending loop on cells for single material
    } // ending loop on materials;

    if(deleteSid){
        delete sid;
    }
}

/*!
    Regroup a list of mesh cells by materials name attached to them.
    \param[in] cellList list of cell ids of the polygonal mesh
    \return cell list reordered so that same material cells are clustered together.
*/
std::unordered_map<std::string, std::vector<long>>
IOWavefrontOBJ::regroupCellsByMaterials(const livector1D & cellList){

    std::unordered_map<std::string, std::vector<long>> map;
    if(m_extData == nullptr){
         map["NONE"] = cellList;
         return map;
    }


    for(long id: cellList){
        if(m_extData->materials.exists(id)){
            map[m_extData->materials[id]].push_back(id);
        }else{
            map["NONE"].push_back(id);
        }
    }

    return map;
}

/*!
    convert a key entry useful for reading stage in a int
*/
int IOWavefrontOBJ::convertKeyEntryToInt(const std::string & key){

    int res = -1;
    if(key == "v")  res=0;
    if(key == "vt")  res=1;
    if(key == "vn")  res=2;
    if(key == "f")  res=3;
    if(key == "usemtl")  res=4;
    if(key == "g")  res=5;
    if(key == "s")  res=6;
    if(key == "vp")  res=7;
    if(key == "l")  res=8;
    if(key == "p")  res=9;

    return res;
}




}
