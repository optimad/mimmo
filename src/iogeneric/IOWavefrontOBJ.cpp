/*
    Get a License Header here- Temporary POC for Click-Ins
*/
#include "IOWavefrontOBJ.hpp"
#include "bitpit_common.hpp"
#include "customOperators.hpp"
#include <algorithm>

namespace mimmo{


/*!
 Constructor
*/
WavefrontObjData::WavefrontObjData(){
    materialfile= "UnknownFile.mtl";
    materials.setName("Material");
    smoothids.setName("SmooothingGroupId");
    materials.setDataLocation(MPVLocation::CELL);
    smoothids.setDataLocation(MPVLocation::CELL);
    vertexTexture.setName("Texture");
    vertexNormal.setName("Normal");
    vertexTexture.setDataLocation(MPVLocation::POINT);
    vertexNormal.setDataLocation(MPVLocation::POINT);
}

/*!
    Swap operator. Swap content with a twin class.
    \param[in] other WavefrontObjData object.
*/
void
WavefrontObjData::swap(WavefrontObjData & x) noexcept{
    materials.swap(x.materials);
    smoothids.swap(x.smoothids);
    vertexTexture.swap(x.vertexTexture);
    vertexNormal.swap(x.vertexNormal);
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

    // dump texture vertex
    location  = static_cast<int>(MPVLocation::POINT);
    name = vertexTexture.getName();
    totSize = vertexTexture.size();
    bitpit::utils::binary::write(out, location);
    bitpit::utils::binary::write(out, name);
    bitpit::utils::binary::write(out, totSize);

    for(auto it=vertexTexture.begin(); it!=vertexTexture.end(); ++it){
        bitpit::utils::binary::write(out, it.getId());
        bitpit::utils::binary::write(out, *it);
    }

    // dump normal vertex
    location  = static_cast<int>(MPVLocation::POINT);
    name = vertexNormal.getName();
    totSize = vertexNormal.size();
    bitpit::utils::binary::write(out, location);
    bitpit::utils::binary::write(out, name);
    bitpit::utils::binary::write(out, totSize);

    for(auto it=vertexNormal.begin(); it!=vertexNormal.end(); ++it){
        bitpit::utils::binary::write(out, it.getId());
        bitpit::utils::binary::write(out, *it);
    }

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

    vertexTexture.clear();
    //restore vertex ntexture
    {
    	bitpit::utils::binary::read(in, location);
    	bitpit::utils::binary::read(in, name);
    	bitpit::utils::binary::read(in, totSize);

    	vertexTexture.setName(name);
    	vertexTexture.setDataLocation(static_cast<MPVLocation>(location));
    	vertexTexture.reserve(totSize);

    	std::array<double,3> data;
    	for(int i=0; i<totSize; ++i){
    		bitpit::utils::binary::read(in, id);
    		bitpit::utils::binary::read(in, data);
    		vertexTexture.insert(id, data);
    	}
    }

    vertexNormal.clear();
    //restore vertex ntexture
    {
    	bitpit::utils::binary::read(in, location);
    	bitpit::utils::binary::read(in, name);
    	bitpit::utils::binary::read(in, totSize);

    	vertexNormal.setName(name);
    	vertexNormal.setDataLocation(static_cast<MPVLocation>(location));
    	vertexNormal.reserve(totSize);

    	std::array<double,3> data;
    	for(int i=0; i<totSize; ++i){
    		bitpit::utils::binary::read(in, id);
    		bitpit::utils::binary::read(in, data);
    		vertexNormal.insert(id, data);
    	}
    }

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
    m_texturetol = 1.0E-8;
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
    m_tol = 1.0E-8;
    m_texturetol = 1.0E-8;

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
    std::swap(m_texturetol, x.m_texturetol);

    std::swap(m_intPatch, x.m_intPatch);
    std::swap(m_intData, x.m_intData);
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
    built = (built && createPortIn<WavefrontObjData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setData, X_WDATA));
    built = (built && createPortIn<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setMaterials, M_STRINGFIELD, false, 2));
    built = (built && createPortIn<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setSmoothIds, M_LONGFIELD, false, 2));
    built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setNormals, M_VECTORFIELD, false, 2));
    built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setTexture, M_VECTORFIELD2, false, 2));
    built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setGeometryDisplacements, M_GDISPLS));
    built = (built && createPortIn<std::string, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setMaterialFile, M_NAME, false, 2));

    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getGeometry, M_GEOM));
    built = (built && createPortOut<WavefrontObjData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getData, X_WDATA));
    built = (built && createPortOut<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getMaterials, M_STRINGFIELD));
    built = (built && createPortOut<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getSmoothIds, M_LONGFIELD));
    built = (built && createPortOut<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getNormals, M_VECTORFIELD));
    built = (built && createPortOut<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getTexture, M_VECTORFIELD2));
    built = (built && createPortOut<std::string, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getMaterialFile, M_NAME));
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
    Return string with filename with materials lib.
    \return texture data
*/
std::string   IOWavefrontOBJ::getMaterialFile(){
    if(m_extData != nullptr) return (m_extData->materialfile);
    else if(m_intData != nullptr) return (m_intData->materialfile);
    else return nullptr;
};

/*!
    Return pointer to WavefrontObjData texture. Texture is treated as
    an external surface geometry who shares the same cell ids then polygonal geoemetry
    which refers. The cell type is identical between mesh.
    \return texture data
*/
MimmoPiercedVector<std::array<double, 3>>*   IOWavefrontOBJ::getTexture(){
    if(m_extData != nullptr) return &(m_extData->vertexTexture);
    else if(m_intData != nullptr) return &(m_intData->vertexTexture);
    else return nullptr;
};

/*!
    Return pointer to WavefrontObjData (materials, group, smoothids) attached to polygonal mesh.
    Data are referred always to cell.
    \return data attached to the mesh
*/
WavefrontObjData*   IOWavefrontOBJ::getData(){
	if(m_extData != nullptr) return m_extData;
	else if(m_intData != nullptr) return m_intData.get();
	else return nullptr;

};

/*!
    Return pointer to WavefrontObjData materials attached to polygonal mesh.
    Data are referred always to cell.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<std::string>*   IOWavefrontOBJ::getMaterials(){
	if(m_extData != nullptr) return &(m_extData->materials);
	else if(m_intData != nullptr) return &(m_intData->materials);
	else return nullptr;
};

/*!
    Return pointer to WavefrontObjData normals attached to polygonal mesh.
    Data are referred always to vertex.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<std::array<double, 3>>*   IOWavefrontOBJ::getNormals(){
	if(m_extData != nullptr) return &(m_extData->vertexNormal);
	else if(m_intData != nullptr) return &(m_intData->vertexNormal);
	else return nullptr;
};


/*!
    Return pointer to WavefrontObjData smoothids attached to polygonal mesh.
    Data are referred always to cell.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<long>*   IOWavefrontOBJ::getSmoothIds(){
	if(m_extData != nullptr) return &(m_extData->smoothids);
	else if(m_intData != nullptr) return &(m_intData->smoothids);
	else return nullptr;
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
    Set the geometry meant to be written.  Does nothing in read/restore mode.
    \param[in] geo MimmoObject surface mesh of type 1
*/
void    IOWavefrontOBJ::setGeometry(MimmoObject * geo){
    if(!geo || whichModeInt()<2) return;
    if(geo->getType() != 1) return;
    m_geometry = geo;
}

/*!
    Set the material file linked to the object to be written.
    \param[in] materialfile Name of the material file related to obj file
 */
void    IOWavefrontOBJ::setMaterialFile(std::string materialfile){
    if(whichModeInt()<2) return;
    // I can test extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData != nullptr){
        *(m_log)<<m_name<<" : external data structure already linked before set texture calling. Use linked data structure."<<std::endl;
        return;
    }
    // Initialize internal data structure if not yet initialized
    if (m_intData == nullptr)
    	m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
    m_intData->materialfile = materialfile;
}

/*!
    Set texture mesh related to surface polygonal mesh. They should share
    the same cell-ids.
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] texture mesh.
*/
void    IOWavefrontOBJ::setTexture(MimmoPiercedVector<std::array<double, 3>>* data){
    if(!data || whichModeInt()<2) return;
    // I can test extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData != nullptr){
        *(m_log)<<m_name<<" : external data structure already linked before set texture calling. Use linked data structure."<<std::endl;
        return;
    }
    // Initialize internal data structure if not yet initialized
    if (m_intData == nullptr)
    	m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
    m_intData->vertexTexture = *data;
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
    // I can test extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData != nullptr){
        *(m_log)<<m_name<<" : external data structure already linked before set texture calling. Use linked data structure."<<std::endl;
        return;
    }
    // Initialize internal data structure if not yet initialized
    if (m_intData == nullptr)
    	m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
    m_intData->materials = *data;
};

/*!
    Set normals data attached to surface polygonal mesh, i.e the materials associated to vertices
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] data.
*/
void    IOWavefrontOBJ::setNormals(MimmoPiercedVector<std::array<double, 3>>* data){
    if(!data || whichModeInt()<2) return;
    // I can test extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData != nullptr){
        *(m_log)<<m_name<<" : external data structure already linked before set texture calling. Use linked data structure."<<std::endl;
        return;
    }
    // Initialize internal data structure if not yet initialized
    if (m_intData == nullptr)
    	m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
    m_intData->vertexNormal = *data;
};
/*!
    Set materials data attached to surface polygonal mesh, i.e the smoothing group ids associated to cells
    Meaningful only in write mode. Does nothing in read mode.
    \param[in] data.
*/
void    IOWavefrontOBJ::setSmoothIds(MimmoPiercedVector<long>* data){
    if(!data || whichModeInt()<2) return;
    // I can test extData pointer during port communications in an execution chain
    // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
    if (m_extData != nullptr){
        *(m_log)<<m_name<<" : external data structure already linked before set texture calling. Use linked data structure."<<std::endl;
        return;
    }
    // Initialize internal data structure if not yet initialized
    if (m_intData == nullptr)
    	m_intData = std::unique_ptr<WavefrontObjData>(new WavefrontObjData());
    m_intData->smoothids = *data;
};

/*!
    Set the vertex displacements of the geometry. If filled the normals of the moved vertices are recomputed during write.
    \param[in] geometry displacements.
*/
void    IOWavefrontOBJ::setGeometryDisplacements(MimmoPiercedVector<std::array<double, 3>>* data){
    if(!data)  return;
	m_displacements = *data;
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
    set the geometric tolerance for duplicated vertices collapsing
    \param[in] tolerance
*/
void    IOWavefrontOBJ::setGeometryTolerance(double tolerance){
    m_tol= std::max(std::numeric_limits<double>::min(), tolerance);
}

/*!
    set the texture tolerance for duplicated vertices collapsing
    \param[in] tolerance
*/
void    IOWavefrontOBJ::setTextureTolerance(double tolerance){
    m_texturetol= std::max(std::numeric_limits<double>::min(), tolerance);
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

    if(slotXML.hasOption("GeometryTolerance")){
        input = slotXML.get("GeometryTolerance");
        double value = 1.0E-8;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setGeometryTolerance(value);
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
    slotXML.set("GeometryTolerance", std::to_string(m_tol));
    slotXML.set("TextureTolerance", std::to_string(m_texturetol));

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

        // Update member counters
        m_totalVertexCount = nVertTot;
        m_totalCellCount = nCellTot;

        //reserve geoemtry vertices and cells
        m_intPatch->getVertices().reserve(nVertTot);
        m_intPatch->getCells().reserve(nCellTot);

        //reserve materials and smoothids cells
        m_intData->materials.reserve(nCellTot);
        m_intData->smoothids.reserve(nCellTot);

        //compile safely data options
        m_intData->materials.setGeometry(m_intPatch.get());
        m_intData->smoothids.setGeometry(m_intPatch.get());
        m_intData->vertexTexture.setGeometry(m_intPatch.get());
        m_intData->vertexNormal.setGeometry(m_intPatch.get());
        m_intData->materials.setDataLocation(MPVLocation::CELL);
        m_intData->smoothids.setDataLocation(MPVLocation::CELL);
        m_intData->vertexTexture.setDataLocation(MPVLocation::POINT);
        m_intData->vertexNormal.setDataLocation(MPVLocation::POINT);
        m_intData->materials.setName("Material");
        m_intData->smoothids.setName("SmooothingGroupId");
        m_intData->vertexTexture.setName("Texture");
        m_intData->vertexNormal.setName("Normal");


        long pidObject(0);
        long vOffset(1), vnOffset(1), vTxtOffset(1), cOffset(1);
        //read file object by object
        for(std::streampos & pos : objectPositions){
//            (*m_log)<<m_name<<" : reading object "<<objectNames[pidObject]<<"...";
            readObjectData(in, pos, pidObject,vOffset, vnOffset, vTxtOffset, cOffset);
            m_intPatch->setPIDName(pidObject, objectNames[pidObject]);
            ++pidObject;
//            (*m_log)<<"done"<<std::endl;

        }

        in.close();

    }else{
        *(m_log)<<m_name<<" : impossible to read from obj file "<<filename<<std::endl;
        throw std::runtime_error("IOWavefrontOBJ::read(), impossible reading obj file");
    }



    //shrink to fit all the reserve stuffs;
    m_intPatch->getVertices().shrinkToFit();
    m_intPatch->getCells().shrinkToFit();

    m_intData->materials.shrinkToFit();
    m_intData->smoothids.shrinkToFit();
    m_intData->vertexTexture.shrinkToFit();
    m_intData->vertexNormal.shrinkToFit();

    // sync the lists of intData
    m_intData->syncListsOnData();

    //collapse duplicated vertices if any
    m_intPatch->getPatch()->setTol(m_tol);
    m_intPatch->getPatch()->deleteCoincidentVertices();
    m_intPatch->getPatch()->deleteOrphanVertices();
    //resyncPID
    m_intPatch->resyncPID();

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

    // Use internal or external data object
    WavefrontObjData* objData;
    if (m_extData != nullptr)
    	objData = m_extData;
    else if (m_intData != nullptr)
    	objData = m_intData.get();
    else{
        (*m_log)<<"WARNING: no data in "<<m_name<<". Nothing to write on file "<<filename<<std::endl;
        return;
    }

    std::ofstream out(filename);
    if(out.is_open()){

        out<<"# Optimad's mimmo : "<<m_name<<" OBJ file"<<std::endl;
        out<<"# www.github.com/optimad/mimmo"<<std::endl;
        std::string materialfile = objData->materialfile;
        out<<"mtllib "<< materialfile<<std::endl;

        auto mapParts = getSubParts();
        long vOffset(1), vnOffset(1), vTxtOffset(1), cOffset(1);
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
            if (!objData->vertexTexture.empty())
            	txtVertList = objData->vertexTexture.getIds();

            livector1D vnVertList;
            if (!objData->vertexNormal.empty()){
            	vnVertList = objData->vertexNormal.getIds();
            	if (!m_displacements.empty()){
            		computeMovedNormals();
            	}
            }

//            (*m_log)<<m_name<<" : writing object "<<entryPart.second<<"...";
            writeObjectData(objData, out,vertList,txtVertList, vnVertList, cellList, vOffset, vnOffset, vTxtOffset, cOffset);
//            (*m_log)<<"done"<<std::endl;

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
    bitpit::utils::binary::read(in, m_texturetol);

    bool geoMark, txtMark, dataMark;

    bitpit::utils::binary::read(in, geoMark);
    if(geoMark){
        m_intPatch->restore(in);
    }else{
        m_intPatch = nullptr;
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
    bitpit::utils::binary::write(out, m_texturetol);

    bool geoMark = (m_geometry != nullptr);
    bitpit::utils::binary::write(out, geoMark);
    if(geoMark) m_geometry->dump(out);

    // Use internal or external data object
    WavefrontObjData* objData = nullptr;
    if (m_extData != nullptr)
    	objData = m_extData;
    else if (m_intData != nullptr)
    	objData = m_intData.get();

    bool dataMark = (objData != nullptr);
    bitpit::utils::binary::write(out, dataMark);
    if(dataMark) objData->dump(out);

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
    \param[in, out] vnOffset i/o current id of the next-to-be inserted mesh vertex normal.
    \param[in, out] vTxtOffset i/o current id of the next-to-be inserted texture vertex.
    \param[in, out] cOffset i/o current id of the next-to-be inserted mesh cell (same for mesh and texture)
*/
void IOWavefrontOBJ::readObjectData(std::ifstream & in, const std::streampos &begObjectStream, const long &PID,
                     long &vOffset, long &vnOffset, long &vTxtOffset, long &cOffset)
{
    if (!in.good()){
        in.clear();
    }
    in.seekg(begObjectStream);
    std::string line, key(" ");

    std::string activeMaterial = "NONE", stringsmooth;
    long activeSmooth = 0;
    std::vector<long> locConn, txtConn, vnormConn;

    MimmoPiercedVector<std::array<double,3>> texture; /**< temporary vertex texture of the object */
    std::unordered_map<long, std::vector<long>> temporaryCellTextureConnectivity; /**< temporary map with texture connectivity for each geometry cell Id */

    MimmoPiercedVector<std::array<double,3>> normals; /**< temporary vertex normals of the object */
    std::unordered_map<long, std::vector<long>> temporaryCellNormalConnectivity; /**< temporary map with normal connectivity for each geometry cell Id */

	texture.reserve(m_totalVertexCount);
	normals.reserve(m_totalVertexCount);

    // TODO Make facet definition as a class member and use it during read/write
    // Facet definition label.
    // -1 : not yet defined
    //  0 : only vertex 				v
    //  1 : vertex and texture			v/vt
    //  2 : vertex and normal			v//vn
    //  3 : vertex, texture and normal	v/vt/vn
    int facetdefinition = -1;

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
        {
        	temp.fill(0.0);
        	ss>>temp[0]>>temp[1]>>temp[2];
        	m_intPatch->addVertex(temp, vOffset);
        	++vOffset;
        }
        break;
        case 1: //texture vt
        {
        	// Reserve if not already done
        	if (facetdefinition == -1){
        	}

        	temp.fill(0.0);
        	int c = 0;
        	while(ss.good() && c<3){
        		ss>>temp[c];
        		++c;
        	}
        	texture.insert(vTxtOffset, temp);
        	++vTxtOffset;
        }
        break;
        case 2: //normals vn,
        {
        	temp.fill(0.0);
        	ss>>temp[0]>>temp[1]>>temp[2];
        	normals.insert(vnOffset, temp);
        	++vnOffset;
        }
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
        	vnormConn.resize(tempstring.size());

        	int counter(0);
        	bool textureOn = false;
        	for(std::string & str: tempstring){

        		// Define the type of facet definition if not already defined
        		if (facetdefinition == -1){

        			int doubleSlashOccurence = countSubstring(str, "//");
        			if (doubleSlashOccurence == 1){
        				facetdefinition = 2;
        			}
        			else{
        				int singleSlashOccurence = countSubstring(str, "/");
        				if (singleSlashOccurence == 2){
        					facetdefinition = 3;
        				}
        				else if (singleSlashOccurence == 1){
        					facetdefinition = 1;
        				}
        				else{
        					facetdefinition = 0;
        				}
        			}

        			if (facetdefinition != 0){
        				if (facetdefinition != 2){
        					// Reserve texture cells and vertices
        					// Check if already reserved
        					if (m_intData->vertexTexture.capacity() != m_totalVertexCount){
        						m_intData->vertexTexture.reserve(m_totalVertexCount);
        						m_intData->vertexTexture.reserve(m_totalCellCount);
        					}
        				}
        				if (facetdefinition != 1){
        					// Reserve normal vertices
        					m_intData->vertexNormal.reserve(m_totalVertexCount);
        				}
        			}

        		} // end facet definition

        		std::replace(str.begin(), str.end(), '/', ' ');
        		std::stringstream ss(str);
        		ss>>locConn[counter];

        		//handling negative vertex indexing rule of obj file
        		if(locConn[counter] < 0) locConn[counter] += vOffset;

        		// If facet defined by only vertex skip texture and normal
        		if (facetdefinition != 0){

        			// If facet defined by only vertex and normal skip texture
        			if (facetdefinition != 2){
        				if(ss.good()){
        					ss>>txtConn[counter];
        					if(txtConn[counter] < 0) txtConn[counter] += vTxtOffset;
        				}
        			}

        			// If facet defined by only vertex and texture skip normal
        			if (facetdefinition != 1){
        				if(ss.good()){
        					ss>>vnormConn[counter];
        					if(vnormConn[counter] < 0) vnormConn[counter] += vnOffset;
        				}
        			}

        		} // end if only vertex facet definition

        		//increment counter;
        		++counter;
        	}

        	//adding cell desuming cell type from locConn.
        	switch(long(locConn.size())){
        	case 3: //triangles
        		m_intPatch->addConnectedCell(locConn, bitpit::ElementType::TRIANGLE, PID, cOffset, -1);
        		break;
        	case 4: //quads
        		m_intPatch->addConnectedCell(locConn, bitpit::ElementType::QUAD, PID, cOffset, -1);
        		break;
        	default: //polygons
        	{
        		std::vector<long> tt(1,locConn.size());
        		tt.insert(tt.end(), locConn.begin(), locConn.end());
        		m_intPatch->addConnectedCell(tt, bitpit::ElementType::POLYGON, PID, cOffset, -1);
        	}
        	break;
        	}

        	// Insert item in texture connectivity for cells
        	temporaryCellTextureConnectivity[cOffset] = txtConn;

        	// Insert item in normal connectivity for cells
        	temporaryCellNormalConnectivity[cOffset] = vnormConn;

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

	//Fill vertex normals and texture collocated on vertices of the geometry
	for (bitpit::Cell & cell : m_intPatch->getCells()){

		{
			//Recover temporary cell texture connectivity
			long cellId = cell.getId();
			std::vector<long> & textureConnectivity = temporaryCellTextureConnectivity[cellId];

			// Insert vertex normals
			const bitpit::ConstProxyVector<long> vertices = cell.getVertexIds();
			for (std::size_t i = 0; i < cell.getVertexCount(); i++){
				long vertexId = vertices[i];
				if (!m_intData->vertexTexture.exists(vertexId))
					m_intData->vertexTexture.insert(vertexId, texture[textureConnectivity[i]]);
			}
		}

		{
			//Recover temporary cell normal connectivity
			long cellId = cell.getId();
			std::vector<long> & normalConnectivity = temporaryCellNormalConnectivity[cellId];

			// Insert vertex normals
			const bitpit::ConstProxyVector<long> vertices = cell.getVertexIds();
			for (std::size_t i = 0; i < cell.getVertexCount(); i++){
				long vertexId = vertices[i];
				if (!m_intData->vertexNormal.exists(vertexId))
					m_intData->vertexNormal.insert(vertexId, normals[normalConnectivity[i]]);
			}
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
    \param[in] objData source object data structure
    \param[in] out writing stream
    \param[in] vertList list of mesh vertex ids belonging to current object
    \param[in] txtVertList list of texture vertex ids belonging to current object. If list is empty texture is ignored.
    \param[in] cellList list of mesh(and texture) cell ids belonging to current object
    \param[in, out] vOffset i/o current id of the next-to-be inserted mesh vertex.
    \param[in, out] vnOffset i/o current id of the next-to-be inserted mesh vertex normal.
    \param[in, out] vTxtOffset i/o current id of the next-to-be inserted texture vertex.
    \param[in, out] cOffset i/o current id of the next-to-be inserted mesh cell (same for mesh and texture)
*/
void IOWavefrontOBJ::writeObjectData(WavefrontObjData* objData, std::ofstream & out, const livector1D &vertList,
                                     const livector1D & txtVertList, const livector1D & vnVertList, const livector1D & cellList,
                                    long &vOffset, long &vnOffset, long &vTxtOffset, long &cOffset)
{

    // create map for current mesh vertex id, and order of insertion counter
    bool shutTexture = txtVertList.empty();
    std::unordered_map<long, long> mapVID_index, mapVTxtID_index, mapVnID_index;
    darray3E temp;

    //write vertices and create the insertion map mapVID_index
    int count(0);
    for(long id: vertList){
        temp = m_geometry->getVertexCoords(id);
        out<<"v "<<std::fixed<<std::setprecision(6)<<temp[0]<<" "<<temp[1]<<" "<<temp[2]<<std::endl;
        mapVID_index[id] = vOffset+count;
        ++count;
    }
    vOffset += count;

    //texture vertices with its own insertion map
	bitpit::KdTree<3,std::array<double,3>,long> kdTreeTexture;
    count = 0;
    for(long idTexture: txtVertList){
    	darray3E & texturePoint = objData->vertexTexture[idTexture];

    	// Check if the texture point is already written
    	// Set a texture tolerance
		int inode = kdTreeTexture.hNeighbor(&texturePoint, m_texturetol, false);

		// If returned node index is < 0 no coincident texture points are found
		if (inode < 0){
			long index = vTxtOffset+count;
			kdTreeTexture.insert(&texturePoint, index);
			out<<"vt "<<std::fixed<<std::setprecision(6)<<texturePoint[0]<<" "<<texturePoint[1]<<" "<<texturePoint[2]<<std::endl;
			mapVTxtID_index[idTexture] = index;
			++count;
		}
		else{
			mapVTxtID_index[idTexture] = kdTreeTexture.nodes[inode].label;
		}
    }
    vTxtOffset += count;

    // vertex normals
    count = 0;
    for(long idNormal: vnVertList){

    	// Access to normal
    	std::array<double,3> normal = objData->vertexNormal[idNormal];

    	// If not normalized normalize it
    	double value = norm2(normal);
    	if (!bitpit::utils::DoubleFloatingEqual()(value, 1.)){
    		normal /= value;
    	}

    	// Write normal
    	out<<"vn "<<std::fixed<<std::setprecision(6)<< normal[0]<<" "<< normal[1]<<" "<< normal[2]<<std::endl;
    	mapVnID_index[idNormal] = vnOffset+count;
    	++count;
    }
    vnOffset += count;

    //connectivity time. Use insertion maps and regroup by materials in m_extData.
    //pay attention also to smoothing ids.
    std::unordered_map<std::string, std::vector<long>> matgCells = regroupCellsByMaterials(objData, cellList);
    MimmoPiercedVector<long>* sid;
    bool deleteSid;

    sid = &(objData->smoothids);
    deleteSid = false;
    sid->completeMissingData(0);

//    std::cout<<" material g cells "<<matgCells.size()<<std::endl;

    for(auto & matEntry : matgCells){
        // write material
        if(matEntry.first != "NONE"){
            out<<"usemtl "<<matEntry.first<<std::endl;
        }

        long activeSID = 0;
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

            std::vector<std::string> mconn(countV, ""), tconn(countV, ""),  vnconn(countV, "");

            bitpit::ConstProxyVector<long> mvIds= meshCell.getVertexIds();
            for(int i=0; i<countV; ++i ){
                mconn[i] = std::to_string(mapVID_index[mvIds[i]]);
            }

            // Fill texture connectivity
            for(int i=0; i<countV; ++i ){
            	tconn[i] = std::to_string(mapVTxtID_index[mvIds[i]]);
            }

            // Fill normals connectivity
            for(int i=0; i<countV; ++i ){
            	vnconn[i] = std::to_string(mapVnID_index[mvIds[i]]);
            }

            //write facet line v1/vt1/vn1 v2/vt2,vn2 ...;
            out<<"f ";
            for(int i=0; i<countV; ++i){
                out<<mconn[i]<<"/"<<tconn[i]<<"/"<<vnconn[i]<<" ";
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
    \param[in] objData source object data structure
    \param[in] cellList list of cell ids of the polygonal mesh
    \return cell list reordered so that same material cells are clustered together.
*/
std::unordered_map<std::string, std::vector<long>>
IOWavefrontOBJ::regroupCellsByMaterials(const WavefrontObjData* objData, const livector1D & cellList){

    std::unordered_map<std::string, std::vector<long>> map;
    for(long id: cellList){
        if(objData->materials.exists(id)){
            map[objData->materials[id]].push_back(id);
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

/*!
 * Compute the vertex normals of the vertices involved by geometry displacements.
 * Use geometry tolerance as threshold on the norm of the displacement vector
 * to define if a vertex is moved.
 */
void IOWavefrontOBJ::computeMovedNormals(){

	if (m_displacements.getDataLocation() != MPVLocation::POINT){
        (*m_log)<<m_name + " : geometry displacements not on point location"<<std::endl;
        return;
	}

	if (m_displacements.getGeometry() != getGeometry()){
        (*m_log)<<m_name + " : geometry linked in displacements different ffrom target geometry "<<std::endl;
        return;
	}

    // Use internal or external data object
    WavefrontObjData* objData;
    if (m_extData != nullptr)
    	objData = m_extData;
    else if (m_intData != nullptr)
    	objData = m_intData.get();
    else{
        (*m_log)<<"WARNING: no data in "<<m_name<<". Nothing to do in normals computing. "<<std::endl;
        return;
    }

	// Compute vertex normals
	std::map<long, std::array<double,3>> vertexNormals;
	for (const bitpit::Cell & cell : getGeometry()->getCells()){
		long cellId = cell.getId();
		for (std::size_t ivertex=0; ivertex<cell.getVertexCount(); ivertex++){
			long vertexId = cell.getVertexId(ivertex);
			if (!vertexNormals.count(vertexId)){
				vertexNormals[vertexId] = static_cast<bitpit::SurfaceKernel*>(getGeometry()->getPatch())->evalVertexNormal(cellId, ivertex);
			}
		}
	}

	// Define a local tolerance
//	double tolerance = 1.0e-05;
	double tolerance = m_tol;
	for (long vertexId : getGeometry()->getPatch()->getVertices().getIds()){
		if (norm2(m_displacements[vertexId]) > tolerance){
			std::array<double,3> vertexNormal = vertexNormals[vertexId];
			objData->vertexNormal[vertexId] = vertexNormal;
		}
	}
}


}
