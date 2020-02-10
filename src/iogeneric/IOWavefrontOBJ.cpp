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
#include "IOWavefrontOBJ.hpp"
#include "bitpit_common.hpp"
#include "customOperators.hpp"
#include <algorithm>

namespace mimmo{

////////////////////////////////////////////////////////////////////////////////
//////////////        WAVEFRONTOBJDATA IMPLEMENTATION               ////////
////////////////////////////////////////////////////////////////////////////////

/*!
 Constructor
*/
WavefrontOBJData::WavefrontOBJData(){
    materialfile= "UnknownFile.mtl";
    materials.setName("Material");
    smoothids.setName("SmooothingGroupId");
    materials.setDataLocation(MPVLocation::CELL);
    smoothids.setDataLocation(MPVLocation::CELL);
}

/*!
    Swap operator. Swap content with a twin class.
    \param[in] other WavefrontOBJData OBJect.
*/
void
WavefrontOBJData::swap(WavefrontOBJData & x) noexcept{
    materials.swap(x.materials);
    smoothids.swap(x.smoothids);
    std::swap(materialsList, x.materialsList);
    std::swap(smoothidsList, x.smoothidsList);
    std::swap(inv_materialsList, x.inv_materialsList);
    std::swap(inv_smoothidsList, x.inv_smoothidsList);
    std::swap(materialfile, x.materialfile);
    std::swap(refGeometry, x.refGeometry);
    std::swap(textures, x.textures);
    std::swap(normals, x.normals);
}

/*!
    Synchronize the lists with the inner data currently present in the class.
    Following convention is used:
     - id = 0 of materials is string '' (ie no material is expressed for the cell)
     - id = 0 of smooth ids is 'off' key for smoothing (no smoothing is applied to the cell)
*/
void
WavefrontOBJData::syncListsOnData(){
    //materials
    materialsList.clear();
    long id = 1;
    std::size_t mapsize(0);
    for(auto it = materials.begin(); it!= materials.end(); ++it){
        if (it->empty()){
            materialsList.insert({{*it, long(0)}});
            mapsize = materialsList.size();
        }else{
            materialsList.insert({{*it, id}});
            id += long(mapsize != materialsList.size());
            mapsize = materialsList.size();
        }
    }

    //smoothids
    smoothidsList.clear();
    for(auto it = smoothids.begin(); it!= smoothids.end(); ++it){
        if(*it == 0){
            smoothidsList.insert({{"off", *it}});
        }else{
            smoothidsList.insert({{std::to_string(*it), *it}});
        }
    }

    //compute the inverse maps
    inv_materialsList.clear();
    for(auto & entry: materialsList){
        inv_materialsList[entry.second] = entry.first;
    }
    inv_smoothidsList.clear();
    for(auto & entry: smoothidsList){
        inv_smoothidsList[entry.second] = entry.first;
    }
}
/*!
    autocomplete cell field information (material, smoothids, and cellgroups)) with
    default values ('', 0, '' respectively).
*/
void
WavefrontOBJData::autoCompleteCellFields(){
    smoothids.completeMissingData(0);
    materials.completeMissingData(std::string());
}

/*!
    Dump class contents to an output stream in binary format
    \param[in] out binary output stream
*/
void
WavefrontOBJData::dump(std::ostream & out){
    int location  = static_cast<int>(MPVLocation::CELL);
    std::string name;
    std::size_t totSize, listSize;

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
    // ref geometry cannot BE DUMPED

    // dump textures
    bool textureOn = textures.get() != nullptr;
    bitpit::utils::binary::write(out, textureOn);
    if(textureOn)   textures->dump(out);
    // dump normals
    bool normalOn = normals.get() != nullptr;
    bitpit::utils::binary::write(out, normalOn);
    if(normalOn)   normals->dump(out);

}
/*!
    Restore class contents from an input stream in binary format
    \param[in] in binary input stream
*/
void
WavefrontOBJData::restore(std::istream & in){

    int location, status;
    std::string name;
    std::size_t totSize, listSize;
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
    // refGeometry not restorable. You need to connect in another moment.
    refGeometry = nullptr;

    //restore textures
    bool textureOn;
    bitpit::utils::binary::read(in, textureOn);
    if(textureOn){
        textures = std::unique_ptr<MimmoObject>(new MimmoObject(1));
        textures->restore(in);
    }
    //restore normals
    bool normalOn;
    bitpit::utils::binary::read(in, normalOn);
    if(normalOn){
        normals = std::unique_ptr<MimmoObject>(new MimmoObject(1));
        normals->restore(in);
    }

    // resync lists
    syncListsOnData();
}
////////////////////////////////////////////////////////////////////////////////
//////////////         MANIPULATEWFOBJDATA IMPLEMENTATION               ////////
////////////////////////////////////////////////////////////////////////////////

/*!
    Constructor
*/
ManipulateWFOBJData::ManipulateWFOBJData(){
    m_name = "mimmo.ManipulateWFOBJData";
    m_extData = nullptr;
    m_checkNormalsMag = false;
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ManipulateWFOBJData::ManipulateWFOBJData(const bitpit::Config::Section & rootXML){
    m_name = "mimmo.ManipulateWFOBJData";
    m_extData = nullptr;
    m_checkNormalsMag = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);

    input = bitpit::utils::string::trim(input);

    if(input == "mimmo.ManipulateWFOBJData"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
    Destructor
*/
ManipulateWFOBJData::~ManipulateWFOBJData(){};

/*!
    Swap utility
    \param[in] x other object to swap from
*/
void ManipulateWFOBJData::swap(ManipulateWFOBJData & x) noexcept{
    std::swap(m_extData, x.m_extData);
    m_displacements.swap(x.m_displacements);
    std::swap(m_annotations, x.m_annotations);
    std::swap(m_checkNormalsMag, x.m_checkNormalsMag);

    BaseManipulation::swap(x);
}

/*!
    Building class ports
*/
void ManipulateWFOBJData::buildPorts(){
    bool built = true;
    built = (built && createPortIn<WavefrontOBJData*, ManipulateWFOBJData>(this, &ManipulateWFOBJData::setData, M_WAVEFRONTDATA, true));
    built = (built && createPortIn<MimmoPiercedVector<long>*, ManipulateWFOBJData>(this, &ManipulateWFOBJData::addAnnotation, M_LONGFIELD));
    built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, ManipulateWFOBJData>(this, &ManipulateWFOBJData::setGeometryDisplacements, M_VECTORFIELD));

    built = (built && createPortOut<WavefrontOBJData*, ManipulateWFOBJData>(this, &ManipulateWFOBJData::getData, M_WAVEFRONTDATA));
    m_arePortsBuilt = built;
}

/*!
    Return pointer to WavefrontOBJData (materials, smoothids, textures and normals, etc...)
    linked.
    \return data attached to the Wavefront mesh
*/
WavefrontOBJData*   ManipulateWFOBJData::getData(){
    return m_extData;
};


/*!
    Set the vertex displacements of the geometry. If filled, the vertex normals of the
    object sub-parts are recomputed.
*/
void    ManipulateWFOBJData::setGeometryDisplacements(MimmoPiercedVector<std::array<double, 3>>* data){
    if(!data)  return;
    m_displacements = *data;
};

/*!
    Link the target WavefrontOBJData (materials, smoothids, textures and normals, etc...)
    to be manipulated.
    \param[in] data attached to the Wavefront mesh
*/
void  ManipulateWFOBJData::setData(WavefrontOBJData * data){
    if(!data)   return;
    if(!data->refGeometry){
        *(m_log)<<"Error in "<<m_name<<" : cannot connect WavefrontOBJData with nullptr refGeometry member"<<std::endl;
    }
    m_extData = data;

};

/*!
    Add a field of cell ids as annotation for the Wavefront mesh. Annotated cell ids will
    be redistributed in new PIDs for the refGeometry of the WavefrontOBJData.
    F.E. If you want to annotate 2 mesh cells with the name "marked", and this cells belong
    on two distinct mesh sub-parts/PIDs part1,part2, then 2 new PIDs will be created with the
    name part1_marked and part2_marked.
    BEWARE: Annotation name must be specified as name of the MimmoPiercedVector structure. Multiple
    annotations must be added with a unique name, otherwise they will be discarded.
*/
void    ManipulateWFOBJData::addAnnotation(MimmoPiercedVector<long>* data){
    if(!data)  return;
    std::string name = data->getName();
    for(MimmoPiercedVector<long> & ann: m_annotations){
        if(name == ann.getName())  return;
    }
    m_annotations.push_back(*data);
};

/*!
    If any vertex normals data are present in linked Wavefront data object, check
    their magnitudes to be always unitary, and if not fix them.
    \param[in] flag true/false to activate check & fix of normals magnitudes.
*/

void    ManipulateWFOBJData::setCheckNormalsMagnitude(bool flag){
    m_checkNormalsMag = flag;
};

/*!
    Absorb class parameters from a section xml
    \param[in] slotXML xml section
    \param[in] name unused string
*/
void    ManipulateWFOBJData::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name ){
    BITPIT_UNUSED(name);
    std::string input;

    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("CheckNormalsMagnitude")){
        input = slotXML.get("CheckNormalsMagnitude");
        bool value = false;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setCheckNormalsMagnitude(value);
    };

}

/*!
    Flush class parameters to a section xml
    \param[in] slotXML xml section to be filled
    \param[in] name unused string
*/
void    ManipulateWFOBJData::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    slotXML.set("CheckNormalsMagnitude", std::to_string(m_checkNormalsMag));
}

/*!
    Class workflow execution
*/
void   ManipulateWFOBJData::execute(){

    if(!m_extData) {
        *(m_log)<<"WARNING in "<<m_name<<" : no WavefrontOBJData linked to the class. Execution skipped"<<std::endl;
    }

    if(m_extData->normals.get() && m_checkNormalsMag){
        checkNormalsMagnitude();
    }

    computeAnnotations();

    computeMovedNormals();

}

/*!
    Check and fix normals magnitudes. They must be always unitary.
    It assumes m_extData has already a valid normals MimmoObject in it,
    no check is perfomed internally.
*/
void ManipulateWFOBJData::checkNormalsMagnitude(){

    std::array<double,3> temp;
    double norm_temp;
    for(auto it = m_extData->normals->getPatch()->vertexBegin(); it!=m_extData->normals->getPatch()->vertexEnd(); ++it){
        temp = it->getCoords();
        norm_temp = norm2(temp);
        if(!bitpit::utils::DoubleFloatingEqual()(norm_temp, 1.0) && norm_temp > std::numeric_limits<double>::min()){
            temp /= norm_temp;
            m_extData->normals->modifyVertex(temp, it.getId());
        }
    }

}
/*!
    Compute annotations and send it as new pids to the reference geometry.
    New Created PID will be named as the oldPid (if any) + _name of annotation.
*/
void ManipulateWFOBJData::computeAnnotations(){

    MimmoObject * geo  = m_extData->refGeometry;
    livector1D geoCellsIds = geo->getCells().getIds();

    //find starting pid from m_extData->refgeometry
    std::unordered_set<long> & meshPidList = geo->getPIDTypeList();
    std::unordered_map<long, std::string> & meshPidListNames = geo->getPIDTypeListWNames();

    long pid_current = *(std::max_element(meshPidList.begin(), meshPidList.end()));
    //increment by 1
    ++pid_current;

    //check for annotations:
    bitpit::PiercedVector<bitpit::Cell> & geoCells = geo->getCells();

    for(MimmoPiercedVector<long> & data : m_annotations){

        //clean up ids not in the current geo.
        data.squeezeOutExcept(geoCellsIds, false);
        //iterate on contents and set Pids.
        std::unordered_map<long,long> oldPidNewPid;
        long oldPid;
        for(auto it=data.begin(); it!= data.end(); ++it){

            oldPid = geoCells.at(it.getId()).getPID();

            if(oldPidNewPid.count(oldPid) == 0){

                oldPidNewPid.insert({{oldPid, pid_current}});
                geo->setPIDCell(it.getId(), pid_current);
                ++pid_current;

            }else{

                geo->setPIDCell(it.getId(), oldPidNewPid[oldPid]);
            }
        }

        //adjust PID lists informations on geo.
        std::string ann_suffix = std::string("_") + data.getName();
        for(auto & entry: oldPidNewPid){
            std::string pidName = meshPidListNames.at(entry.first) + ann_suffix;
            geo->setPIDName(entry.second, pidName);
        }

    }// end of annotations.


    //resyncronize PID to erase empty cell pid if any
    geo->resyncPID();

    //reverse this new pidding inside textures and normals too, since they must share
    // the same set of cells
    MimmoObject * txt  = m_extData->textures.get();
    MimmoObject * normals = m_extData->normals.get();

    std::unordered_map<long, long> completepids;

    if(txt){
        completepids = geo->getPID();
        txt->setPID(completepids);
    }
    if(normals){
        if(completepids.empty())    completepids = geo->getPID();
        normals->setPID(completepids);
    }

}

/*!
    Computing normals for mesh with a displacement attached.
    DO NOTHING FOR NOW>
*/
void ManipulateWFOBJData::computeMovedNormals(){
    //do nothing
}


////////////////////////////////////////////////////////////////////////////////
//////////////              IOWAVEFRONTOBJ IMPLEMENTATION               ////////
////////////////////////////////////////////////////////////////////////////////

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
    m_tol = std::numeric_limits<double>::min();
    m_cleanDoubleVertices = false;
    m_textureUVMode = false;
    m_extData = nullptr;
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
    m_tol = std::numeric_limits<double>::min();
    m_cleanDoubleVertices = false;
    m_textureUVMode = false;
    m_extData = nullptr;

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
    std::swap(m_cleanDoubleVertices, x.m_cleanDoubleVertices);
    std::swap(m_textureUVMode, x.m_textureUVMode);

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
    bool mandatoryWrite = false;
    switch(m_mode){
        case IOMode::WRITE:
        case IOMode::DUMP:
            mandatoryWrite = true;
            break;
        default: //all other cases, read, restore
            break;
    }

    built = (built && createPortIn<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setGeometry, M_GEOM, mandatoryWrite));
    built = (built && createPortIn<WavefrontOBJData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setData, M_WAVEFRONTDATA));
    // built = (built && createPortIn<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setMaterials, M_STRINGFIELD, false, 2));
    // built = (built && createPortIn<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setSmoothIds, M_LONGFIELD, false, 2));
    // built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::addNormals, M_VECTORFIELD, false, 3));
    // built = (built && createPortIn<MimmoPiercedVector<std::array<double,3>>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::addTexture, M_VECTORFIELD2, false,3));
    // built = (built && createPortIn<std::vector<MimmoPiercedVector<std::array<double,3>>*>, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setAllNormals, M_VECVECTORFIELDS, false, 2));
    // built = (built && createPortIn<std::vector<MimmoPiercedVector<std::array<double,3>>*>, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setAllTextures, M_VECVECTORFIELDS2, false, 2));
    // built = (built && createPortIn<std::string, IOWavefrontOBJ>(this, &IOWavefrontOBJ::setMaterialFile, M_NAME, false, 2));

    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getGeometry, M_GEOM));
    built = (built && createPortOut<WavefrontOBJData*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getData, M_WAVEFRONTDATA));
    built = (built && createPortOut<MimmoPiercedVector<std::string>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getMaterials, M_STRINGFIELD));
    built = (built && createPortOut<MimmoPiercedVector<long>*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getSmoothIds, M_LONGFIELD));
    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getNormals, M_GEOM3));
    built = (built && createPortOut<MimmoObject*, IOWavefrontOBJ>(this, &IOWavefrontOBJ::getTextures, M_GEOM2));
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
    Return pointer to WavefrontOBJData (materials, smoothids, textures and normals, etc...)
    attached to polygonal mesh.
    Data are referred to cell (materials, smoothids ...) or points(texture, normals, ...).
    \return data attached to the mesh
*/
WavefrontOBJData*   IOWavefrontOBJ::getData(){
	if(m_extData != nullptr)      return m_extData;
	else if(m_intData != nullptr) return m_intData.get();
	else return nullptr;

};

/*!
    Return the map with PID and name of the object which the mesh is subdivided in.
    In read/restore mode, this subdivision is extrapolated from file obj, in write/dump mode this refers to
    subdivision hold by the geometry set with setGeometry method.
    \return pid/name mesh internal subdivision.
*/
std::unordered_map<long, std::string> IOWavefrontOBJ::getSubParts(){
    if(!getGeometry())          return std::unordered_map<long, std::string>();

    return getGeometry()->getPIDTypeListWNames();
}

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
    if(m_extData != nullptr)      return (m_extData->materialfile);
    else if(m_intData != nullptr) return (m_intData->materialfile);
    else return nullptr;
};

/*!
    Return pointer to WavefrontOBJData materials attached to polygonal mesh.
    Data are referred always to cell.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<std::string>*   IOWavefrontOBJ::getMaterials(){
	if(m_extData != nullptr)       return &(m_extData->materials);
	else if(m_intData != nullptr)  return &(m_intData->materials);
	else return nullptr;
};


/*!
    Return pointer to WavefrontOBJData smoothids attached to polygonal mesh.
    Data are referred always to cell.
    \return materials data attached to the mesh
*/
MimmoPiercedVector<long>*   IOWavefrontOBJ::getSmoothIds(){
	if(m_extData != nullptr)       return &(m_extData->smoothids);
	else if(m_intData != nullptr)  return &(m_intData->smoothids);
	else return nullptr;
};

/*!
    Return texture object
    \return texture data
*/
MimmoObject *
IOWavefrontOBJ::getTextures()
{
    if(m_extData != nullptr)       return (m_extData->textures.get());
	else if(m_intData != nullptr)  return (m_intData->textures.get());
	else return nullptr;
};


/*!
    Return normals object
    \return normals data
*/
MimmoObject *
IOWavefrontOBJ::getNormals()
{
    if(m_extData != nullptr)       return (m_extData->normals.get());
	else if(m_intData != nullptr)  return (m_intData->normals.get());
	else return nullptr;
};

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
    Set data attached to surface polygonal mesh, i.e:
    - materials associated to cells
    - cell group labels associated to cells
    - smoothing group ids associated to cells
    - texture and vnormals list if any
    - material file path.
    Meaningful only in write mode. Does nothing in read mode.
    Beware any other data inserted manually (setMaterials, setMaterialFile,
    setAlltextures/addTextures, etc...) are erased.
    \param[in] data.
*/
void    IOWavefrontOBJ::setData(WavefrontOBJData* data){
    if(!data || whichModeInt()<2) return;
    m_extData = data;
    m_intData = nullptr;
};

// /*!
//     Set the material file linked to the object to be written.
//     If setData is already used this method will set nothing.
//     Meaningful only in write mode. Does nothing in read mode.
//     \param[in] materialfile Name of the material file related to obj file
//  */
// void    IOWavefrontOBJ::setMaterialFile(std::string materialfile){
//     if(whichModeInt()<2) return;
//     // I can test extData pointer during port communications in an execution chain
//     // because in the ports this method has a lower priority than setData (i.e. this port is executed after)
//     if (m_extData != nullptr){
//         *(m_log)<<m_name<<" : external data structure already linked before setMaterialFile calling. Using previously linked data structure."<<std::endl;
//         return;
//     }
//     // Initialize internal data structure if not yet initialized
//     if (!m_intData.get())
//     	m_intData = std::unique_ptr<WavefrontOBJData>(new WavefrontOBJData());
//     m_intData->materialfile = materialfile;
// }


/*!
    Set the directory where the read file is searched
    or the written file has to be placed.
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
    Set the geometric tolerance for duplicated vertices/vertex normals collapsing
    \param[in] tolerance
*/
void    IOWavefrontOBJ::setGeometryTolerance(double tolerance){
    m_tol= std::max(std::numeric_limits<double>::min(), tolerance);
}


/*!
    For READ mode only, if true force double mesh vertices collapsing inside geometric
    tolerance m_tol (setGeometryTolerance method). This collapsing will not affect
    wavefront obj data attached to the mesh. Do nothing if false.
    \param[in] clean true/false
*/
void    IOWavefrontOBJ::setCleanDoubleMeshVertices(bool clean){
    m_cleanDoubleVertices= clean;
}

/*!
    For WRITE mode only, if any textures is available in wavefront optional data
    force writing textures with only the first two components U and V (ignoring optional W components.).
    \param[in] UVmode true/false to activate/deactivate texture UV mode.
*/
void    IOWavefrontOBJ::setTextureUVMode(bool UVmode){
    m_textureUVMode= UVmode;
}

/*!
    If set to true print a resume file of the data referenced by the mesh.
    \param[in] print true/false
*/
void    IOWavefrontOBJ::printResumeFile(bool print){
    m_resume= print;
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
        double value = std::numeric_limits<double>::min();
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setGeometryTolerance(value);
    };


    if(slotXML.hasOption("CleanDoubleMeshVertices")){
        input = slotXML.get("CleanDoubleMeshVertices");
        bool value = false;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setCleanDoubleMeshVertices(value);
    };

    if(slotXML.hasOption("TextureUVMode")){
        input = slotXML.get("TextureUVMode");
        bool value = false;
        if(!input.empty()){
            input = bitpit::utils::string::trim(input);
            std::stringstream ss(input);
            ss >> value;
        }
        setTextureUVMode(value);
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
    slotXML.set("CleanDoubleMeshVertices", std::to_string(m_cleanDoubleVertices));
    slotXML.set("TextureUVMode", std::to_string(m_textureUVMode));
}


/*!
    Class workflow execution
*/
void   IOWavefrontOBJ::execute(){
    switch(m_mode) {
        case IOMode::READ:
            //instantiation of a brand new MimmoObject to absorb grid
            m_intPatch = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            m_intData = std::unique_ptr<WavefrontOBJData>(new WavefrontOBJData());
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
            m_intData = std::unique_ptr<WavefrontOBJData>(new WavefrontOBJData());
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
    Read mesh and data attached from file obj.
*/
void IOWavefrontOBJ::read(const std::string & filename){

    std::ifstream in(filename, std::ios::binary);
    if(in.is_open()){

        //search and store the material file fullpath.
        m_intData->materialfile = searchMaterialFile(in);
        //search and store the stream position of the sub-objects
        std::vector<std::streampos> objectPositions; //streampos of each new object
        std::vector<std::string> objectNames; // name of the new object
        std::vector<std::array<long,3>> objectVCounters; //number of v, vt, and vn in each object
        std::array<long,3>  totVCounters; //total number of v, vt, and vn
        long nCellTot; // number of total facet f

        // open the file and fill the information required.
        searchObjectPosition(in, objectPositions, objectNames, objectVCounters, totVCounters, nCellTot);

        //reserve geometry vertices and cells
        m_intPatch->getVertices().reserve(totVCounters[0]);
        m_intPatch->getCells().reserve(nCellTot);

        //compile safely data options
        m_intData->materials.reserve(nCellTot);
        m_intData->smoothids.reserve(nCellTot);

        m_intData->materials.setGeometry(m_intPatch.get());
        m_intData->smoothids.setGeometry(m_intPatch.get());
        m_intData->materials.setDataLocation(MPVLocation::CELL);
        m_intData->smoothids.setDataLocation(MPVLocation::CELL);
        m_intData->materials.setName("Material");
        m_intData->smoothids.setName("SmooothingGroupId");

        // prepare textures and normals
        m_intData->textures = nullptr;
        if(totVCounters[1]> 0){
            m_intData->textures = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            m_intData->textures->getVertices().reserve(totVCounters[1]);
            m_intData->textures->getCells().reserve(nCellTot);
        }
        m_intData->normals = nullptr;
        if(totVCounters[2]> 0){
            m_intData->normals = std::unique_ptr<MimmoObject>(new MimmoObject(1));
            m_intData->normals->getVertices().reserve(totVCounters[2]);
            m_intData->normals->getCells().reserve(nCellTot);
        }

        long pidObject(0);
        long vOffset(1), vnOffset(1), vTxtOffset(1), cOffset(1);
        //read file object by object
        for(std::streampos & pos : objectPositions){

            readObjectData(in, pos, pidObject,objectVCounters[pidObject], vOffset, vnOffset, vTxtOffset, cOffset);
            m_intPatch->setPIDName(pidObject, objectNames[pidObject]);
            ++pidObject;
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

    if(m_intData->textures.get()){
        m_intData->textures->getVertices().shrinkToFit();
        m_intData->textures->getCells().shrinkToFit();
    }
    if(m_intData->normals.get()){
        m_intData->normals->getVertices().shrinkToFit();
        m_intData->normals->getCells().shrinkToFit();
    }

    // sync the lists of intData
    m_intData->syncListsOnData();
    m_intData->refGeometry = m_intPatch.get();

    //original mesh cleaning stuffs
    m_intPatch->getPatch()->setTol(m_tol);
    // delete orphan vertices not in the tessellation.
    m_intPatch->getPatch()->deleteOrphanVertices();
    //force cleaning of double vertices part.
    if(m_cleanDoubleVertices){
        m_intPatch->getPatch()->deleteCoincidentVertices();
    }

}

/*!
    Write mesh and data attached to file obj
*/
void IOWavefrontOBJ::write(const std::string & filename){
    if(!m_geometry){
        (*m_log)<<"WARNING: no geometry linked in "<<m_name<<". Nothing to write on file "<<filename<<'\n';
        return;
    }
    if(!m_geometry->areAdjacenciesBuilt())  m_geometry->buildAdjacencies();

    // Use internal or external data object
    WavefrontOBJData* objData = getData();
    if (objData == nullptr){
        (*m_log)<<"WARNING: no OBJ optional data in "<<m_name<<". Nothing to write on file "<<filename<<'\n';
        return;
    }

    if (objData->refGeometry != m_geometry){
        (*m_log)<<"WARNING: in "<<m_name<<". OBJ optional data linked refers to a different geometry w.r.t that linked in the current class. Nothing to write on "<<filename<<'\n';
        return;
    }


    objData->autoCompleteCellFields();
    objData->syncListsOnData();

    std::array<std::vector<long>,3> vertexLists;
    std::vector<long> cellList;
    long cOffset(1);
    std::array<long,3> vOffsets;
    vOffsets.fill(1);
    std::array<std::unordered_map<long,long>,3> vinsertion_maps; // key long id, written id.
    vinsertion_maps[0].reserve(m_geometry->getNVertices());
    if(objData->textures.get())     vinsertion_maps[1].reserve(objData->textures->getNVertices());
    if(objData->normals.get())     vinsertion_maps[2].reserve(objData->normals->getNVertices());

    long activeMaterial = 0;
    long activeSmoothId = 0;

    std::ofstream out(filename, std::ios::binary);
    if(out.is_open()){

        out<<"# Optimad's mimmo : "<<m_name<<" OBJ file"<<'\n';
        out<<"# www.github.com/optimad/mimmo"<<'\n';
        std::string materialfile = objData->materialfile;
        out<<"mtllib "<< materialfile<<'\n';

        std::map<long, livector1D> pidSubdivision = m_geometry->extractPIDSubdivision();
        std::map<long, std::string> mapParts;
        {
            std::unordered_map<long, std::string> mapParts_temp = getSubParts();
            mapParts.insert(mapParts_temp.begin(), mapParts_temp.end());
        }

        long refPid;

        //write object by object
        for(auto & entryPart : mapParts){
            //write header object;
            out<<"o "<<entryPart.second<<'\n';

            refPid = entryPart.first;
            //select list of vertices and cells referring to this pid.
            cellList = pidSubdivision[refPid];
            vertexLists[0] = m_geometry->getVertexFromCellList(cellList);

            if (objData->textures.get()){
                vertexLists[1] = objData->textures->getVertexFromCellList(cellList);
            }else{
                vertexLists[1].clear();
            }
            if (objData->normals.get()){
                vertexLists[2] = objData->normals->getVertexFromCellList(cellList);
            }else{
                vertexLists[2].clear();
            }

            writeObjectData(objData, out, vertexLists, cellList, vOffsets,
                            cOffset, vinsertion_maps, activeMaterial, activeSmoothId);
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
    // m_intPatch, m_intData are supposed to be
    //initialized

    bitpit::utils::binary::read(in, m_tol);

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
        m_intData->refGeometry = m_intPatch.get();

    }else{
        m_intData = nullptr;
    }

}

/*!
    write mesh and data attached to class own dump file
*/
void IOWavefrontOBJ::dump(std::ostream & out){

    bitpit::utils::binary::write(out, m_tol);

    bool geoMark = (m_geometry != nullptr);
    bitpit::utils::binary::write(out, geoMark);
    if(geoMark) m_geometry->dump(out);

    // Use internal or external data object
    WavefrontOBJData* objData = nullptr;
    if (m_extData != nullptr)
    	objData = m_extData;
    else if (m_intData.get() != nullptr)
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
            std::map<long, std::string> locmap(pidlist.begin(), pidlist.end());
            // for(auto & entry : pidlist){
            //     locmap.insert({{entry.first, entry.second}});
            // }
            for(auto & entry : locmap){
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

    std::string line,key("mtllib"), result("");
    bool breakloop = false;
    while(in && !breakloop){
        std::getline(in, line);
        std::size_t fsplit = line.find(key);
        if(fsplit == std::string::npos) continue;

        result = line.substr(fsplit+key.length());
        result = bitpit::utils::string::trim(result);
        breakloop = true;
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
    \param[in,out] mapVCountObject map pid vs 3elements array containing number of v, vt,and vn for each sub-part.
    \param[out] nVertTot total vertices in the obj file
    \param[out] nCellTot total cells in the obj file
*/
void IOWavefrontOBJ::searchObjectPosition(std::ifstream & in,
                          std::vector<std::streampos> & mapPos,
                          std::vector<std::string>& mapNames,
                          std::vector<std::array<long,3>>& mapVCountObject,
                          std::array<long,3>& mapVCountTotal,
                          long &nCellTot)
{
    mapVCountTotal.fill(0);
    nCellTot = 0;
    if(!in.good()){
        in.clear();
        in.seekg(0);
    };

    std::string line;
    char key;
    std::array<long,3> vCounter({{0,0,0}});
    std::vector<std::array<long,3> > workingCounter;
    while(in){
        std::getline(in, line);
        if(line.length()< 2) continue; //ignore lines with less then 2 characters into.

        key = line.at(0);
        if(key == 'o'){
            mapPos.push_back(in.tellg());
            mapNames.push_back(line.substr(1));
            workingCounter.push_back(vCounter);
            vCounter.fill(0);
        }
        else if(key == 'v' && line.at(1) == ' '){
            ++mapVCountTotal[0];
            ++vCounter[0];
        }
        else if(key == 'v' && line.at(1) == 't'){
            ++mapVCountTotal[1];
            ++vCounter[1];
        }
        else if(key == 'v' && line.at(1) == 'n'){
            ++mapVCountTotal[2];
            ++vCounter[2];
        }
        else if(key == 'f') {
             ++nCellTot;
        }
    }
    workingCounter.push_back(vCounter);

    //stash the first useless entry of working counter;
    std::size_t sizeWC = workingCounter.size();
    if(sizeWC > 1){
        mapVCountObject.resize(sizeWC-1);
        memcpy(mapVCountObject.data(), &workingCounter[1], 3*(sizeWC-1)*sizeof(long));
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
    m_intData(for attached WavefrontOBJData datatypes)
    \param[in] in reading stream
    \param[in] begObjectStream position of the stream to start absorbing data
    \param[in] PID part identifier assigned to the current object.
    \param[in] vCounter 3 element array carryong the number of v, vt and vn currently present in this sub-part.
    \param[in, out] vOffset i/o current id of the next-to-be inserted mesh vertex.
    \param[in, out] vnOffset i/o current id of the next-to-be inserted mesh vertex normal.
    \param[in, out] vTxtOffset i/o current id of the next-to-be inserted texture vertex.
    \param[in, out] cOffset i/o current id of the next-to-be inserted mesh cell (same for mesh and texture)
*/
void IOWavefrontOBJ::readObjectData(std::ifstream & in, const std::streampos &begObjectStream, const long &PID,
                     const std::array<long,3> & vCounter, long &vOffset, long &vnOffset, long &vTxtOffset, long &cOffset)
{
    if (!in.good()){
        in.clear();
    }
    in.seekg(begObjectStream);
    std::string line;
    char key = ' ';

    std::string activeMaterial(""), stringsmooth;
    long activeSmooth = 0;
    //std::vector<long> locConn;
    std::string connEntry, txtEntry, normalEntry;

    MimmoObject * textures = m_intData->textures.get();
    MimmoObject * normals = m_intData->normals.get();


    // the idea is v, vt, vn are always compacted in block and not sparsed in the file.
    // I have already calculated the number of lines for each one and stored in the input vCounter[0,1,2],
    // so i can read the block of v lines (that are vCounter[0] lines) without controlling the typ of line each time.
    // same scheme for vt and vn if any.
    // connectivity info and other stuff are sparse. so i need the switch control on them.
    std::stringstream ss;
    std::string dummy, dummy2, cast_stod;
    long connvalue;

    int facetdefinition = -1;
    // this variable is set once the first occurrence of facet is read in this object chunk.
    // no need to check it every facet line. I'm sure all the facet are written the same.
    // cases are:
    // -1: undefined. Something wrong here.
    //  0: v v v           (facet made by only vertices - classic connectivity)
    //  1: v/vt v/vt v/vt   (texture prop attached. Be aware vt index has uniquely correspondence with v index)
    //  2: v/vt/vn v/vt/vn v/vt/vn (complete form)
    //  3: v//vn v//vn v//vn (normal prop attached. Be aware vn index has uniquely correspondence with v index)

    while(in.good() && key != 'o'){
        std::getline(in, line);

        if(line.length() < 2) continue; //ignore all lines with less then 2 characters into.
        key = line.at(0);

        if(key == 'v' && line.at(1) == ' '){
            //absorb mesh vertices for a number of lines = vCounter[0].
            std::array<double,3> temp;
            for(long i=0; i<vCounter[0]; ++i){
                ss.str(line);
                temp.fill(0.0);
                ss>>dummy;
                for(double & val : temp){
                    ss>>cast_stod;
                    val = std::stod(cast_stod);
                }
                m_intPatch->addVertex(temp, vOffset + i );
                std::getline(in, line);
                key = line.at(0);
                ss.clear();
            }
            vOffset += vCounter[0];
        }

        if(key == 'v' && line.at(1) == 't'){
            //absorb mesh textures for a number of lines = vCounter[1].
            std::array<double,3> temp;
            for(long i=0; i<vCounter[1]; ++i){
                ss.str(line);
                temp.fill(0.0);
                ss>>dummy;
                std::array<double,3>::iterator it_temp = temp.begin();
                while(ss.good() && it_temp != temp.end()){
                    ss>>cast_stod;
                    *it_temp = std::stod(cast_stod);
                    ++it_temp;
                }
                textures->addVertex(temp, vTxtOffset + i );
                std::getline(in, line);
                key = line.at(0);
                ss.clear();
            }
            vTxtOffset += vCounter[1];
        }

        if(key == 'v' && line.at(1) == 'n'){
            //absorb vertex normals for a number of lines = vCounter[2].
            std::array<double,3> temp;
            for(long i=0; i<vCounter[2]; ++i){
                ss.str(line);
                temp.fill(0.0);
                ss>>dummy;
                for(double & val: temp){
                    ss>>cast_stod;
                    val = std::stod(cast_stod);
                }
                normals->addVertex(temp, vnOffset + i );
                std::getline(in, line);
                key = line.at(0);
                ss.clear();
            }
            vnOffset += vCounter[2];
        }

        switch( convertKeyEntryToInt(key) ){
            case 0: //facet cell fn
                {
                    // be sure to have a facet definition.
                    if(facetdefinition < 0){
                        ss.str(line);
                        ss>>dummy>>dummy2;
                        ss.clear();

                        facetdefinition = checkFacetDefinition(dummy2);
                        if(facetdefinition<0) break;
                    }

                    std::replace(line.begin(), line.end(), '/', ' ');
                    ss.str(line);
                    //resize the local connectivity.
                    std::vector<long> locConn, txtConn, normalConn;
                    //read first f and first connectivity entry.
                    ss>>dummy>>connEntry;
                    while (ss.good()){
                        connvalue = std::stol(connEntry);
                        if(connvalue < 0) connvalue += vOffset;
                        locConn.push_back(connvalue);

                        switch(facetdefinition){
                            case 1:
                                ss>> txtEntry;
                                connvalue = std::stol(txtEntry);
                                if(connvalue < 0) connvalue += vTxtOffset;
                                txtConn.push_back(connvalue);
                                break;
                            case 2:
                                ss>> txtEntry >> normalEntry;
                                connvalue = std::stol(txtEntry);
                                if(connvalue < 0) connvalue += vTxtOffset;
                                txtConn.push_back(connvalue);

                                connvalue = std::stol(normalEntry);
                                if(connvalue < 0) connvalue += vnOffset;
                                normalConn.push_back(connvalue);

                                break;
                            case 3:
                                ss >> normalEntry;
                                connvalue = std::stol(normalEntry);
                                if(connvalue < 0) connvalue += vnOffset;
                                normalConn.push_back(connvalue);
                                break;
                            default:
                                //do nothing
                                break;
                        }//end switch
                        //advance stream and get rubbish of eof
                        ss>>connEntry;
                    }//end while
                    ss.clear();

                    //adding cell desuming cell type from locConn.
                    if(pushCell(m_intPatch.get(), locConn, PID, cOffset, -1) == bitpit::Cell::NULL_ID){
                        *(m_log)<<"WARNING "<<m_name<<" : skipping unsupported facet while reading obj file. "<<std::endl;
                    }else{
                        if(textures)    pushCell(textures, txtConn, PID, cOffset, -1);
                        if(normals)     pushCell(normals, normalConn, PID, cOffset, -1);

                        //adjusting data;
                        m_intData->materials.insert(cOffset, activeMaterial);
                        m_intData->smoothids.insert(cOffset, activeSmooth);
                        //increment the final coffset
                        ++cOffset;
                    }
                }
                break;
            case 1: //usemtl material name
                ss.str(line);
                ss>>dummy>> activeMaterial;
                ss.clear();
                break;
            case 3: //s smoothgroup id
                stringsmooth = "";
                ss.str(line);
                ss>>dummy >>stringsmooth;
                stringsmooth = bitpit::utils::string::trim(stringsmooth);
                if(stringsmooth.empty() || stringsmooth == "off"){
                    activeSmooth = 0;
                }else{
                    activeSmooth = std::stol(stringsmooth);
                }
                ss.clear();
                break;
            case 99: // o entry
                //DO NOTHING.
                break;
            case 2: //g cellgroup string labels
            case 4: //l conn line
            case 5: //p conn point
            default: //unsupported flag
                *(m_log)<<"WARNING "<<m_name<<" : unsupported flag "<<key<<" declaration while reading obj file. Ignoring..."<<std::endl;
                break;
        }

    } //end of stream chunk reading.

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
    m_extData(for attached WavefrontOBJData datatypes)
    \param[in] objData source object data structure
    \param[in] out writing stream
    \param[in] vertexLists list of points involved for the object for all 3 type v, vt, vn
    \param[in] cellList list of cells involved for the object
    \param[in,out] vOffsets offset for local v,vt,vn counting
    \param[in,out] cOffset offset for local facet/cell counting
    \param[in,out] vinsertion_maps local map id-written insertion index for v,vt,vn (to avoid rewrite vertices)
    \param[in,out] activeMaterial currently active material for the mesh
    \param[in,out] activeSmoothIds currently active smoothids for the mesh
*/
void IOWavefrontOBJ::writeObjectData(WavefrontOBJData* objData, std::ofstream & out,
                                     const std::array<std::vector<long>,3> & vertexLists,
                                     const std::vector<long> & cellList,
                                     std::array<long,3> &vOffsets,
                                     long &cOffset,
                                     std::array<std::unordered_map<long,long>,3> & vinsertion_maps,
                                     long & activeMaterial, long &activeSmoothId)
{
    int facetType = 0; //only mesh present;
    if(objData->textures.get() && objData->normals.get())   facetType = 2;
    if(objData->textures.get() && !objData->normals.get())   facetType = 1;
    if(!objData->textures.get() && objData->normals.get())   facetType = 3;

    //writing mesh vertices
    for(long id: vertexLists[0]){
        if(vinsertion_maps[0].count(id) > 0) continue;

        bitpit::Vertex & point = m_geometry->getVertices().at(id);
        out<<"v "<<std::fixed<<std::setprecision(6)<<point[0]<<" "<<point[1]<<" "<<point[2]<<'\n';
        vinsertion_maps[0].insert({{id, vOffsets[0]}});
        ++vOffsets[0];
    }
    //write texture vertices
    for(long id: vertexLists[1]){
        if(vinsertion_maps[1].count(id) > 0) continue;

        bitpit::Vertex &point = objData->textures->getVertices().at(id);
        out<<"vt "<<std::fixed<<std::setprecision(6)<<point[0]<<" "<<point[1]<<" ";
        if(!m_textureUVMode)    out<<std::fixed<<std::setprecision(6)<<point[2];
        out<<'\n';
        vinsertion_maps[1].insert({{id, vOffsets[1]}});
        ++vOffsets[1];
    }
    //write vnormals
    for(long id: vertexLists[2]){
        if(vinsertion_maps[2].count(id) > 0) continue;

        bitpit::Vertex &point = objData->normals->getVertices().at(id);
        out<<"vn "<<std::fixed<<std::setprecision(6)<<point[0]<<" "<<point[1]<<" "<<point[2]<<'\n';
        vinsertion_maps[2].insert({{id, vOffsets[2]}});
        ++vOffsets[2];
    }

    //connectivity time. Use insertion maps and regroup by cellgroups-materials-smoothids if any in m_extData.
    MaterialGroups gCells = regroupCells(objData, cellList);

    //enter the material key subtree
    for(auto & matList : gCells){
        //write material
        if(activeMaterial != matList.first){
            activeMaterial = matList.first;
            out<<"usemtl "<<objData->inv_materialsList.at(matList.first)<<'\n';
        }
        // access and fnally write the chunk of facets/cells
        for(long idCell : matList.second){
            //write smoothing id
            long currentsid = objData->smoothids.at(idCell);
            if(activeSmoothId != currentsid){
                activeSmoothId = currentsid;
                out<<"s "<<objData->inv_smoothidsList.at(currentsid)<<'\n';
            }


            // prepare to write facets
            bitpit::ConstProxyVector<long> meshIds = m_geometry->getCells()[idCell].getVertexIds();
            std::size_t countV = meshIds.size(); //it should be the same for  textures and normals also

            switch(facetType){
                case 1: //v and vt
                    {
                        bitpit::ConstProxyVector<long> txtIds= objData->textures->getCells().at(idCell).getVertexIds();
                        out<<"f ";
                        for(int i=0; i<countV; ++i){
                            out<<vinsertion_maps[0].at(meshIds[i])<<"/";
                            out<<vinsertion_maps[1].at(txtIds[i])<<" ";
                        }
                        out<<'\n';
                    }

                    break;
                case 2: //v vt and vn
                    {
                        bitpit::ConstProxyVector<long> txtIds= objData->textures->getCells().at(idCell).getVertexIds();
                        bitpit::ConstProxyVector<long> normalIds= objData->normals->getCells().at(idCell).getVertexIds();
                        out<<"f ";
                        for(int i=0; i<countV; ++i){
                            out<<vinsertion_maps[0].at(meshIds[i])<<"/";
                            out<<vinsertion_maps[1].at(txtIds[i])<<"/";
                            out<<vinsertion_maps[2].at(normalIds[i])<<" ";
                        }
                        out<<'\n';

                    }

                    break;
                case 3: //v and vn
                    {
                        bitpit::ConstProxyVector<long> normalIds= objData->normals->getCells().at(idCell).getVertexIds();
                        out<<"f ";
                        for(int i=0; i<countV; ++i){
                            out<<vinsertion_maps[0].at(meshIds[i])<<"//";
                            out<<vinsertion_maps[2].at(normalIds[i])<<" ";
                        }
                        out<<'\n';
                    }
                    break;
                default: // v only
                        out<<"f ";
                        for(int i=0; i<countV; ++i){
                            out<<vinsertion_maps[0].at(meshIds[i])<<" ";
                        }
                        out<<'\n';
                    break;
            } //end switch
        }//ending leaflist - chunk of cells.
    }  // ending materials

}

/*!
    Regroup a list of mesh cells by materials
    \param[in] objData source object data structure
    \param[in] cellList list of cell ids of the polygonal mesh
    \return tree containing cell sublists reordered as [Material] -> sublist.
*/
IOWavefrontOBJ::MaterialGroups
IOWavefrontOBJ::regroupCells(const WavefrontOBJData* objData, const livector1D & cellList)
{

    // materials and cellgroup need to be fully compiled, even with the default flags.
    MaterialGroups map;
    long mat_key;
    for(long id: cellList){
        mat_key = objData->materialsList.at(objData->materials[id]);
        if(map.count(mat_key) == 0){
            map[mat_key].reserve(cellList.size());
        }
        map[mat_key].push_back(id);
    }

    for(auto & entrymap : map){
        entrymap.second.shrink_to_fit();
    }

    return map;
}

/*!
    convert a key entry in a int related to connectivity.
    The vertex entry keys (v,vn,vt,vp) are not considered.
    useful for reading stage.

*/
int IOWavefrontOBJ::convertKeyEntryToInt(char key){

    int res = -1;
    if(key == 'f')  res=0;
    if(key == 'u')  res=1; //usemtl
    if(key == 'g')  res=2;
    if(key == 's')  res=3;
    if(key == 'l')  res=4;
    if(key == 'p')  res=5;
    if(key == 'o')  res=99;
    return res;
}

/*!
 * Compute the vertex normals of the vertices involved by geometry displacements.
 * Use geometry tolerance as threshold on the norm of the displacement vector
 * to define if a vertex is moved.
 */


/*!
    Internal utility to push a cell in a m_intPatch. Wrapping MimmoObject::addConnectedCell.
    return the id of newly inserted cell if successfull, bitpit::Cell::NULL_ID if failed.
    \param[in] obj pointer of mesh where the cell need to be pushed
    \param[in] conn connectivity array
    \param[in] PID part identifier on the cell
    \param[in] id unique id of the cell
    \param[in] rank MPI only -  proc rank.
*/

long IOWavefrontOBJ::pushCell(MimmoObject * obj, std::vector<long> &conn, long PID, long id, int rank ){
    long markedid = bitpit::Cell::NULL_ID;

    switch(int(conn.size())){
    case 0:
    case 1:
    case 2:
        //not allowed cases
        break;
    case 3: //triangles
        markedid = obj->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, PID, id, -1);
        break;
    case 4: //quads
        markedid = obj->addConnectedCell(conn, bitpit::ElementType::QUAD, PID, id, -1);
        break;
    default: //polygons
        {
            std::vector<long> tt(1,conn.size());
            tt.insert(tt.end(), conn.begin(), conn.end());
            markedid =  obj->addConnectedCell(tt, bitpit::ElementType::POLYGON, PID, id, -1);
        }
        break;
    }

    return markedid;
}

/*!
    \return facetdefinition according to the legend.
    Legend:
    - 1: undefined. Something wrong here.
    - 0: v v v           (facet made by only vertices - classic connectivity)
    - 1: v/vt v/vt v/vt   (texture prop attached. Be aware vt index has uniquely correspondence with v index)
    - 2: v/vt/vn v/vt/vn v/vt/vn (complete form)
    - 3: v//vn v//vn v//vn (normal prop attached. Be aware vn index has uniquely correspondence with v index)
     \param[in] target string read xxx/xxx/xxxx
     \param[out] number of valid elements in the target string.
*/
int IOWavefrontOBJ::checkFacetDefinition(const std::string & str){

    int definition = -1;
    if(str.empty()){
        return definition;
    }

    definition = std::count(str.begin(),str.end(), '/');
    if (definition == 2){
        definition += countSubstring(str, "//");
    }
    return definition;
}


}
