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
 \ *---------------------------------------------------------------------------*/

#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
namespace mimmo{

/*!
 * Basic Constructor. Need to know kind of topology chosen 1-3D surface,
 * 2-VolumeMesh. Other options are not available,
 * if forced trigger default value of 1.
 */
SelectionByMapping::SelectionByMapping(int topo){
    m_name = "mimmo.SelectionByMapping";
    m_type = SelectionType::MAPPING;
    m_tolerance = 1.E-08;

    topo = std::max(1,topo);
    if(topo > 2)    topo=1;
    m_topo = topo;


    m_allowedType.resize(3);
    m_allowedType[1].insert(FileType::STL);
    m_allowedType[1].insert(FileType::STVTU);
    m_allowedType[1].insert(FileType::SQVTU);
    m_allowedType[1].insert(FileType::NAS);
    m_allowedType[2].insert(FileType::VTVTU);
    m_allowedType[2].insert(FileType::VHVTU);

};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByMapping::SelectionByMapping(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SelectionByMapping";
    m_type = SelectionType::MAPPING;
    m_tolerance = 1.E-08;

    std::string fallback_name = "ClassNONE";
    std::string fallback_topo = "-1";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    input_name = bitpit::utils::string::trim(input_name);

    std::string input_topo = rootXML.get("Topology", fallback_topo);
    input_topo = bitpit::utils::string::trim(input_topo);

    int topo = std::stoi(input_topo);
    topo = std::max(1,topo);
    if(topo > 2)    topo=1;
    m_topo = topo;

    m_allowedType.resize(3);
    m_allowedType[1].insert(FileType::STL);
    m_allowedType[1].insert(FileType::STVTU);
    m_allowedType[1].insert(FileType::SQVTU);
    m_allowedType[1].insert(FileType::NAS);
    m_allowedType[2].insert(FileType::VTVTU);
    m_allowedType[2].insert(FileType::VHVTU);


    if(input_name == "mimmo.SelectionByMapping"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor. Non null geometry set also the topology allowed bu the class.
 * \param[in]    geolist    List of external files to read from, for comparison purposes
 * \param[in]    target    Pointer to target geometry
 * \param[in]    tolerance    Proximity criterium tolerance
 */
SelectionByMapping::SelectionByMapping(std::unordered_map<std::string, int> & geolist, MimmoObject * target, double tolerance){
    m_name = "mimmo.SelectionByMapping";
    m_type = SelectionType::MAPPING;
    m_tolerance = 1.E-08;

    m_allowedType.resize(3);
    m_allowedType[1].insert(FileType::STL);
    m_allowedType[1].insert(FileType::STVTU);
    m_allowedType[1].insert(FileType::SQVTU);
    m_allowedType[1].insert(FileType::NAS);
    m_allowedType[2].insert(FileType::VTVTU);
    m_allowedType[2].insert(FileType::VHVTU);

    if(target != NULL){
        m_topo = target->getType();
        m_topo = std::min(1, m_topo);
        if(m_topo > 2)    m_topo = 1;
        setGeometry(target);
        setFiles(geolist);
        m_tolerance = tolerance;
    }

};

/*!
 * Destructor
 */
SelectionByMapping::~SelectionByMapping(){};

/*!
 * copy Constructor
 */
SelectionByMapping::SelectionByMapping(const SelectionByMapping & other):GenericSelection(other){
    m_tolerance = other.m_tolerance;
    m_geolist = other.m_geolist;
    m_mimmolist = other.m_mimmolist;
    m_allowedType = other.m_allowedType;
};

/*!
 * Copy Operator
 */
SelectionByMapping & SelectionByMapping::operator=(SelectionByMapping other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x SelectionByMapping object
 */
void SelectionByMapping::swap(SelectionByMapping & x) noexcept
{
    std::swap(m_tolerance, x.m_tolerance);
    std::swap(m_geolist, x.m_geolist);
    std::swap(m_mimmolist, x.m_mimmolist);
    std::swap(m_allowedType, x.m_allowedType);
    GenericSelection::swap(x);
}


/*!
 * It builds the input/output ports of the object
 */
void
SelectionByMapping::buildPorts(){

    bool built = true;

    GenericSelection::buildPorts();

    built = (built && createPortIn<MimmoObject *, SelectionByMapping>(this, &SelectionByMapping::addMappingGeometry,M_GEOM2));

    m_arePortsBuilt = built;
};

/*!
 * Return tolerance actually set in your class
 * \return tolerance set in your class
 */
double
SelectionByMapping::getTolerance(){
    return    m_tolerance;
};

/*!
 * Set Proximity tolerance.
 * Under this threshold compared geometries are considered coincident.
 * \param[in] tol Desired tolerance.
 */
void
SelectionByMapping::setTolerance(double tol){
    if(tol == 0.0){
        tol = 1.e-8;
    }
    m_tolerance = tol;
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * GenericSelection::setGeometry();
 * \param[in] target Pointer to target geometry.
 */
void
SelectionByMapping::setGeometry( MimmoObject * target){

    if(target->getType() != m_topo){
        (*m_log)<< " " << std::endl;
        (*m_log)<< m_name << " target Topology : "<< target->getType() << std::endl;
        (*m_log)<<" supported Topology : "<< m_topo << std::endl;
        (*m_log)<<" selectionMapping cannot support current geometry. Topology not supported."<<std::endl;
        (*m_log)<< " " << std::endl;
        return;
    }
    m_geometry = target;

};

/*!
 * Return the actual list of external geometry files to read from and compare your target geometry for mapping purposes
 * \return the list of external files to read
 */
const std::unordered_map<std::string, int> &
SelectionByMapping::getFiles() const{
    return    m_geolist;
}

/*!
 * Set a list of external geometry files to read from and compare your target geometry for mapping purposes
 * \param[in] files List external geometries to be read.
 */
void
SelectionByMapping::setFiles(std::unordered_map<std::string, int>  files){
    for(auto && val : files){
        addFile(val);
    }
};

/*!
 * Add a new file to a list of external geometry files
 * \param[in] file External file to be read
 */
void
SelectionByMapping::addFile(std::pair<std::string, int> file){
    int type = m_topo;
    if(m_allowedType[type].find(file.second) != m_allowedType[type].end()){
        m_geolist.insert(file);
    }
};

/*!
 * Add a new mimmo object to a list of external geometries (only surface mesh allowed)
 * \param[in] obj Pointer to MimmoObject of external geometry
 */
void
SelectionByMapping::addMappingGeometry(MimmoObject* obj){
    if(obj->getType() == m_topo){
        m_mimmolist.insert(obj);
    }
};

/*!
 * Remove an existent file to a list of external geometry files. If not in the list, do nothing 
 * \param[in] file Name of the file to be removed from the list
 */
void
SelectionByMapping::removeFile(std::string file){
     if(m_geolist.find(file) != m_geolist.end())    m_geolist.erase(file);
};

/*!
 * Empty your list of file for mapping 
 */
void
SelectionByMapping::removeFiles(){
    m_geolist.clear();
};

/*!
 * Empty your list of geometries for mapping
 */
void
SelectionByMapping::removeMappingGeometries(){
    m_mimmolist.clear();
};

/*!
 * Clear your class
 */
void
SelectionByMapping::clear(){
    m_subpatch.reset(nullptr);
    removeFiles();
    removeMappingGeometries();
    m_topo = 0;
    BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tessellation extracted
 */
livector1D
SelectionByMapping::extractSelection(){

    if(!(getGeometry()->isSkdTreeSync()))    getGeometry()->buildSkdTree();
    std::set<long> cellList;

    for (auto && file : m_geolist){
        livector1D list = getProximity(file);
        cellList.insert(list.begin(), list.end());
    }

    for (auto mimmo : m_mimmolist){
        livector1D list = getProximity(mimmo);
        cellList.insert(list.begin(), list.end());
    }

/* check if dual */
    livector1D result;
    if(m_dual){
        livector1D totID = getGeometry()->getCellsIds();
        result.resize(totID.size() - cellList.size());
        if(result.size() == 0) return result;

        std::sort(totID.begin(), totID.end());
        int counter = 0;
        auto tot_it  = totID.begin();
        auto cand_it = cellList.begin();
        while(tot_it != totID.end()){
            long val = *tot_it;
            if (cand_it == cellList.end() || val != *cand_it) {
                result[counter] = val;
                ++counter;
            } else {
                ++cand_it;
            }
            ++tot_it;
        }
    }else{
        result.insert(result.end(), cellList.begin(), cellList.end());
    }
    return    result;
};

/*!
 * Return portion of target geometry near to an external geometry
 * \param[in] val Pair with file of external geometry to be compared and
 * number of total raw points for level set evaluation
 */
livector1D
SelectionByMapping::getProximity(std::pair<std::string, int> val){

    svector1D info = extractInfo(val.first);

    MimmoGeometry * geo = new MimmoGeometry();
    geo->setIOMode(IOMode::READ);
    geo->setDir(info[0]);
    geo->setFilename(info[1]);
    geo->setFileType(val.second);
    geo->setBuildSkdTree(true);
    geo->execute();

    if(geo->getGeometry()->getNVertex() == 0 || geo->getGeometry()->getNCells() == 0 ){
        m_log->setPriority(bitpit::log::NORMAL);
        (*m_log)<< m_name << " failed to read geometry in SelectionByMapping::getProximity"<<std::endl;
        m_log->setPriority(bitpit::log::DEBUG);
        return livector1D();
    }
    livector1D result = mimmo::skdTreeUtils::selectByPatch(geo->getGeometry()->getSkdTree(), getGeometry()->getSkdTree(), m_tolerance);
    delete geo;
    geo=NULL;

    return    result;
};

/*!
 * Return portion of target geometry near to an external geometry
 * \param[in] obj Pointer to external geometry to be compared.
 */
livector1D
SelectionByMapping::getProximity(MimmoObject* obj){

    obj->buildSkdTree(true);

    if(obj->getNVertex() == 0 || obj->getNCells() == 0 ){
        m_log->setPriority(bitpit::log::NORMAL);
        (*m_log)<< m_name << " failed to read geometry in SelectionByMapping::getProximity"<<std::endl;
        m_log->setPriority(bitpit::log::DEBUG);
        return livector1D();
    }
    livector1D result = mimmo::skdTreeUtils::selectByPatch(obj->getSkdTree(), getGeometry()->getSkdTree(), m_tolerance);

    return    result;
};

/*!
 * Extract root dir/filename/tag from an absolute file pattern
 * \return vector with extracted dir, filename and tag
 */
svector1D
SelectionByMapping::extractInfo(std::string file){

    std::string root, name, tag,temp;
    std::string key1=".", key2="/\\";

    std::size_t found = file.find_last_of(key2);
    root = file.substr(0, found);
    temp = file.substr(found+1);

    found = temp.find_last_of(key1);
    name = temp.substr(0,found);
    tag = temp.substr(found+1);

    svector1D result(3);
    result[0] = root;
    result[1] = name;
    result[2] = tag;

    return     result;
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByMapping::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    //checking topology
    if(slotXML.hasOption("Topology")){
        std::string input = slotXML.get("Topology");
        input = bitpit::utils::string::trim(input);
        int temp = -1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        if(m_topo != temp) {
            throw std::runtime_error (m_name + " : xml absorbing failed.");
        }
    }

    BaseManipulation::absorbSectionXML(slotXML, name);
    //start absorbing
    if(slotXML.hasOption("Dual")){
        std::string input = slotXML.get("Dual");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setDual(value);
    }

    if(slotXML.hasOption("Tolerance")){
        std::string input = slotXML.get("Tolerance");
        input = bitpit::utils::string::trim(input);
        double temp = 1.0E-8;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
            setTolerance(temp);

    }

    std::unordered_map<std::string, int> mapp;
    if(slotXML.hasSection("Files")){

        const bitpit::Config::Section & filesXML = slotXML.getSection("Files");

        for(auto & subfile : filesXML.getSections()){
            std::string path;
            std::string tag;

            if(subfile.second->hasOption("fullpath"))    {
                path = subfile.second->get("fullpath");
                path = bitpit::utils::string::trim(path);
            }
            if(subfile.second->hasOption("tag")){
                tag = subfile.second->get("tag");
                tag = bitpit::utils::string::trim(tag);
                //check tag;
                auto maybe_tag = FileType::_from_string_nothrow(tag.c_str());
                if(!maybe_tag)    tag.clear();
                else    tag = maybe_tag->_to_string();
            }

            if(!path.empty() && !tag.empty()){
                mapp[path] = (int) FileType::_from_string(tag.c_str());
            }
        }

        setFiles(mapp);

    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByMapping::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    
    slotXML.set("Topology", m_topo);

    int value = m_dual;
    slotXML.set("Dual", std::to_string(value));


    if(m_tolerance != 1.E-08){
        std::stringstream ss;
        ss<<std::scientific<<m_tolerance;
        slotXML.set("Tolerance", ss.str());
    }

    bitpit::Config::Section & filesXML = slotXML.addSection("Files");

    int counter = 0;
    for(auto & file : m_geolist){
        std::string name = "file"+std::to_string(counter);
        bitpit::Config::Section & local = filesXML.addSection(name);
        local.set("fullpath", file.first);
        std::string typetag = (FileType::_from_integral(file.second))._to_string();
        local.set("tag", typetag);
        ++counter;
    }

};

}
