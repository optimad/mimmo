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
 * Basic Constructor
 */
SelectionByPID::SelectionByPID(){
    m_name = "mimmo.SelectionByPID";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByPID::SelectionByPID(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.SelectionByPID";

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.SelectionByPID"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Custom Constructor
 * \param[in] pidlist List of pid to be included into selection
 * \param[in] target    Pointer to target geometry
 */
SelectionByPID::SelectionByPID(shivector1D & pidlist, MimmoObject * target){
    m_name = "mimmo.SelectionByPID";
    setGeometry(target);
    setPID(pidlist);
};

/*!
 * Destructor
 */
SelectionByPID::~SelectionByPID(){};

/*!
 * Copy constructor
 */
SelectionByPID::SelectionByPID(const SelectionByPID & other):GenericSelection(other){
    m_activePID = other.m_activePID;
};

/*!
 * Copy Operator
 */
SelectionByPID & SelectionByPID::operator=(const SelectionByPID & other){
    *(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection * >(&other));
    m_activePID = other.m_activePID;
    return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void
SelectionByPID::buildPorts(){

    bool built = true;

    GenericSelection::buildPorts();

    built = (built && createPortIn<short, SelectionByPID>(this, &SelectionByPID::setPID, PortType::M_VALUESI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::SHORT));
    built = (built && createPortIn<std::vector<short>, SelectionByPID>(this, &SelectionByPID::setPID, PortType::M_VECTORSI, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::SHORT));

    m_arePortsBuilt = built;
};


/*!
 * Return all available PID for current geometry linked. Void list means not linked target geometry
 * \return vector with PIDs for current geometry
 */
shivector1D
SelectionByPID::getPID(){
    shivector1D result(m_activePID.size());
    int counter = 0;
    for(auto && val : m_activePID){
        result[counter] = val.first;
        ++counter;
    }
    return     result;
};
/*!
 * Return active/inactive PIDs for current geometry linked.
 * \param[in] active boolean to get active-true, inactive false PID for selection
 * \return active/inactive PIDs for current geometry linked.
 */
shivector1D
SelectionByPID::getActivePID(bool active){
    shivector1D result;
    for(auto && val : m_activePID){
        if(val.second == active)    result.push_back(val.first);
    }
    return(result);
};


/*!
 * Set pointer to your target geometry. Reimplemented from mimmo::BaseManipulation::setGeometry()
 * \param[in] target Pointer to MimmoObject target geometry.
 */
void
SelectionByPID::setGeometry(MimmoObject * target ){
    if(target == NULL) return;
    m_geometry = target;

    std::unordered_set<short> & pids =     target->getPIDTypeList();
    std::unordered_set<short>::iterator it, itE = pids.end();

    for(it = pids.begin(); it!=itE; ++it){
        m_activePID.insert(std::make_pair(*it,false));
    }
};

/*!
 * Activate flagged PID i. If i<0, activates all PIDs available.
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] i PID to be activated.
 */
void
SelectionByPID::setPID(short i){
    if(m_setPID.count(-1) >0 || (!m_setPID.empty() && i==-1))    m_setPID.clear();
    m_setPID.insert(i);

};

/*!
 * Activate a list of flagged PID i. SelectionByPID::setPID(short i).
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] list    list of PID to be activated
 */
void
SelectionByPID::setPID(shivector1D list){
    for(auto && index : list){
        setPID(index);
    }
};
/*!
 * Deactivate flagged PID i. If i<0, deactivates all PIDs available.
 * if i > 0 but does not exist in PID available list, do nothing
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] i PID to be deactivated.
 */

void
SelectionByPID::removePID(short i){
    if(i>0){
        if(m_setPID.count(i) >0) m_setPID.erase(i);
    }else{
        m_setPID.clear();
    }
};

/*!
 * Deactivate a list of flagged PID i. SelectionByPID::removePID(short i = -1).
 * Modifications are available after applying them with syncPIDList method; 
 * \param[in] list    list of PID to be deactivated
 */
void
SelectionByPID::removePID(shivector1D list){
    for(auto && index : list){
        removePID(index);
    }
};

/*!
 * Clear your class
 */
void
SelectionByPID::clear(){
    m_subpatch.reset(nullptr);
    m_activePID.clear();
    BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tessellation extracted
 */
livector1D
SelectionByPID::extractSelection(){
    livector1D result;
    std::set<long> extraction;

    syncPIDList();
    shivector1D pids = getActivePID();

    for(auto && pid : pids){
        livector1D local = getGeometry()->extractPIDCells(pid);
        extraction.insert(local.begin(),local.end());
    }

    //check if dual selection is triggered

    if(m_dual){
        livector1D totID = getGeometry()->getMapCell();
        result.resize(totID.size() - extraction.size());
        if(result.size() == 0) return result;

        std::sort(totID.begin(), totID.end());
        int counter = 0;
        auto tot_it  = totID.begin();
        auto cand_it = extraction.begin();
        while(tot_it != totID.end()){
            long val = *tot_it;
            if (cand_it == extraction.end() || val != *cand_it) {
                result[counter] = val;
                ++counter;
            } else {
                ++cand_it;
            }
            ++tot_it;
        }
    }else{
        result.insert(result.end(), extraction.begin(), extraction.end());
    }
    return    result;
};

/*!
 * Checks & synchronizes user given pid list w/ current pid list available by linked geometry.
 * Activate positive matching PIDs and erase negative matches from user list. 
 * if geometry is not currently linked does nothing.
 */
void
SelectionByPID::syncPIDList(){
    if(getGeometry() == NULL)    return;

    if(m_setPID.count(-1) == 0){

        std::unordered_set<short>::iterator itU;
        shivector1D negative;
        for(itU = m_setPID.begin(); itU != m_setPID.end(); ++itU){
            short value = *itU;
            if(m_activePID.count(value) > 0)    m_activePID[value] = true;
            else    negative.push_back(value);
        }

        for(auto &val: negative)    m_setPID.erase(val);


    }else{
        m_setPID.clear();
        for(auto &val: m_activePID ){
            m_setPID.insert(val.first);
            val.second = true;

        }
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByPID::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    //start absorbing
    BaseManipulation::absorbSectionXML(slotXML, name);

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

    int nPID = 0;
    if(slotXML.hasOption("nPID")){
        std::string input = slotXML.get("nPID");
        input = bitpit::utils::string::trim(input);
        nPID = 0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>nPID;
            std::max(0, nPID);
        }
    }
    shivector1D pidlist(nPID, -1);

    if(slotXML.hasOption("PID")){
        std::string input = slotXML.get("PID");
        input = bitpit::utils::string::trim(input);
        if(!input.empty()){
            std::stringstream ss(input);
            for(auto &val: pidlist){
                ss>>val;
            }
        }
    }

    setPID(pidlist);

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByPID::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::flushSectionXML(slotXML, name);
    int value = m_dual;
    slotXML.set("Dual", std::to_string(value));


    shivector1D selected = getActivePID(true);
    int size = selected.size();


    if(size != 0){
        slotXML.set("nPID", std::to_string(size));

        std::stringstream ss;
        for(int i=0; i<size; ++i){
            ss<<selected[i]<<'\t';
        }

        slotXML.set("PID", ss.str());
    }

};

}

