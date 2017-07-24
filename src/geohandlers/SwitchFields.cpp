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

#include "SwitchFields.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!Default constructor of SwitchField.
 */
SwitchField::SwitchField(){
    m_mapping = false;
}

/*!
 * Default destructor of SwitchField.
 */
SwitchField::~SwitchField(){
    clear();
};

/*!Copy constructor of SwitchField.Soft Copy of MimmoObject;
 */
SwitchField::SwitchField(const SwitchField & other):BaseManipulation(){
    *this = other;
};

/*!
 * Assignement operator of SwitchField. Soft copy of MimmoObject
 */
SwitchField & SwitchField::operator=(const SwitchField & other){
    clear();
    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_mapping = other.m_mapping;
    return *this;
};

/*!
 * Build the ports of the class;
 */
void
SwitchField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, SwitchField>(this, &mimmo::SwitchField::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    m_arePortsBuilt = built;
}


/*!
 * Set the target geometry where your target switched field is defined.
 * \param[in] geo  Pointer to MimmoObject
 */
void
SwitchField::setGeometry(MimmoObject* geo){
    if(geo->isEmpty()) return;
    m_geometry = geo;
};

/*!
 * Set if force the research by mapping after the first guess.
 * \param[in] flag  True force the second research
 */
void
SwitchField::setMapping(bool flag){
    m_mapping = flag;
};

/*!
 * Clear all stuffs in your class
 */
void
SwitchField::clear(){
    BaseManipulation::clear();
};

/*!Execution command.
 * Switch a field over a target geometry from a vector of fields.
 */
void
SwitchField::execute(){

    bool check = mswitch();
    if(!check){
        (*m_log)<<"error in class "<<m_name<<". Field cannot be switch"<<std::endl;
        (*m_log)<<"This could be due to not correct setting of geometries or division maps"<<std::endl;
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SwitchField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);


    if(slotXML.hasOption("Mapping")){
        std::string input = slotXML.get("Mapping");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setMapping(value);
    }

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SwitchField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

    if(m_mapping){
        slotXML.set("Mapping", std::to_string(1));
    }

};

}
