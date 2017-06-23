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

#include "ExtractFields.hpp"

using namespace std;
using namespace bitpit;
namespace mimmo{

/*!Default constructor of ExtractField.
 */
ExtractField::ExtractField(){
    m_mode = ExtractMode::ID;
}

/*!
 * Default destructor of ExtractField.
 */
ExtractField::~ExtractField(){
    clear();
};

/*!Copy constructor of ExtractField.Soft Copy of MimmoObject;
 */
ExtractField::ExtractField(const ExtractField & other):BaseManipulation(){
    *this = other;
};

/*!
 * Assignement operator of ExtractField. Soft copy of MimmoObject
 */
ExtractField & ExtractField::operator=(const ExtractField & other){
    clear();
    *(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
    m_mode = other.m_mode;
    return *this;
};

/*!
 * Build the ports of the class;
 */
void
ExtractField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, ExtractField>(this, &mimmo::ExtractField::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_, true));
    m_arePortsBuilt = built;
}

/*!
 * Set the target geometry where your target extracted field is defined.
 * \param[in] geo  Pointer to MimmoObject
 */
void
ExtractField::setGeometry(MimmoObject* geo){
    if(geo->isEmpty()) return;
    m_geometry = geo;
};

/*!
 * Set mode of extraction.
 * \param[in] mode Extraction mode
 */
void
ExtractField::setMode(ExtractMode mode){
    setMode(static_cast<int>(mode));
};

/*!
 * Set mode of extraction.
 * \param[in] mode Extraction mode
 */
void
ExtractField::setMode(int mode){
    if(mode <1 ||mode > 3)    return;
    m_mode = static_cast<ExtractMode>(mode);
};

/*!
 * Clear all stuffs in your class
 */
void
ExtractField::clear(){
    BaseManipulation::clear();
};

/*!Execution command.
 * Extract a field over a target geometry from an input field.
 */
void
ExtractField::execute(){

    bool check = extract();
    if(!check){
        (*m_log)<<"error in class "<<m_name<<". Field cannot be extract"<<std::endl;
        (*m_log)<<"This could be due to not correct setting of geometries or division maps"<<std::endl;
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ExtractField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);

    std::string input;

    if(slotXML.hasOption("ExtractMode")){
        input = slotXML.get("ExtractMode");
        input = bitpit::utils::trim(input);
        int value = 1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::min(std::max(1, value),4);
        }
    };

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ExtractField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);


    int value = static_cast<int>(m_mode);
    slotXML.set("ExtractMode", std::to_string(value));

};

}
