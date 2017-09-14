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
    m_tol  = 1.0e-08;
}

/*!
 * Default destructor of ExtractField.
 */
ExtractField::~ExtractField(){
    clear();
};

/*!
 * Swap function of ExtractField
 * \param[in] x object to be swapped.
 */
void ExtractField::swap(ExtractField & x ) noexcept
{
    std::swap(m_mode,x.m_mode);
    std::swap(m_tol, x.m_tol);
    BaseManipulation::swap(x);
};

/*!
 * Build the ports of the class;
 */
void
ExtractField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<MimmoObject*, ExtractField>(this, &mimmo::ExtractField::setGeometry, M_GEOM, true));
    m_arePortsBuilt = built;
}

/*!
 * Set the target geometry where your target extracted field will be defined.
 * \param[in] geo  Pointer to MimmoObject
 */
void
ExtractField::setGeometry(MimmoObject* geo){
    if(geo == NULL)     return;
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
    if(mode < 1 ||mode > 3)    return;
    m_mode = static_cast<ExtractMode>(mode);
};

/*!
 * Set tolerance for extraction by patch.
 * \param[in] tol Tolerance
 */
void
ExtractField::setTolerance(double tol){
    m_tol = std::max(1.0e-12, tol);
};

/*!
 * Get tolerance for extraction by patch actually set in the class.
 * \return tolerance
 */
double
ExtractField::getTolerance(){
    return m_tol;
};

/*!
 * Get current extraction mode set on the class.
 * \return type of extraction in enum ExtractMode
 */
ExtractMode
ExtractField::getMode(){
    return m_mode;
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
        input = bitpit::utils::string::trim(input);
        int value = 1;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
            value = std::min(std::max(1, value),4);
        }
    };

    if(slotXML.hasOption("Tolerance")){
        std::string input = slotXML.get("Tolerance");
        input = bitpit::utils::string::trim(input);
        double temp = 0.0;
        if(!input.empty()){
            std::stringstream ss(input);
            ss>>temp;
        }
        setTolerance(temp);
    }

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
    if (m_mode == ExtractMode::MAPPING)
        slotXML.set("Tolerance", std::to_string(m_tol));

};

}
