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

#include "SelectField.hpp"

namespace mimmo{

/*!Default constructor of SelectField.
 */
SelectField::SelectField(){
	m_mode = SelectType::NAME;
	m_fieldname = "data";
    m_tol = 1.0e-12;
}

/*!
 * Default destructor of SelectField.
 */
SelectField::~SelectField(){
    clear();
};

/*!Copy constructor of SelectField.
 */
SelectField::SelectField(const SelectField & other):BaseManipulation(other){
	m_mode = other.m_mode;
	m_fieldname = other.m_fieldname;
    m_tol = other.m_tol;
};

/*!
 * Swap function.
 * \param[in] x object to be swapped
 */
void SelectField::swap(SelectField & x) noexcept
{
    std::swap(m_mode, x.m_mode);
    std::swap(m_fieldname, x.m_fieldname);
    std::swap(m_tol, x.m_tol);
    BaseManipulation::swap(x);
}

/*!
 * Build the ports of the class;
 */
void
SelectField::buildPorts(){
    bool built = true;
    built = (built && createPortIn<mimmo::MimmoSharedPointer<MimmoObject>, SelectField>(this, &mimmo::SelectField::setGeometry, M_GEOM));
    built = (built && createPortIn<std::string, SelectField>(this, &mimmo::SelectField::setFieldName, M_NAME));
    built = (built && createPortOut<std::string, SelectField>(this, &mimmo::SelectField::getFieldName, M_NAME));
    built = (built && createPortOut<mimmo::MimmoSharedPointer<MimmoObject>, SelectField>(this, &mimmo::BaseManipulation::getGeometry, M_GEOM));
    m_arePortsBuilt = built;
}

/*!
 * Set the name of your target Selected field.
 * Note. This is the name used as result field after selection; in case of SelectType by geoemtry,
 * the name of the field is replaced by the selected one.
 * \param[in] fieldname  Input field name
 */
void
SelectField::setFieldName(std::string fieldname){
	m_fieldname = fieldname;
};

/*!
 * Set selection method.
 * \param[in] mode  Selection mode
 */
void
SelectField::setMode(SelectType mode){
    m_mode = mode;
};

/*!
 * Set selection method.
 * \param[in] mode  Selection mode
 */
void
SelectField::setMode(int mode){
	setMode(static_cast<SelectType>(mode));
};

/*!
 * Set tolerance for extraction by patch.
 * \param[in] tol Tolerance
 */
void
SelectField::setTolerance(double tol){
    m_tol = std::max(1.0e-12, tol);
};

/*!
 * Get the name of your target Selected field.
 * \return Target field name
 */
std::string
SelectField::getFieldName(){
	return m_fieldname;
};
/*!
 * Clear all stuffs in your class
 */
void
SelectField::clear(){
    BaseManipulation::clear();
};

/*!Execution command.
 * Select a field over a target geometry from a vector of fields.
 */
void
SelectField::execute(){

	bool check = mSelect();
    if(!check){
        (*m_log)<<"Error in "<<m_name<<". Field cannot be Selected"<<std::endl;
        (*m_log)<<"This could be due to not correct settings:"<<std::endl;
        switch(m_mode){
            case SelectType::GEOMETRY :
                (*m_log)<<"    missing or unfound ref Geometry to select field ->see method setGeometry//M_GEOM port"<<std::endl;
                break;
            case SelectType::MAPPING :
                (*m_log)<<"    empty mapping performed or missing mapping Geometry to select field ->see method setGeometry//M_GEOM port"<<std::endl;
                break;
            case SelectType::NAME :
                (*m_log)<<"    missing or unfound name to select field ->see method setFieldName//M_NAME port"<<std::endl;
                break;
            default:
                //never been reached
                break;
        }
        (*m_log)<<"Return empty field"<<std::endl;
    }
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectField::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    BaseManipulation::absorbSectionXML(slotXML, name);


    if(slotXML.hasOption("SelectType")){
        std::string input = slotXML.get("SelectType");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setMode(value);
    }

    if(slotXML.hasOption("FieldName")){
        std::string input = slotXML.get("FieldName");
        input = bitpit::utils::string::trim(input);
        if(!input.empty()){
            setFieldName(input);
        }
    }
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
void SelectField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("SelectType", std::to_string(static_cast<int>(m_mode)));
    if (m_mode == SelectType::MAPPING)
    	slotXML.set("Tolerance", std::to_string(m_tol));

    slotXML.set("FieldName", m_fieldname);

};

}
