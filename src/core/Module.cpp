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

#include "Module.hpp"

namespace mimmo{

/*!Default constructor of Module.
 */
Module::Module(){
    m_name = "mimmo.Module";
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Module::Module(const bitpit::Config::Section & rootXML){

    std::string fallback_name = "ClassNONE";
    std::string input_name = rootXML.get("ClassName", fallback_name);
    input_name = bitpit::utils::string::trim(input_name);

    m_name = "mimmo.Module";

    if(input_name == "mimmo.Module"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*!
 * Default destructor of Module.
 */
Module::~Module(){};

/*!
 * Copy Constructor.
 * \param[in] other class of type Module
 */
Module::Module(const Module & other):BaseManipulation(other){
	m_field = other.m_field;
	m_result = other.m_result;
}

/*!
 * Swap function of Module
 * \param[in] x object to be swapped.
 */
void Module::swap(Module & x ) noexcept
{
    m_field.swap(x.m_field);
    m_result.swap(x.m_result);
    BaseManipulation::swap(x);
};

/*!
 * Assignement operator of Module.
 * \param[in] other class of type Module
 */
Module & Module::operator=(Module other){
    this->swap(other);
	return *this;
};

/*!
 * Build the ports of the class;
 */
void
Module::buildPorts(){
    bool built = true;
    built = (built && createPortIn<dmpvecarr3E*, Module>(this, &mimmo::Module::setField, M_VECTORFIELD));
    built = (built && createPortOut<dmpvector1D*, Module>(this, &mimmo::Module::getResult, M_SCALARFIELD));
    m_arePortsBuilt = built;
}

/*!
 * Clear all stuffs in your class
 */
void
Module::clear(){
    m_field.clear();
    m_result.clear();
    BaseManipulation::clear();
};

/*!
 * Set input field for computing.
 * \param[in] field input field
 */
void
Module::setField(dmpvecarr3E *field){
    if(!field)  return;
    m_field = *field;
}

/*!
 * Get resulting scalar field.
 * \return magnitude field
 */
dmpvector1D *
Module::getResult(){
    return &m_result;
}

/*!Execution command.
 * Compute the module field of an input vector field.
 */
void
Module::execute(){

    mimmo::MPVLocation loc = m_field.getDataLocation();

    m_result.clear();
    m_result.setDataLocation(loc);
    m_result.setGeometry(m_field.getGeometry());

    for (long id : m_field.getIds()){
    	double val = norm2(m_field[id]);
    	m_result.insert(id, val);
    }
}


/*!
 * Plot computed field over the linked geometry
 */
void
Module::plotOptionalResults(){

	m_result.setName("magnitude");
	write(m_result.getGeometry(), m_result);

}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Module::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void Module::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

};

}
