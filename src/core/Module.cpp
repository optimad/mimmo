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

    if (m_result.size() == 0 || m_result.getGeometry() == NULL) return;

    bitpit::VTKLocation loc = bitpit::VTKLocation::UNDEFINED;
    switch(m_result.getDataLocation()){
        case MPVLocation::POINT :
            loc = bitpit::VTKLocation::POINT;
        break;
        case MPVLocation::CELL :
            loc = bitpit::VTKLocation::CELL;
        break;
        default:
            (*m_log)<<"warning: Undefined Reference Location in plotOptionalResults of "<<m_name<<std::endl;
            (*m_log)<<"Interface or Undefined locations are not supported in VTU writing." <<std::endl;
            return;
        break;
    }

    //check size of field and adjust missing values to zero for writing purposes only.
    dmpvector1D field_supp = m_result;
    if(!field_supp.completeMissingData(0.0)) return;
    dvector1D field = field_supp.getDataAsVector();


    if(m_result.getGeometry()->getType() != 3){
        m_result.getGeometry()->getPatch()->getVTK().setDirectory(m_outputPlot+"/");
        m_result.getGeometry()->getPatch()->getVTK().setName(m_name+std::to_string(getId()));
        m_result.getGeometry()->getPatch()->getVTK().addData("magnitude", bitpit::VTKFieldType::SCALAR, loc, field);
        m_result.getGeometry()->getPatch()->write();
        m_result.getGeometry()->getPatch()->getVTK().removeData("magnitude");
    }else{
        if(loc == bitpit::VTKLocation::CELL){
            (*m_log)<<"Warning: Attempt writing Cell data field on cloud of points in plotOptionalResults of "<<m_name<<std::endl;
            throw std::runtime_error("Attempt writing Cell data field on cloud of points");
        }

        liimap mapDataInv;
        dvecarr3E points = m_result.getGeometry()->getVerticesCoords(&mapDataInv);
        int size = points.size();
        ivector2D connectivity(size, ivector1D(1));
        for(int i=0; i<size; ++i)    connectivity[i][0]=i;
        bitpit::VTKUnstructuredGrid output(m_outputPlot,m_name+std::to_string(getId()),bitpit::VTKElementType::VERTEX);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
        output.setDimensions(connectivity.size(), points.size());
        output.addData("magnitude", bitpit::VTKFieldType::SCALAR, loc, field);
        output.setCodex(bitpit::VTKFormat::APPENDED);
        output.write();
    }
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
