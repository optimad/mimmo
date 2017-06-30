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
SelectionByBoxWithScalar::SelectionByBoxWithScalar(){
    m_name = "mimmo.SelectionByBoxWithScalar";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(const bitpit::Config::Section & rootXML){
	
	m_name = "mimmo.SelectionByBoxWithScalar";
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::string::trim(input);
	if(input == "mimmo.SelectionByBoxWithScalar"){
		absorbSectionXML(rootXML);
	}else{	
        warningXML(m_log, m_name);
	};
}

/*!
 * Custom Constructor
 * \param[in] origin Origin of the box->baricenter.
 * \param[in] span     Span of the box, width/height/depth.
 * \param[in] target    Pointer to MimmoObject target geometry.
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(darray3E origin, darray3E span, MimmoObject * target){
	m_name = "mimmo.SelectionByBoxWithScalar";
	m_type = SelectionType::BOX;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
};

/*!
 * Destructor
 */
SelectionByBoxWithScalar::~SelectionByBoxWithScalar(){};

/*!
 * Copy Constructor
 */
SelectionByBoxWithScalar::SelectionByBoxWithScalar(const SelectionByBoxWithScalar & other):SelectionByBox(other){
};

/*!
 * Copy operator
 */
SelectionByBoxWithScalar & SelectionByBoxWithScalar::operator=(const SelectionByBoxWithScalar & other){
    *(static_cast<SelectionByBox * >(this)) = *(static_cast<const SelectionByBox *>(&other));
    return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void
SelectionByBoxWithScalar::buildPorts(){

    bool built = true;

    SelectionByBox::buildPorts();

    built = (built && createPortIn<dmpvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::MPVECTOR, mimmo::pin::dataTAG::FLOAT));

    built = (built && createPortOut<dmpvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::getField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::MPVECTOR, mimmo::pin::dataTAG::FLOAT));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void
SelectionByBoxWithScalar::clear(){
    m_field.clear();
    SelectionByBox::clear();
};

/*! It sets the starting scalar field attached to the whole patch.
 * \param[in] field Scalar field.
 */
void
SelectionByBoxWithScalar::setField(dmpvector1D field){
    m_field = field;
}

/*! It gets the scalar field attached to the extracted patch.
 * \return Scalar field.
 */
dmpvector1D
SelectionByBoxWithScalar::getField(){
    return (m_field);
}


/*!
 * Execute your object. A selection is extracted and trasferred in
 * an indipendent MimmoObject structure pointed by m_subpatch member.
 * The extracted field attached to the selection is built starting from the
 * intial whole scalar field given as input and stored in member
 * m_field (modified after the execution).
 */
void
SelectionByBoxWithScalar::execute(){

    SelectionByBox::execute();

    if (m_field.size() != 0){
        MimmoPiercedVector<double> temp;
        temp.setGeometry(m_subpatch.get());
        temp.setName(m_field.getName());
        bitpit::PiercedVector<bitpit::Vertex> vertices = m_subpatch->getVertices();
        for (const auto & vertex : vertices){
            temp[vertex.getId()] = m_field[vertex.getId()];
        }
        m_field.clear();
        m_field = temp;
    }
}


/*!
 * Plot optional result of the class in execution. It plots the selected patch
 * as standard vtk unstructured grid and the related scalar field.
 */
void
SelectionByBoxWithScalar::plotOptionalResults(){
    if(getPatch()->isEmpty()) return;
    std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";

    dvector1D temp;
    for (const auto & vertex : getPatch()->getVertices()){
        temp.push_back(m_field[vertex.getId()]);
    }

    getPatch()->getPatch()->getVTK().addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, temp);

    getPatch()->getPatch()->write(name);
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByBoxWithScalar::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    SelectionByBox::absorbSectionXML(slotXML, name);
};

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
SelectionByBoxWithScalar::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
    
    BITPIT_UNUSED(name);
    SelectionByBox::flushSectionXML(slotXML, name);
};

}

