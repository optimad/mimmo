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
#include "ExtractFields.hpp"
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
    m_field = other.m_field;
};

/*!
 * Copy operator
 */
SelectionByBoxWithScalar & SelectionByBoxWithScalar::operator=(SelectionByBoxWithScalar other){
    swap(other);
    return *this;
};

/*!
 * Swap function. Assign data of this class to another of the same type and vice-versa.
 * Resulting patch of selection is not swapped, ever.
 * \param[in] x SelectionByBoxWithScalar object
 */
void SelectionByBoxWithScalar::swap(SelectionByBoxWithScalar & x) noexcept
{
    //std::swap(m_field, x.m_field);
    m_field.swap(x.m_field);
    SelectionByBox::swap(x);
    
}

/*!
 * It builds the input/output ports of the object
 */
void
SelectionByBoxWithScalar::buildPorts(){

    bool built = true;

    SelectionByBox::buildPorts();

    built = (built && createPortIn<dmpvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::setField, M_SCALARFIELD));

    built = (built && createPortOut<dmpvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::getField, M_SCALARFIELD));

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

    if(m_field.isEmpty()){
        (*m_log)<<"warning in "<<m_name<<" : empty scalar field found"<<std::endl;
        return;
    }
    //check m_field in input if coherent with the linked geometry
    if(m_field.getGeometry() != getGeometry())  {
        throw std::runtime_error(m_name+" : linked scalar field is not referred to target geometry");
    }
    if(getPatch()->getType()==3 && m_field.getDataLocation()!= MPVLocation::POINT){
        (*m_log)<<"warning in "<<m_name<<" : Attempting to extract a non POINT located field on a Point Cloud target geometry. Do Nothing."<<std::endl;
        return;
    }
    
    ExtractScalarField * extractorField = new ExtractScalarField();
    extractorField->setGeometry(getPatch());
    extractorField->setMode(ExtractMode::ID);
    extractorField->setField(m_field);
    extractorField->execute();
    m_field = extractorField->getExtractedField();
    delete extractorField;
}


/*!
 * Plot optional result of the class in execution. It plots the selected patch
 * as standard vtk unstructured grid and the related scalar field.
 */
void
SelectionByBoxWithScalar::plotOptionalResults(){
    if(getPatch() == NULL)      return;
    if(getPatch()->isEmpty()) return;
    if(m_field.getGeometry() != getPatch()) return;
    
    bitpit::VTKLocation loc = bitpit::VTKLocation::UNDEFINED;
    switch(m_field.getDataLocation()){
        case MPVLocation::POINT :
            loc = bitpit::VTKLocation::POINT;
            break;
        case MPVLocation::CELL :
            loc = bitpit::VTKLocation::CELL;
            break;
        default:
            (*m_log)<<"warning: Undefined Reference Location in plotOptionalResults of "<<m_name<<std::endl;
            (*m_log)<<"Interface or Undefined locations are not supported in VTU writing." <<std::endl;
            break;   
    }
    if(loc == bitpit::VTKLocation::UNDEFINED)  return;
    
    if(getPatch()->getType()==3 && m_field.getDataLocation()!= MPVLocation::POINT){
        (*m_log)<<"warning in "<<m_name<<" : Attempting to plot a non POINT located field on a Point Cloud target geometry. Do Nothing."<<std::endl;
        return;
    }
    
    std::string dir  = m_outputPlot;
    std::string name = m_name + "_Patch."+ std::to_string(getId());
        
    //check size of field and adjust missing values to zero for writing purposes only.
    dmpvector1D field_supp = m_field;
    bool checkField = field_supp.completeMissingData(0.0);
    
    if(getPatch()->getType() != 3){
        if(checkField){
            dvector1D field = field_supp.getDataAsVector();
            getPatch()->getPatch()->getVTK().addData("field", bitpit::VTKFieldType::SCALAR, loc, field);
        }
        std::string totpath = dir +"/"+name;
        getPatch()->getPatch()->write(totpath);
        if(checkField)  getPatch()->getPatch()->getVTK().removeData("field");
    }else{
        liimap mapDataInv;
        dvecarr3E points = getPatch()->getVerticesCoords(&mapDataInv);
        ivector2D connectivity;
        bitpit::VTKElementType cellType = bitpit::VTKElementType::VERTEX;
        
        int np = points.size();
        connectivity.resize(np);
        for (int i=0; i<np; i++){
            connectivity[i].resize(1);
            connectivity[i][0] = i;
            
        }
        bitpit::VTKUnstructuredGrid output(dir,name,cellType);
        output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
        output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
        output.setDimensions(connectivity.size(), points.size());
        output.setCodex(bitpit::VTKFormat::APPENDED);
        if(checkField){
            dvector1D field = field_supp.getDataAsVector();
            output.addData("field", bitpit::VTKFieldType::SCALAR, loc, field);
        }
        output.write();
    }
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

