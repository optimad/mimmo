/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

//------------------------------------------------------------------------
//SELECTION BY BOX WITH SCALAR class    **********************************
//------------------------------------------------------------------------

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
	input = bitpit::utils::trim(input);
	if(input == "mimmo.SelectionByBoxWithScalar"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::SelectionByBoxWithScalar constructor. No valid xml data found"<<std::endl;
	};
}

/*!
 * Custom Constructor
 * \param[in] origin of the box->baricenter
 * \param[in] span   of the box, width/height/depth
 * \param[in] target    pointer to a target geometry
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
SelectionByBoxWithScalar::SelectionByBoxWithScalar(const SelectionByBoxWithScalar & other):SelectionByBox(){
    *this = other;
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
void SelectionByBoxWithScalar::buildPorts(){

    bool built = true;

    //inheritance
    SelectionByBox::buildPorts();

    //input
    built = (built && createPortIn<dvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

    //output
    built = (built && createPortOut<dvector1D, SelectionByBoxWithScalar>(this, &SelectionByBoxWithScalar::getField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));

    m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByBoxWithScalar::clear(){
    m_field.clear();
    SelectionByBox::clear();
};

/*! It sets the starting scalar field attached to the whole patch.
 * \param[in] field Scalar field.
 */
void
SelectionByBoxWithScalar::setField(dvector1D field){
    m_field = field;
}

/*! It gets the scalar field attached to the extracted patch.
 * \return Scalar field.
 */
dvector1D
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
void SelectionByBoxWithScalar::execute(){

    SelectionByBox::execute();

    if (m_field.size() != 0){
        m_field.resize(getGeometry()->getNVertex(), 0.0);
        bitpit::PiercedVector<bitpit::Vertex> vertices = m_subpatch->getVertices();
        dvector1D field_tmp(vertices.size());
        for (auto vertex : vertices){
            field_tmp[m_subpatch->getMapDataInv(vertex.getId())] = m_field[getGeometry()->getMapDataInv(vertex.getId())];
        }
        m_field = field_tmp;
    }
}


/*!
 * Plot optional result of the class in execution, that is the selected patch
 * as standard vtk unstructured grid and the related scalar field.
 */
void SelectionByBoxWithScalar::plotOptionalResults(){
    if(getPatch()->isEmpty()) return;
    std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";

    getPatch()->getPatch()->getVTK().addData("field", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, m_field);

    getPatch()->getPatch()->write(name);
}


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 *  --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual
 * - <B>Origin</B>: array of 3 doubles identifying origin
 * - <B>Span</B>: span of the box (width, height, depth)
 * - <B>RefSystem</B>: reference system of the box;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByBoxWithScalar::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	//start absorbing
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
	
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}
	
	if(slotXML.hasOption("Origin")){
		std::string input = slotXML.get("Origin");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setOrigin(temp);
		}else{
			setOrigin(temp);
		}	
	}
	
	if(slotXML.hasOption("Span")){
		std::string input = slotXML.get("Span");
		input = bitpit::utils::trim(input);
		darray3E temp = {{1.0,1.0,1.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setSpan(temp);
		}else{
			setSpan(temp);
		}	
	}	
	
	if(slotXML.hasSection("RefSystem")){
		
		const bitpit::Config::Section & axesXML = slotXML.getSection("RefSystem");
		dmatrix33E axes;
		for(int i=0; i<3; ++i){
			axes[i].fill(0.0);
			axes[i][i] = 1.0;
		}
		
		if(axesXML.hasOption("axis0")){
			std::string input = axesXML.get("axis0");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
			}
		}
		
		if(axesXML.hasOption("axis1")){
			std::string input = axesXML.get("axis1");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
			}
		}
		
		if(axesXML.hasOption("axis2")){
			std::string input = axesXML.get("axis2");
			input = bitpit::utils::trim(input);
			if(!input.empty()){
				std::stringstream ss(input);
				ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
			}
		}
		
		setRefSystem(axes);
	}	
	
	if(slotXML.hasOption("PlotInExecution")){
		std::string input = slotXML.get("PlotInExecution");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setPlotInExecution(value);
	}
	
	if(slotXML.hasOption("OutputPlot")){
		std::string input = slotXML.get("OutputPlot");
		input = bitpit::utils::trim(input);
		std::string temp = ".";
		if(!input.empty())	setOutputPlot(input);
		else			  	setOutputPlot(temp);
	}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 *  --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.SelectionByBoxWithScalar"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual
 * - <B>Origin</B>: array of 3 doubles identifying origin
 * - <B>Span</B>: span of the box (width, height, depth)
 * - <B>RefSystem</B>: reference system of the box;
 * 					<RefSystem>
 * 						<axis0>	1.0 0.0 0.0 </axis0>
 * 						<axis1>	0.0 1.0 0.0 </axis1>
 * 						<axis2>	0.0 0.0 1.0 </axis2>
 * 					</RefSystem>
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByBoxWithScalar::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	{
		darray3E org = getOrigin();
		std::stringstream ss;
		ss<<std::scientific<<org[0]<<'\t'<<org[1]<<'\t'<<org[2];
		slotXML.set("Origin",ss.str());
	}
	
	{
		darray3E span = getSpan();
		std::stringstream ss;
		ss<<std::scientific<<span[0]<<'\t'<<span[1]<<'\t'<<span[2];
		slotXML.set("Span",ss.str());
	}
	
	{
		dmatrix33E axes = getRefSystem();
		bitpit::Config::Section & axesXML = slotXML.addSection("RefSystem");
		
		for(int i=0; i<3; ++i){
			std::string name = "axis"+std::to_string(i);
			std::stringstream ss;
			ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
			axesXML.set(name, ss.str());
		}
	}
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
	return;
};

}

