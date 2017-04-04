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
//SELECTION	BY BOX class 	**********************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByBox::SelectionByBox(){
	m_name = "mimmo.SelectionByBox";
	m_type = SelectionType::BOX;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByBox::SelectionByBox(const bitpit::Config::Section & rootXML){
	
	m_name = "mimmo.SelectionByBox";
	m_type = SelectionType::BOX;

	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "mimmo.SelectionByBox"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::SelectionByBox constructor. No valid xml data found"<<std::endl;
	};
}

/*!
 * Custom Constructor
 * \param[in] origin of the box->baricenter
 * \param[in] span	 of the box, width/height/depth
 * \param[in] target	pointer to a target geometry
 */
SelectionByBox::SelectionByBox(darray3E origin, darray3E span, MimmoObject * target){
	m_name = "mimmo.SelectionByBox";
	m_type = SelectionType::BOX;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
};

/*!
 * Destructor
 */
SelectionByBox::~SelectionByBox(){};

/*!
 * Copy Constructor
 */
SelectionByBox::SelectionByBox(const SelectionByBox & other):GenericSelection(), Cube(){
	*this = other;
};

/*!
 * Copy operator
 */
SelectionByBox & SelectionByBox::operator=(const SelectionByBox & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection *>(&other));
	*(static_cast<Cube * >(this)) = *(static_cast<const Cube *>(&other));
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByBox::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByBox>(this, &SelectionByBox::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, SelectionByBox>(this, &SelectionByBox::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));

	m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByBox::clear(){
	m_subpatch.reset(nullptr);
	BaseManipulation::clear();
};



/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByBox::extractSelection(){
	switch(m_topo){
		case 3:
			if(m_dual)	return  excludeCloudPoints(getGeometry());
			else		return	includeCloudPoints(getGeometry());
			break;
		default:
			if(m_dual)	return  excludeGeometry(getGeometry());
			else		return	includeGeometry(getGeometry());
			break;
	}
	
	return livector1D(0); //never been reached
};


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
void SelectionByBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
	//start absorbing
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
 * - <B>ClassName</B>: name of the class as "mimmo.SelectionByBox"
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
void SelectionByBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
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