/*---------------------------------------------------------------------------*\
 * 
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
namespace mimmo{


//------------------------------------------------------------------------
//SELECTION	BY CYLINDER class 	**********************************************
//------------------------------------------------------------------------

/*!
 * Basic Constructor
 */
SelectionByCylinder::SelectionByCylinder(){
	m_name = "MiMMO.SelectionByCylinder";
	m_type = SelectionType::CYLINDER;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByCylinder::SelectionByCylinder(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.SelectionByCylinder";
	m_type = SelectionType::CYLINDER;

	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.SelectionByCylinder"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::SelectionByCylinder constructor. No valid xml data found"<<std::endl;
	};
}

/*!
 * Custom Constructor. Pay attention span of angular coordinate must be at most 2*pi.
 * \param[in] origin of the cylinder->baricenter
 * \param[in] span	 of the cylinder, basis radius/ span of angular coord in radians/height
 * \param[in] infLimTheta	starting origin of the angular coordinate. default is 0 radians.
 * \param[in] mainAxis	orientation of the cylinder height axis 
 * \param[in] target	pointer to a target geometry 
 */
SelectionByCylinder::SelectionByCylinder(darray3E origin, darray3E span, double infLimTheta, darray3E mainAxis, MimmoObject * target){
	m_name = "MiMMO.SelectionByCylinder";
	m_type = SelectionType::CYLINDER;
	setGeometry(target);
	setOrigin(origin);
	setSpan(span[0],span[1],span[2]);
	setInfLimits(infLimTheta,1);
	setRefSystem(2, mainAxis);
};

/*!
 * Destructor
 */
SelectionByCylinder::~SelectionByCylinder(){};

/*!
 * Copy Constructor
 */
SelectionByCylinder::SelectionByCylinder(const SelectionByCylinder & other):GenericSelection(), Cylinder(){
	*this = other;
};

/*!
 * Copy operator
 */
SelectionByCylinder & SelectionByCylinder::operator=(const SelectionByCylinder & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection *>(&other));
	*(static_cast<Cylinder * >(this)) = *(static_cast<const Cylinder *>(&other));
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByCylinder::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	//input
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dmatrix33E, SelectionByCylinder>(this, &SelectionByCylinder::setRefSystem, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, SelectionByCylinder>(this, &SelectionByCylinder::setInfLimits, PortType::M_INFLIMITS, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));

	m_arePortsBuilt = built;
};

/*!
 * Clear content of the class
 */
void SelectionByCylinder::clear(){
	m_subpatch.reset(nullptr);
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByCylinder::extractSelection(){
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
 * - <B>Span</B>:span of the cylinder (base radius, azimuthal span(in radians), height)
 * - <B>RefSystem</B>: reference system of the cylinder (z/2 is the cylinder height axis);
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * - <B>InfLimits</B>: inferior limits for span dimensioning
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByCylinder::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
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
		darray3E temp = {{1.0,2.0*M_PI,1.0}};
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
	
	if(slotXML.hasOption("InfLimits")){
		std::string input = slotXML.get("InfLimits");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp[0]>>temp[1]>>temp[2];
			setInfLimits(temp);
		}else{
			setInfLimits(temp);
		}	
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
 *   --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "MiMMO.SelectionByCylinder"
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual
 * - <B>Origin</B>: array of 3 doubles identifying origin
 * - <B>Span</B>:span of the cylinder (base radius, azimuthal span(in radians), height)
 * - <B>RefSystem</B>: reference system of the cylinder (z/2 is the cylinder height axis);
 * 					<RefSystem>
 * 						<axis0>	1.0 0.0 0.0 </axis0>
 * 						<axis1>	0.0 1.0 0.0 </axis1>
 * 						<axis2>	0.0 0.0 1.0 </axis2>
 * 					</RefSystem> 
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * - <B>InfLimits</B>: inferior limits for span dimensioning
 * 
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByCylinder::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
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
	
	{
		darray3E inflim = getInfLimits();
		std::stringstream ss;
		ss<<std::scientific<<inflim[0]<<'\t'<<inflim[1]<<'\t'<<inflim[2];
		slotXML.set("InfLimits",ss.str());
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
