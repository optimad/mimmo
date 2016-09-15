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
\*---------------------------------------------------------------------------*/
#include "ControlDeformMaxDistance.hpp"
#include <valgrind/callgrind.h>

using namespace mimmo;

/*!Default constructor of ControlDeformMaxDistance
*/
ControlDeformMaxDistance::ControlDeformMaxDistance(){
	m_name = "MiMMO.ControlDeformMaxDistance";
	m_maxDist= 0.0 ;
	m_violation = -1.E18;
};

/*!Default destructor of ControlDeformMaxDistance
 */
ControlDeformMaxDistance::~ControlDeformMaxDistance(){};

/*!Copy constructor of ControlDeformMaxDistance.
 */
ControlDeformMaxDistance::ControlDeformMaxDistance(const ControlDeformMaxDistance & other){
	*this = other;
};

/*!
 * Assignement operator of ControlDeformMaxDistance. Create an exact copy of the class,
 * except for the deformation field referred to the target geometry.
 */
ControlDeformMaxDistance & ControlDeformMaxDistance::operator=(const ControlDeformMaxDistance & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_maxDist = other.m_maxDist;
	m_violation = other.m_violation;
	//deformation field is not copied
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ControlDeformMaxDistance::buildPorts(){
	bool built = true;
	
	built = (built && createPortIn<dvecarr3E, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setDefField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setLimitDistance, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<double, ControlDeformMaxDistance>(this, &mimmo::ControlDeformMaxDistance::getViolation, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*! Return the value of violation of deformed geometry, after class execution. 
 *  If value is positive or at least zero, constraint violation occurs, penetrating or touching at least in one point the 
 *  constraint sourface. Otherwise, returning negative values means that no violation occurs 
 *  delimited by the plane (where plane normal pointing), true the exact opposite.
 * \return violation value
 */
double 
ControlDeformMaxDistance::getViolation(){
	return(m_violation);
};

/*!
 * Set the deformative field associated to each point of the target geometry. 
 * Field resize occurs in execution, if point dimension between field and geoemetry does not match.
 * \param[in]	field of deformation
 */
void
ControlDeformMaxDistance::setDefField(dvecarr3E field){
	m_defField.clear();
	m_defField = field;
};

/*! Set limit distance d of the constraint surface. Must be a positive definite value (>= 0).
 *  Given a target geometry surface (open or closed), its constraint surface is intended 
 *  as the iso-level surface where every points are at distance d from the target surface.
 * \param[in] dist limit distance value
 */
void 
ControlDeformMaxDistance::setLimitDistance(double dist){
	m_maxDist = std::fmax(0.0,dist);
};

/*!Execution command. Calculate violation value and store it in the class member m_violation
 */
void
ControlDeformMaxDistance::execute(){

	MimmoObject * geo = getGeometry();
	if(geo == NULL || geo->isEmpty() || !(geo->isBvTreeSupported())) return;
	
	m_defField.resize(getGeometry()->getNVertex(),darray3E{{0.0,0.0,0.0}});
	dvector1D violationField;
	violationField.reserve(m_defField.size());
	
	if(!(geo->isBvTreeBuilt()))	geo->buildBvTree();
	
	CALLGRIND_START_INSTRUMENTATION;
	dvecarr3E points = geo->getVertexCoords();
	points+= m_defField;
	
	double radius = m_maxDist;
	double rate = 0.1;
	int kmax = 150;
	int kiter;
	int count=0;
	long id;
	double dist;
	for(auto &p : points){
		dist =1.0E+18;
		kiter = 0;
		while(dist == 1.0E+18 && kiter < kmax){
			radius *=(1.0 + rate);
			dist = bvTreeUtils::distance(&p, geo->getBvTree(), id, radius);
			kiter++;
			count++;
		}
		
		if(kiter == kmax)	dist = m_maxDist - dist;
		violationField.push_back(dist);
	}
	
	std::cout<<"called bvTreeUtils::distance "<<count<<" times"<<std::endl;
	
	m_violation = -1.E+18;
	for(auto & val : violationField){
		m_violation = std::fmax(m_violation, (val-m_maxDist));
	}
	
	CALLGRIND_STOP_INSTRUMENTATION;
	CALLGRIND_DUMP_STATS;
	return;
};



/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * 1) LimitDistance	- contraint surface distance from target geometry
 * 
 * Geometry and its deformation fiels are mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ControlDeformMaxDistance::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	//start absorbing
	if(slotXML.hasOption("LimitDistance")){
		std::string input = slotXML.get("LimitDistance");
		input = bitpit::utils::trim(input);
		double value = 0.0;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setLimitDistance(value);
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
 * 1) LimitDistance	- contraint surface distance from target geometry
 * 
 * Geometry and its deformation fiels are mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ControlDeformMaxDistance::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	slotXML.set("LimitDistance", std::to_string(m_maxDist));

	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	return;
};




