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
#include "ProjectCloud.hpp"

using namespace mimmo;

/*!Default constructor of ProjectCloud
*/
ProjectCloud::ProjectCloud(){
	m_name = "MiMMO.ProjectCloud";
};

/*!Default destructor of ProjectCloud
 */
ProjectCloud::~ProjectCloud(){};

/*!Copy constructor of ProjectCloud.
 */
ProjectCloud::ProjectCloud(const ProjectCloud & other){
	*this = other;
};

/*!Assignement operator of ProjectCloud.
 */
ProjectCloud & ProjectCloud::operator=(const ProjectCloud & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_points = other.m_points;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
ProjectCloud::buildPorts(){
	bool built = true;
	
	built = (built && createPortIn<dvecarr3E, ProjectCloud>(this, &mimmo::ProjectCloud::setCoords, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject * , ProjectCloud>(this, &mimmo::ProjectCloud::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<dvecarr3E, ProjectCloud>(this, &mimmo::ProjectCloud::getCloudResult, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It gets the coordinates of points stored in the object.
 * \return Coordinates of points stored in the object.
 */
dvecarr3E
ProjectCloud::getCoords(){
	return(m_points);
};

/*!It gets the resulting points after class functionality execution
 * \return Projected points stored in the object.
 */
dvecarr3E
ProjectCloud::getCloudResult(){
	return(m_proj);
};

/*!It sets the coordinates of original points to be processed.
 * \param[in] coords Coordinates of points to be used .
 */
void
ProjectCloud::setCoords(dvecarr3E coords){
	m_points = coords;
};

/*!
 * clear contents of the class
 */
void ProjectCloud::clear(){
	BaseManipulation::clear();
	m_points.clear();
	m_proj.clear();
}


/*!Execution command. Project list of points on the target geometry.
 */
void
ProjectCloud::execute(){
	
	if(getGeometry() == NULL || getGeometry()->isEmpty())	return;
	
	if(!getGeometry()->isBvTreeBuilt())	getGeometry()->buildBvTree();
	
	//project points on surface.
	int counter = 0;
	m_proj.resize(m_points.size());
	for(auto &val : m_points){
		 m_proj[counter]= bvTreeUtils::projectPoint(&val, getGeometry()->getBvTree()); 
		++counter;
	}
	return;
};

/*!
 * Plot optional result of the class in execution, that is the projected cloud
 * as standard vtk unstructured grid.
 */
void ProjectCloud::plotOptionalResults(){
	bitpit::VTKFormat codex = bitpit::VTKFormat::APPENDED;
	
	int size = m_proj.size();
	ivector1D conn(size);
	for(int i=0; i<size; i++){
		conn[i] = i;
	}
	std::string dir = "./";
	std::string file = m_name + "_" + std::to_string(getClassCounter());
	
	bitpit::VTKUnstructuredGrid vtk(dir, file, bitpit::VTKElementType::VERTEX);
	vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, m_proj) ;
	vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, conn) ;
	vtk.setDimensions(size, size);
	vtk.setCodex(codex);
	
	vtk.write();
}



/*!
* Get infos from a XML bitpit::Config::section. 
* Coordinates are mandatorely passed through ports
* 
* \param[in] slotXML 	bitpit::Config::Section of XML file
* \param[in] name   name associated to the slot
*/
void ProjectCloud::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name ){
	
	//start absorbing

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
 * Plot infos from a XML bitpit::Config::section.
 * Coordinates are mandatorely passed through ports
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void ProjectCloud::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
	
	return;
};

