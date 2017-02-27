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

#include "SplitFields.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


/*!Default constructor of SplitField.
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 */
SplitField::SplitField(int topo){
	m_topo = topo;
}

/*!
 * Default destructor of SplitField.
 */
SplitField::~SplitField(){
	clear();
};

/*!Copy constructor of SplitField.Soft Copy of MimmoObject;
 */
SplitField::SplitField(const SplitField & other){
	*this = other;
};	

/*!
 * Assignement operator of SplitField. Soft copy of MimmoObject
 */
SplitField & SplitField::operator=(const SplitField & other){
	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_topo = other.m_topo;
	m_originals = other.m_originals;
	m_mapCellDivision = other.m_mapCellDivision;
	m_mapVertDivision = other.m_mapVertDivision;

	return *this;
};

/*!
 * Build the ports of the class;
 */
void
SplitField::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, SplitField>(this, &mimmo::SplitField::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<std::vector<MimmoObject*>, SplitField>(this, &mimmo::SplitField::setSplittedGeometries, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));

	built = (built && createPortIn<std::unordered_map<long,std::pair<int, long> >, SplitField>(this, &mimmo::SplitField::setCellDivisionMap, PortType::M_MAPDCELL, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));
	built = (built && createPortIn<std::unordered_map<long,std::pair<int, long> >, SplitField>(this, &mimmo::SplitField::setVertDivisionMap, PortType::M_MAPDVERT, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGPAIRINTLONG));	
	m_arePortsBuilt = built;
}


/*!
 * Set the target geometry where your field is defined.Topology of the geometry must be coeherent
 * with topology of the class; Reimplementation of BaseManipulation::setGeometry
 * \param[in] geo  pointer to MimmoObject
 */
void
SplitField::setGeometry(MimmoObject* geo){
		if(geo->isEmpty()) return;
		if(geo->getType() != m_topo)	return;
		
		m_geometry = geo;
};

/*!
 * Set original splitted geometries which refers to target geometry. List will be saved as is, replacing any other saved list. 
 */
void
SplitField::setSplittedGeometries(std::vector<MimmoObject *> originals){
	
	if(originals.empty()) return;
	m_originals.resize(originals.size());
	
	int counter = 0;
	for(auto & obj: originals){
		if(obj->isEmpty()) continue;
		if(obj->getType() == m_topo){
			m_originals[counter]= obj;
			++counter;
		}
	}
	m_originals.resize(counter);	
};

/*!
 * It sets the Cell division map relative to the target geometry w.r.t the original splitted geeometries
 * The class interprets this command as an explicit declaration that your current field to be splitted is
 * referred to geometry-cells/simplicies. If any previous compiled Vertex Division map is set, it will be erased.
 * If clas topology is a point cloud, the method do nothing. Please use setVertDivisionMap instead;
 * \param[in] map cell division map of the target geometry w.r.t the original splitted geometries
 */
void
SplitField::setCellDivisionMap(std::unordered_map<long, std::pair<int,long > > map){
	if(m_topo == 3) return;
	m_mapCellDivision = map;
	m_mapVertDivision.clear();
}

/*!
 * It sets the Vertex division map relative to the target geometry w.r.t the original splitted geeometries
 * The class interprets this command as an explicit declaration that your current field to be splitted is
 * referred to geometry-vertices/points. If any previous compiled Cell Division map is set, it will be erased.
 * \param[in] map vertex division map of the target geometry w.r.t the original splitted geometries
 */
void
SplitField::setVertDivisionMap(std::unordered_map<long, std::pair<int,long > > map){
	
	m_mapVertDivision = map;
	m_mapCellDivision.clear();
}
/*!
 * Check if target geometry and its splitted originals are present or not.
 * True - no geometry present, False otherwise.
 */
bool 
SplitField::isEmpty(){
	return (m_geometry == NULL || m_originals.empty());
}

/*!
 * Return current topology type set for your class geometries.
 *  1-surface, 2-volume, 3-pointcloud
 */
int 
SplitField::getTopo(){
	return m_topo;
}

/*!
 * Clear all stuffs in your class
 */
void
SplitField::clear(){
	m_originals.clear();
	m_mapCellDivision.clear();
	m_mapVertDivision.clear();
	BaseManipulation::clear();
};

/*!Execution command.
 * stitch together multiple geoemetry in the same object. 
 */
void
SplitField::execute(){

	bool check = split();
	if(!check){
		std::cout<<"Error in class "<<m_name<<". Field cannot be splitted"<<std::endl;
		std::cout<<"This could be due to not correct setting of geometries or division maps"<<std::endl;
	}
}

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameters and field (which is instantiated internally
 * or passed by port linking), the class reads no other siginificant parameters.
 * 
 * * --> Absorbing data:
 * 		Topology: info on admissible topology format 1-surface, 2-volume, 3-pointcloud
 * 		PlotInExecution : boolean 0/1 print optional results of the class.
 * 		OutputPlot : target directory for optional results writing. 
 * 
 *  
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void SplitField::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	std::string input; 
	
	//checking topology
	if(slotXML.hasOption("Topology")){
		std::string input = slotXML.get("Topology");
		input = bitpit::utils::trim(input);
		int temptop = -1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temptop;
		}
		if(m_topo != temptop)	return;
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
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML. Except of geometry parameters and fields (which are instantiated internally
 * or passed by port linking), the class writes no siginificant parameters:
 * 
 * --> Flushing data// how to write it on XML:
 * 		ClassName : name of the class as "MiMMO.Split<Scalar/Vector>Fields"
 * 		ClassID	  : integer identifier of the class	
 * 		Topology: info on admissible topology format 1-surface, 2-volume, 3-pointcloud
 * 		PlotInExecution : boolean 0/1 print optional results of the class.
 * 		OutputPlot : target directory for optional results writing. 
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot 
 */
void SplitField::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	slotXML.set("Topology", m_topo);

	std::string output;
	
	if(isPlotInExecution()){
		slotXML.set("PlotInExecution", std::to_string(1));
	}
	
	if(m_outputPlot != "."){
		slotXML.set("OutputPlot", m_outputPlot);
	}
return;	
};	

/*!
 * Desume Element type from passed typeGeom and connectivity. Return undefined type for unexistent 
 * or unsupported element, or mixed element type connectivity. NEED TO BE MOVED IN MimmoObject
 */
bitpit::VTKElementType	
SplitField::desumeElement(ivector2D & conn){
	bitpit::VTKElementType result = bitpit::VTKElementType::UNDEFINED;
	if(conn.empty() && m_topo != 3)	return	result;
	
	switch(m_topo){
		case	1:
			if(conn[0].size() == 3)		result = bitpit::VTKElementType::TRIANGLE;
			if(conn[0].size() == 4)		result = bitpit::VTKElementType::QUAD;
			break;
		case	2:
			if(conn[0].size() == 4)		result = bitpit::VTKElementType::TETRA;
			if(conn[0].size() == 8)		result = bitpit::VTKElementType::HEXAHEDRON;
			break;
		case	3:
			result=  bitpit::VTKElementType::VERTEX;
			break;
		default : //never been reached
			break;
	}
	
	return result;
};


