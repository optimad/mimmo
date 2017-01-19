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
 * BaseManipulation::absorbSectionXML. Except of geometry parameters and filed (which is instantiated internally
 * or passed by port linking), the class reads no other siginificant parameters.
 *  
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void SplitField::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	std::string input; 
	int counter;
	
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



//SPLITSCALARFIELD IMPLEMENTATION//

/*!
 * Default constructor. Requires topo flag. 1-surface, 2-volume, 3-pointcloud.
 */
SplitScalarField::SplitScalarField(int topo):SplitField(topo){
	m_name = "MiMMO.SplitScalarField";
}


/*!
 * Default destructor
 */
SplitScalarField::~SplitScalarField(){
	m_field.clear();
	m_result.clear();
}

/*!
 * Build the ports of the class;
 */
void
SplitScalarField::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvector1D, SplitScalarField>(this, &mimmo::SplitScalarField::setField, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<std::unordered_map<MimmoObject*, dvector1D* >, SplitScalarField>(this, &mimmo::SplitScalarField::getSplittedData, PortType::M_UMGEOSFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECFLOAT_));	
	
	SplitField::buildPorts();
	m_arePortsBuilt = built;
}

/*!
 * Get splitted field as a map having pointers to splitted original geometries 
 * and the portions of field relative to them. 
 */
std::unordered_map<MimmoObject*, dvector1D *> 
SplitScalarField::getSplittedData(){
	
	std::unordered_map<MimmoObject*, dvector1D *> map;
	if(m_originals.empty() || m_geometry==NULL) return map;
	if(m_result.size() != m_originals.size())	return map;
		
	int counter = 0;
	for(auto & val : m_originals){
		map[val] = &(m_result[counter]);
		++counter;
	}
	
	return map;
}

/*!
 * Set Field associated to the target geometry and that need to splitted.
 * If the field is associated to the cells or to points of the target geometry,
 * please set this info, choosing the correct division map between setCellDivisionMap or 
 * setVertDivisionMap methods.  
 * \param[in]	field scalar field fo doubles 
 */
void
SplitScalarField::setField(dvector1D field){
	m_field = field;
}

/*!
 * Clear content of the class
 */
void
SplitScalarField::clear(){
	m_field.clear();
	m_result.clear();
	SplitField::clear();
}

/*!
 * Plot splitted field alongside its geometries ;
 */

void 
SplitScalarField::plotOptionalResults(){
	if(isEmpty()) return;
	if(m_originals.size() != m_result.size())	return;
	
	bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
	if(m_mapVertDivision.empty())	loc = bitpit::VTKLocation::CELL;
	
	int counter=0;
	for(auto & geo : m_originals){
		dvecarr3E points = geo->getVertexCoords();
		
		if(getTopo() != 3){
			ivector2D connectivity = geo->getCompactConnectivity();
			bitpit::VTKElementType cellType = desumeElement(connectivity);
			bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),cellType);
			output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
			output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
			output.setDimensions(connectivity.size(), points.size());
			output.addData("field", bitpit::VTKFieldType::SCALAR, loc, m_result[counter]);
			output.setCounter(counter);
			output.setCodex(bitpit::VTKFormat::APPENDED);
			output.write();
		}else{
			int size = points.size();
			ivector1D connectivity(size);
			for(int i=0; i<size; ++i)	connectivity[i]=i;
			bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),	bitpit::VTKElementType::VERTEX);
			output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
			output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
			output.setDimensions(connectivity.size(), points.size());
			output.addData("field", bitpit::VTKFieldType::SCALAR, loc, m_result[counter]);
			output.setCounter(counter);
			output.setCodex(bitpit::VTKFormat::APPENDED);
			output.write();
		}

		++counter;
	}
}

/*!
 * Split your original field along the splitted original geometries provided
 */
bool
SplitScalarField::split(){
	if(m_mapCellDivision.empty() && m_mapVertDivision.empty())	return false;
	if(isEmpty())	return false;
	
	//check original field;
	bool loc = m_mapVertDivision.empty(); //true by cells, false by vert.
	
	int nGeo = m_originals.size();
	m_result.resize(nGeo);
	
	std::unordered_map<long, std::pair<int, long>> * mapp;
	livector1D * mapIdTarget;
	std::vector< std::map<long, int> * > invIdLoc(nGeo);
	
	if(loc)	{
		m_field.resize(getGeometry()->getNCells(), 0.0);
		mapp = &m_mapCellDivision;	
		mapIdTarget = &(getGeometry()->getMapCell());
		for(int i=0; i<nGeo; ++i){
			invIdLoc[i] = &(m_originals[i]->getMapCellInv());
		}
	}
	else	{
		m_field.resize(getGeometry()->getNVertex(), 0.0);
		mapp = &m_mapVertDivision;	
		mapIdTarget = &(getGeometry()->getMapData());
		for(int i=0; i<nGeo; ++i){
			invIdLoc[i] = &(m_originals[i]->getMapDataInv());
		}
	}	
	
	if(m_field.size() != mapp->size())	return false;
	
	//allocate memory for results;	
	for(int i=0; i<nGeo; ++i){
		
		if(loc)	m_result[i].resize(m_originals[i]->getNCells(),0.0);
		else	m_result[i].resize(m_originals[i]->getNVertex(),0.0);
	}
	
	int counter = 0;
	long IDtarget;
	std::pair<int,long> douple;
	int locali;
	for(auto & val: m_field){
		
		IDtarget = (*mapIdTarget)[counter];
		if(mapp->count(IDtarget)> 0 ){
			douple = (*mapp)[IDtarget]; 
			if(douple.first > nGeo)	continue;
			if((invIdLoc[douple.first])->count(douple.second) > 0){
				locali = (*invIdLoc[douple.first])[douple.second];
				m_result[douple.first][locali] = val;
			}	
		}
		++counter;
	}
	return true;
}

//SPLITVECTORFIELD IMPLEMENTATION//

/*!
 * Default constructor. Requires topo flag. 1-surface, 2-volume, 3-pointcloud.
 */
SplitVectorField::SplitVectorField(int topo):SplitField(topo){
	m_name = "MiMMO.SplitVectorField";
}


/*!
 * Default destructor
 */
SplitVectorField::~SplitVectorField(){
	m_field.clear();
	m_result.clear();
}

/*!
 * Build the ports of the class;
 */
void
SplitVectorField::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, SplitVectorField>(this, &mimmo::SplitVectorField::setField, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<std::unordered_map<MimmoObject*, dvecarr3E* >, SplitVectorField>(this, &mimmo::SplitVectorField::getSplittedData, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));	
	SplitField::buildPorts();
	m_arePortsBuilt = built;
}

/*!
 * Get splitted field as a map having pointers to splitted original geometries 
 * and the portions of field relative to them. 
 */
std::unordered_map<MimmoObject*, dvecarr3E *> 
SplitVectorField::getSplittedData(){
	
	std::unordered_map<MimmoObject*, dvecarr3E *> map;
	if(m_originals.empty() || m_geometry==NULL) return map;
	if(m_result.size() != m_originals.size())	return map;

	int counter = 0;
	for(auto & val : m_originals){
		map[val] = &(m_result[counter]);
		++counter;
	}

	return map;
}

/*!
 * Set Field associated to the target geometry and that need to splitted.
 * If the field is associated to the cells or to points of the target geometry,
 * please set this info, choosing the correct division map between setCellDivisionMap or 
 * setVertDivisionMap methods.  
 * \param[in]	field vector field of array at 3 double elements 
 */
void
SplitVectorField::setField(dvecarr3E field){
	m_field = field;
}

/*!
 * Clear content of the class
 */
void
SplitVectorField::clear(){
	m_field.clear();
	m_result.clear();
	SplitField::clear();
}

/*!
 * Plot splitted field alongside its geometries ;
 */

void 
SplitVectorField::plotOptionalResults(){
	if(isEmpty()) return;
	if(m_originals.size() != m_result.size())	return;
	
	bitpit::VTKLocation loc = bitpit::VTKLocation::POINT;
	if(m_mapVertDivision.empty())	loc = bitpit::VTKLocation::CELL;
	
	int counter=0;
	for(auto & geo : m_originals){
		dvecarr3E points = geo->getVertexCoords();
		
		if(getTopo() != 3){
			ivector2D connectivity = geo->getCompactConnectivity();
			bitpit::VTKElementType cellType = desumeElement(connectivity);
			bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),cellType);
			output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
			output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
			output.setDimensions(connectivity.size(), points.size());
			output.addData("field", bitpit::VTKFieldType::VECTOR, loc, m_result[counter]);
			output.setCounter(counter);
			output.setCodex(bitpit::VTKFormat::APPENDED);
			output.write();
		}else{
			int size = points.size();
			ivector1D connectivity(size);
			for(int i=0; i<size; ++i)	connectivity[i]=i;
			bitpit::VTKUnstructuredGrid output(".",m_name+std::to_string(getClassCounter()),	bitpit::VTKElementType::VERTEX);
			output.setGeomData( bitpit::VTKUnstructuredField::POINTS, points);
			output.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, connectivity);
			output.setDimensions(connectivity.size(), points.size());
			output.addData("field", bitpit::VTKFieldType::VECTOR, loc, m_result[counter]);
			output.setCounter(counter);
			output.setCodex(bitpit::VTKFormat::APPENDED);
			output.write();
		}
		
		++counter;
	}
}

/*!
 * Split your original field along the splitted original geometries provided
 */
bool
SplitVectorField::split(){
	if(m_mapCellDivision.empty() && m_mapVertDivision.empty())	return false;
	if(isEmpty())	return false;
	
	//check original field;
	bool loc = m_mapVertDivision.empty(); //true by cells, false by vert.
	
	int nGeo = m_originals.size();
	m_result.resize(nGeo);
	
	std::unordered_map<long, std::pair<int, long>> * mapp;
	livector1D * mapIdTarget;
	std::vector< std::map<long, int> * > invIdLoc(nGeo);
	
	if(loc)	{
		m_field.resize(getGeometry()->getNCells(), {{0.0,0.0,0.0}});
		mapp = &m_mapCellDivision;	
		mapIdTarget = &(getGeometry()->getMapCell());
		for(int i=0; i<nGeo; ++i){
			invIdLoc[i] = &(m_originals[i]->getMapCellInv());
		}
	}
	else	{
		m_field.resize(getGeometry()->getNVertex(), {{0.0,0.0,0.0}});
		mapp = &m_mapVertDivision;	
		mapIdTarget = &(getGeometry()->getMapData());
		for(int i=0; i<nGeo; ++i){
			invIdLoc[i] = &(m_originals[i]->getMapDataInv());
		}
	}	
	
	if(m_field.size() != mapp->size())	return false;
	
	//allocate memory for results;	
	for(int i=0; i<nGeo; ++i){
		
		if(loc)	m_result[i].resize(m_originals[i]->getNCells(),{{0.0,0.0,0.0}});
		else	m_result[i].resize(m_originals[i]->getNVertex(),{{0.0,0.0,0.0}});
	}
	
	int counter = 0;
	long IDtarget;
	std::pair<int,long> douple;
	int locali;
	for(auto & val: m_field){
		
		IDtarget = (*mapIdTarget)[counter];
		if(mapp->count(IDtarget)> 0 ){
			douple = (*mapp)[IDtarget]; 
			if(douple.first > nGeo)	continue;
			if((invIdLoc[douple.first])->count(douple.second) > 0){
				locali = (*invIdLoc[douple.first])[douple.second];
				m_result[douple.first][locali] = val;
			}	
		}
		++counter;
	}
	
	return true;
}
