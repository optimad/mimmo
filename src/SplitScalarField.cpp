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


/*!
 * Default constructor. Requires topo flag. 1-surface, 2-volume, 3-pointcloud.
 */
SplitScalarField::SplitScalarField(int topo):SplitField(topo){
	m_name = "MiMMO.SplitScalarField";
}

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SplitScalarField::SplitScalarField(const bitpit::Config::Section & rootXML){
	
	std::string fallback_name = "ClassNONE";	
	std::string fallback_topo = "-1";	
	std::string input_name = rootXML.get("ClassName", fallback_name);
	std::string input_topo = rootXML.get("Topology", fallback_topo);
	input_name = bitpit::utils::trim(input_name);
	input_topo = bitpit::utils::trim(input_topo);
	
	int topo = std::stoi(input_topo);
	m_topo = std::max(1,topo);
	if (m_topo >3) m_topo = 1;
	
	m_name = "MiMMO.SplitScalarField";
	
	
	if(input_name == "MiMMO.SplitScalarField"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::SplitScalarField constructor. No valid xml data found"<<std::endl;
	};
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

