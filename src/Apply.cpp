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
#include "Apply.hpp"


using namespace mimmo;

REGISTER_MANIP(BaseManipulation, Apply, "MiMMO.Apply")

/*!Default constructor of Apply
 */
Apply::Apply():BaseManipulation(){
	m_name = "MiMMO.Apply";
	m_force = false;
	buildPorts();
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
Apply::Apply(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.Apply";
	buildPorts();
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "MiMMO.Apply"){
			absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::Apply constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of Apply
 */
Apply::~Apply(){};

/*!Copy constructor of Apply.
 */
Apply::Apply(const Apply & other){
	*this = other;
};

/*!Assignement operator of Apply.
 */
Apply & Apply::operator=(const Apply & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_force = other.m_force;
	return(*this);
};


/*! It builds the input/output ports of the object
 */
void
Apply::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, Apply>(this, &Apply::setInput, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, Apply>(this, &BaseManipulation::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<MimmoObject*, Apply>(this, &BaseManipulation::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};

/*!
 * Return true, if rebuilding of search trees of your target geometry of class MimmoObject is forced by the User
 */
bool	Apply::getRefreshGeometryTrees(){
	return m_force;
}

/*!
 * If set true, forces rebuilding of search trees of your target geometry of class MimmoObject
 */
void	Apply::setRefreshGeometryTrees(bool force){
	m_force = force;
}

/*!It sets the displacements input.
 * \param[in] input Input displacements of the geometry vertices.
 */
void
Apply::setInput(dvecarr3E input){
	m_input = input;
}

/*!Execution command.
 * It applies the deformation stored in the input of base class (casting the input
 * for apply object to dvecarr3E) to the linked geometry.
 * After exec() the original geometry will be permanently modified.
 */
void
Apply::execute(){
	if (getGeometry() == NULL) return;
	
	dvecarr3E vertex = getGeometry()->getVertexCoords();
	long nv = getGeometry()->getNVertex();
	nv = long(std::min(int(nv), int(m_input.size())));
	livector1D & idmap = getGeometry()->getMapData();
	for (long i=0; i<nv; i++){
		vertex[i] += m_input[i];
		getGeometry()->modifyVertex(vertex[i], idmap[i]);
	}
	
	if(m_force){
		getGeometry()->buildBvTree();
		getGeometry()->buildKdTree();
	}
	
	return;
};


/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read only RefreshGeometryTrees parameter, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 * --> Absorbing data:
 * 		Priority  : uint marking priority in multi-chain execution; 
 * 		RefreshGeometryTrees: boolean 0,1 activate rebuilding of search trees on target geometry
 * 		
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
 void Apply::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	 
	std::string input; 
	
	if(slotXML.hasOption("Priority")){
		input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	

	if(slotXML.hasOption("RefreshGeometryTrees")){
		input = slotXML.get("RefreshGeometryTrees");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setRefreshGeometryTrees(value);
	}; 
	
};

/*!
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML;
 * The class write only RefreshGeometryTrees parameter, if it is different from its default value, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 *  * --> Flushing data// how to write it on XML:
 * 		ClassName : name of the class as "MiMMO.Apply"
 * 		Priority  : uint marking priority in multi-chain execution; 
 * 		RefreshGeometryTrees: boolean 0,1 activate rebuilding of search trees on target geometry
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void Apply::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
 	slotXML.set("Priority", std::to_string(getPriority()));
	
	bool value = getRefreshGeometryTrees();
	
	std::string towrite = std::to_string(value);
	
	slotXML.set("RefreshGeometryTrees", towrite);
};	





