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
\*---------------------------------------------------------------------------*/
#include "MultiApply.hpp"

namespace mimmo{

/*!
 * Default constructor of Apply
 */
MultiApply::MultiApply():BaseManipulation(){
	m_name = "mimmo.MultiApply";
	m_force = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
MultiApply::MultiApply(const bitpit::Config::Section & rootXML){
	
	m_name = "mimmo.MultiApply";
	m_force = false;
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "mimmo.MultiApply"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::MultiApply constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of Apply
 */
MultiApply::~MultiApply(){};

/*!Copy constructor of Apply.
 */
MultiApply::MultiApply(const MultiApply & other):BaseManipulation(){
	*this = other;
};

/*!Assignement operator of Apply.
 */
MultiApply & MultiApply::operator=(const MultiApply & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_force = other.m_force;
	m_input = other.m_input;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
MultiApply::buildPorts(){
	bool built = true;
	built = (built && createPortIn<std::unordered_map<MimmoObject*, dvecarr3E*>, MultiApply>(this, &MultiApply::setInputList, PortType::M_UMGEOVFD, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::MIMMO_VECARR3FLOAT_));
	m_arePortsBuilt = built;
};

/*!
 * Return true, if rebuilding of search trees of your target geometries of class MimmoObject are forced by the User
 */
bool	MultiApply::getRefreshGeometryTrees(){
	return m_force;
}



/*!
 * If set true, forces rebuilding of search trees of your target geometries of class MimmoObject
 */
void	MultiApply::setRefreshGeometryTrees(bool force){
	m_force = force;

}

/*! Set the input of the class as std::pair of target MimmoObject * and
 * relative displacement field dvecarr3E *. If target geometry already exist or it is NULL
 * no input is added to the class.
 */
void
MultiApply::addInput(std::pair<MimmoObject*, dvecarr3E*> input){

	if (input.first == NULL || m_input.count(input.first) > 0 ) return;
	m_input[input.first] = input.second;
}


/*! Set a whole list of inputs of the class as an unorderd map of target MimmoObject * as key and
 * relative displacement field dvecarr3E * as argument. If a target geometry already exists in the member list or it is NULL
 * no input is added to the class. Multiple insertion of list of inputs are appended in the unique member list.
 */
void
MultiApply::setInputList(std::unordered_map<MimmoObject*, dvecarr3E*> input){
	
	for(auto & val : input){
		
		addInput(val);
	}
}

/*!
 * Clear the input list stored into the class
 */
void
MultiApply::clearList(){
	m_input.clear();
}

/*!Execution command.
 * It applies the deformations stored in the class  to their relative linked geometry.
 * After exec() the original geometries will be permanently modified.
 */
void
MultiApply::execute(){
	

	if (m_input.empty()) return;


	for(auto val: m_input){	
		

		
		dvecarr3E vertex = val.first->getVertexCoords();
		long nv = val.first->getNVertex();
		nv = long(std::min(int(nv), int(val.second->size())));

		livector1D & idmap = val.first->getMapData();
		for (long i=0; i<nv; i++){

			vertex[i] += (*(val.second))[i];
			val.first->modifyVertex(vertex[i], idmap[i]);
		}	

		if(m_force){
			val.first->buildBvTree();
			val.first->buildKdTree();
		}
	}
	return;
};


/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read only RefreshGeometryTrees parameter, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 * --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>RefreshGeometryTrees</B>: 0/1 to force update of trees for current linked geometries
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
 void MultiApply::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	 BITPIT_UNUSED(name);
	 
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
		std::string input = slotXML.get("RefreshGeometryTrees");
	}; 
	
	bool value = false;
	if(!input.empty()){
		std::stringstream ss(bitpit::utils::trim(input));
		ss >> value;
	}
	
	setRefreshGeometryTrees(value);
};

/*!
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML;
 * The class write only RefreshGeometryTrees parameter, if it is different from its default value, 
 * while Input and Geometry parameters are meant to be passed only through Port linking.
 * 
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.MultiApply"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>RefreshGeometryTrees</B>: 0/1 to force update of trees for current linked geometries
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * * \param[in] name   name associated to the slot
 */
void MultiApply::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	bool value = getRefreshGeometryTrees();
	
	std::string towrite = std::to_string(value);
	
	slotXML.set("RefreshGeometryTrees", towrite);
};	


}


