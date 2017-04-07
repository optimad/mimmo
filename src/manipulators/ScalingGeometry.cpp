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
#include "ScalingGeometry.hpp"

namespace mimmo{


/*!
 * Default constructor of ScalingGeometry
 */
ScalingGeometry::ScalingGeometry(darray3E scaling){
    m_scaling = scaling;
	m_name = "mimmo.ScalingGeometry";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
ScalingGeometry::ScalingGeometry(const bitpit::Config::Section & rootXML){
	
    m_scaling.fill(1.0);
	m_name = "mimmo.ScalingGeometry";
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "mimmo.ScalingGeometry"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::ScalingGeometry constructor. No valid xml data found"<<std::endl;
	};
}

/*!Default destructor of ScalingGeometry
 */
ScalingGeometry::~ScalingGeometry(){};

/*!Copy constructor of ScalingGeometry.
 */
ScalingGeometry::ScalingGeometry(const ScalingGeometry & other):BaseManipulation(other){
    m_scaling = other.m_scaling;
};

/*!Assignement operator of ScalingGeometry.
 */
ScalingGeometry & ScalingGeometry::operator=(const ScalingGeometry & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_scaling = other.m_scaling;
	return(*this);
};


/*! It builds the input/output ports of the object
 */
void ScalingGeometry::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray3E, ScalingGeometry>(&m_scaling, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<dvector1D, ScalingGeometry>(this, &mimmo::ScalingGeometry::setFilter, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<MimmoObject*, ScalingGeometry>(&m_geometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
    built = (built && createPortOut<dvecarr3E, ScalingGeometry>(this, &mimmo::ScalingGeometry::getDisplacements, PortType::M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	m_arePortsBuilt = built;
};

/*!It sets the scaling factors of each axis.
 * \param[in] scaling Scaling factor values for x, y and z absolute axis.
 */
void
ScalingGeometry::setScaling(darray3E scaling){
    m_scaling = scaling;
}

/*!It sets the filter field to modulate the displacements of the vertices
 * of the target geometry.
 * \param[in] filter Filter field defined on geometry vertices.
 */
void
ScalingGeometry::setFilter(dvector1D filter){
	m_filter = filter;
}

/*!
 * Return actual computed displacements field (if any) for the geometry linked.
 * @return  The computed deformation field on the vertices of the linked geometry
 */
dvecarr3E
ScalingGeometry::getDisplacements(){
    return m_displ;
};

/*!Execution command. It perform the scaling by computing the displacements
 * of the points of the geometry. It applies a filter field eventually set as input.
 */
void
ScalingGeometry::execute(){

    int nV = m_geometry->getNVertex();
    m_displ.resize(nV);
    m_filter.resize(nV, 1.0);


    long ID;
    int idx;
    liimap mapID = m_geometry->getMapDataInv();

    //computing centroid
    darray3E center;
    center.fill(0.0);
    for (auto vertex : m_geometry->getVertices()){
        ID = vertex.getId();
        idx = mapID[ID];
        center += vertex.getCoords() / double(nV);
    }

    for (auto vertex : m_geometry->getVertices()){
        ID = vertex.getId();
        idx = mapID[ID];

        darray3E coords = vertex.getCoords();
        m_displ[idx] = ( m_scaling*(coords - center) + center ) * m_filter[idx] - coords;
    }
    return;
};

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform scaling
 * a geometry. Filter field, geometry and resulting displacements are passed mandatorily through ports
 * 
 * --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution;	
 * - <B>Scaling</B>: scaling factor values for each cartesian axis.
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void ScalingGeometry::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
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
	
	if(slotXML.hasOption("Scaling")){
		std::string input = slotXML.get("Scaling");
		input = bitpit::utils::trim(input);
		darray3E temp = {{0.0,0.0,0.0}};
		if(!input.empty()){
			std::stringstream ss(input);
			for(auto &val : temp) ss>>val;
		}
		setScaling(temp);
	} 
	
};	
/*!
 * Write settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read essential parameters to perform scaling of a
 * geometry.Filter field, geometry and resulting displacements are passed mandatorily through ports
 * 
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.ScalingGeometry"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;	
 * - <B>Scaling</B>: scaling factor values for each cartesian axis.
 *
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void ScalingGeometry::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	{
		std::stringstream ss;
		ss<<std::scientific<<m_scaling[0]<<'\t'<<m_scaling[1]<<'\t'<<m_scaling[2];
		slotXML.set("Scaling", ss.str());
	}

};	

}
