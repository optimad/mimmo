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

#include "RBFBox.hpp"
#include "LinearAlgebra.hpp"
#include "lapacke.h"

#include <chrono>

using namespace std::chrono;

using namespace std;
using namespace mimmo;

/*! Basic Constructor. Doing nothing.*/
RBFBox::RBFBox(){
	m_name = "MiMMO.RBFBox";
	m_origin.fill(0.0);
	m_span.fill(1.0);
	int counter = 0;
	for(auto &val : m_axes)	{
		val.fill(0.0);
		val[counter] = 1.0;
		++counter;
	}	
    m_nodes.clear();
    m_suppR = 0.0;
};

/*! Destructor */
RBFBox::~RBFBox(){};

/*! Copy Constructor
 *\param[in] other RBFBox where copy from
 */
RBFBox::RBFBox(const RBFBox & other){
	*this = other;
};

/*! Copy Operator
 * \param[in] other RBFBox where copy from
 */
RBFBox & RBFBox::operator=(const RBFBox & other){

	*(static_cast<BaseManipulation *>(this))  = *(static_cast<const BaseManipulation *>(&other));
	m_origin = other.m_origin;
	m_span   = other.m_span;
    m_axes = other.m_axes;
    m_nodes = other.m_nodes;
    m_suppR = other.m_suppR;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void RBFBox::buildPorts(){

	bool built = true;
//creating input ports	
    built = (built && createPortIn<dvecarr3E, RBFBox>(this, &mimmo::RBFBox::setNode, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
    built = (built && createPortIn<double, RBFBox>(this, &mimmo::RBFBox::setSupportRadius, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
// creating output ports
	built = (built && createPortOut<darray3E, RBFBox>(this, &mimmo::RBFBox::getOrigin, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dmatrix33E, RBFBox>(this, &mimmo::RBFBox::getAxes, PortType::M_AXES, mimmo::pin::containerTAG::ARR3ARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, RBFBox>(this, &mimmo::RBFBox::getSpan, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	
	m_arePortsBuilt = built;
};

/*!Clean all stuffs in your class */
void RBFBox::clearRBFBox(){
	clear(); //base manipulation stuff clear
	m_origin.fill(0.0);
	m_span.fill(1.0);
	int counter = 0;
	for(auto &val : m_axes)	{
		val.fill(0.0);
		val[counter] = 1.0;
		++counter;
	}	
	m_nodes.clear();
	m_suppR = 0.0;
	
};

/*! 
 * Return the origin of the RBFBox.
 * \return Number of control nodes
 */
darray3E	RBFBox::getOrigin(){
	return(m_origin);
}

/*! 
 * Return the span of the RBFBox.
 * \return Number of control nodes
 */
darray3E	RBFBox::getSpan(){
	return(m_span);
}


/*! 
 * Return the oriented axes of the RBFBox.
 * \return Number of control nodes
 */
dmatrix33E	RBFBox::getAxes(){
	return(m_axes);
}

/*!Set a list of RBF points as control nodes.
 * @param[in] node coordinates of control points.
 */
void RBFBox::setNode(dvecarr3E nodes){
    m_nodes = nodes;
    return;
};

/*! Set the value of the support radius R of RBF kernel functions.
 * @param[in] suppR_ new value of support radius.
 */
void
RBFBox::setSupportRadius(double suppR_){
    m_suppR = std::fmax(0.0,suppR_);
}

/*!Execute your object, calculate the RBFBox of the input set of RBFs (+1% of support radius).
 *
 */
void 		RBFBox::execute(){

    int np = m_nodes.size();

    darray3E min, max;
    min.fill(1.0e+18);
    max.fill(-1.0e+18);

    for(int i=0; i<np; i++){
        for(int j=0; i<3; j++){
            if (m_nodes[i][j] - m_suppR*1.01 < min[j]) min[j] = m_nodes[i][j];
            if (m_nodes[i][j] + m_suppR*1.01 > max[j]) max[j] = m_nodes[i][j];
        }
    }

    m_span = max - min;
    m_origin = min + m_span * 0.5;

};

