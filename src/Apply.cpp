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

/*!Default constructor of Apply
 */
Apply::Apply():BaseManipulation(){
	m_name = "MiMMO.Apply";
	buildPorts();
};

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
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void
Apply::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, Apply>(this, &Apply::setInput, M_GDISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<MimmoObject*, Apply>(this, &BaseManipulation::setGeometry, M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
};

/*!
 * If set true, forces rebuilding of search trees of your target geometry of class MimmoObject
 */
void	Apply::setRefreshGeometryTrees(bool force){
	if(getGeometry() == NULL) return;
	getGeometry()->buildBvTree();
	getGeometry()->buildKdTree();
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
	return;
};
