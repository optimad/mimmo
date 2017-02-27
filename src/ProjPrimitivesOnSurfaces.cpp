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

#include "ProjPrimitivesOnSurfaces.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!Default constructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(){
	m_name 		= "";
	m_nC = 1000;
	
}

/*!
 * Default destructor of ProjPrimitivesOnSurfaces.
 */
ProjPrimitivesOnSurfaces::~ProjPrimitivesOnSurfaces(){
	clear();
};

/*!Copy constructor of ProjPrimitivesOnSurfaces.Soft Copy of MimmoObject;
 */
ProjPrimitivesOnSurfaces::ProjPrimitivesOnSurfaces(const ProjPrimitivesOnSurfaces & other){
	*this = other;
};	

/*!
 * Assignement operator of ProjPrimitivesOnSurfaces. Soft copy of the target class (no resulting projected 
 * elements will be copied)
 */
ProjPrimitivesOnSurfaces & ProjPrimitivesOnSurfaces::operator=(const ProjPrimitivesOnSurfaces & other){
	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	
	m_topo = other.m_topo;
	m_nC = other.m_nC;
	m_buildBvTree = other.m_buildBvTree;
	m_buildKdTree = other.m_buildKdTree;
	
	//warning the internal data structure of projected element is not copied. Relaunch the execution eventually to fill it.
	return *this;
};

void
ProjPrimitivesOnSurfaces::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	
	built = (built && createPortOut<MimmoObject*, ProjPrimitivesOnSurfaces>(this, &mimmo::ProjPrimitivesOnSurfaces::getProjectedElement, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	m_arePortsBuilt = built;
}


/*!
 * Return kind of topology the primitive element support. 1-1D, 2-2D
 */
int 
ProjPrimitivesOnSurfaces::getTopology(){
	return m_topo;
}

/*!
 * Get the actual number of the cells that will be used to represent a discrete mesh of the primitive projection.
 */
int
ProjPrimitivesOnSurfaces::getProjElementTargetNCells(){
	return m_nC;
}


/*!
 * Get your current projected primitive as a 3D mesh
 */
MimmoObject *
ProjPrimitivesOnSurfaces::getProjectedElement(){
	return m_patch.get();
}



/*!
 * Set a target external surface mesh where primitive projection need to be performed.
 * Topology of the geometry must be of superficial type, so that MimmoObject::getType() must return 1, otherwise 
 * nothing will be set.
 * \param[in] geo  pointer to MimmoObject
 */
void
ProjPrimitivesOnSurfaces::setGeometry(MimmoObject* geo){
		if(geo->isEmpty()) return;
		if(geo->getType() != 1)	return;
		m_geometry = geo;
};

/*!It sets if the BvTree of the projected primitive element data structure has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
ProjPrimitivesOnSurfaces::setBuildBvTree(bool build){
	m_buildBvTree = build;
}

/*!It sets if the KdTree the projected primitive element data structure has to be built during execution.
* \param[in] build If true the KdTree is built in execution and stored in
* the related MimmoObject member.
*/
void
ProjPrimitivesOnSurfaces::setBuildKdTree(bool build){
	m_buildKdTree = build;
}

/*!It sets the number of the cells to be used to represent in a discrete mesh the primitive projection.
 * \param[in] nC number of cells
 */
void
ProjPrimitivesOnSurfaces::setProjElementTargetNCells(int nC){
	m_nC = nC;
}

/*!
 * Check if resulting primitive element is present or not.
 * True - no geometry evaluated, False otherwise.
 */
bool 
ProjPrimitivesOnSurfaces::isEmpty(){
	return (m_patch.get() == NULL);
}

/*!
 * Clear all stuffs in your class
 */
void
ProjPrimitivesOnSurfaces::clear(){
	m_patch.reset(nullptr);
	m_nC = 1000;
	BaseManipulation::clear();
};

/*!Execution command.
 * stitch together multiple geoemetry in the same object. 
 */
void
ProjPrimitivesOnSurfaces::execute(){
	if(m_geometry->isEmpty())	return;
	projection();
	
	if(!isEmpty()){
		if(m_buildBvTree)	m_patch->buildBvTree();
		if(m_buildKdTree)	m_patch->buildKdTree();
	}
};	

/*!
 * Plot resulting projected element in a vtu mesh file;
 */
void 
ProjPrimitivesOnSurfaces::plotOptionalResults(){
	if(isEmpty()) return;
	if(m_patch->getNCells() < 1)	return;
	std::string name = m_name + "_" + std::to_string(getClassCounter()) +  "_Patch";
	m_patch->getPatch()->write(name);
}


