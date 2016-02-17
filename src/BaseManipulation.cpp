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
#include "BaseManipulation.hpp"

using namespace std;

/*!Default constructor of BaseManipulation.
 * It sets to zero/null each member/pointer.
 */
BaseManipulation::BaseManipulation(){
	m_ndeg 			= 0;
	m_ndegout 		= 0;
	m_parent 		= NULL;
	m_geometry 		= NULL;
	m_child.clear();
	m_displ.clear();
	m_displout.clear();
};

/*!Custom constructor of BaseManipulation.
 * \param[in] geometry Pointer to target geometry to be linked.
 * \param[in] child Pointer to reference manipulator object to be linked (default value = NULL).
 */
BaseManipulation::BaseManipulation(MimmoObject* geometry, BaseManipulation* child){
	m_ndeg 			= 0;
	m_ndegout 		= 0;
	m_parent 		= NULL;
	if (child != NULL)	addChild(child);
	m_geometry 		= geometry;
	m_displ.clear();
	m_displout.clear();
};

/*!Custom constructor of BaseManipulation.
 * \param[in] child Pointer to reference manipulator object to be linked.
 */
BaseManipulation::BaseManipulation(BaseManipulation* child){
	m_ndeg 			= 0;
	m_ndegout 		= 0;
	m_parent 		= NULL;
	addChild(child);
	m_geometry 		= NULL;
	m_displ.clear();
	m_displout.clear();
};

/*!Default destructor of BaseManipulation.
 */
BaseManipulation::~BaseManipulation(){
	clear();
	m_ndeg = 0;
	m_ndegout = 0;
};

/*!Copy constructor of BaseManipulation.
 */
BaseManipulation::BaseManipulation(const BaseManipulation & other){
	m_ndeg 			= other.m_ndeg;
	m_displ 		= other.m_displ;
	m_ndegout 		= other.m_ndegout;
	m_displout 		= other.m_displout;
	m_child 		= other.m_child;
	m_parent 		= other.m_parent;
	m_geometry 		= other.m_geometry;
};

/*!Assignement operator of BaseManipulation.
 */
BaseManipulation & BaseManipulation::operator=(const BaseManipulation & other){
	m_ndeg 			= other.m_ndeg;
	m_displ 		= other.m_displ;
	m_ndegout 		= other.m_ndegout;
	m_displout 		= other.m_displout;
	m_parent 		= other.m_parent;
	m_child 		= other.m_child;
	m_geometry 		= other.m_geometry;
};

/*!It gets the number of degrees of freedom of the manipulator object.
 * \return #Degrees of freedom.
 */
uint32_t
BaseManipulation::getNDeg(){
	return m_ndeg;
};

/*!It gets the number of degrees of freedom output of the manipulator object.
 * \return #Degrees of freedom in output.
 */
uint32_t
BaseManipulation::getNDegOut(){
	return m_ndegout;
};

/*!It gets the displacement of the degree of freedom currently stored in the object.
 * \return Displacements of the degrees of freedom.
 */
dvecarr3E*
BaseManipulation::getDisplacements(){
	return &m_displ;
};

/*!It gets the displacement of the degree of freedom output of the object.
 * \return Displacements of the degrees of freedom in output.
 */
dvecarr3E*
BaseManipulation::getDisplacementsOut(){
	return &m_displout;
};

/*!It gets the manipulator object linked by this object.
 * \return Pointer to parent manipulator object.
 */
BaseManipulation*
BaseManipulation::getParent(){
	return m_parent;
};

/*!It gets the number of children linked to the manipulator object.
 * \return #Children.
 */
int
BaseManipulation::getNChild(){
	return m_child.size();
};

/*!It gets one child object linked by this object.
 * \return Pointer to i-th child manipulator object.
 */
BaseManipulation*
BaseManipulation::getChild(int i){
	if (i>=m_child.size()) return NULL;
	return m_child[i];
};

/*!It gets the geometry linked by the manipulator object.
 * \return Pointer to geometry to be deformed by the manipulator object.
 */
MimmoObject*
BaseManipulation::getGeometry(){
	return m_geometry;
};

/*!It sets the number of degrees of freedom of the manipulator object.
 * \param[in] ndeg #Degrees of freedom.
 */
void
BaseManipulation::setNDeg(uint32_t ndeg){
	m_ndeg = ndeg;
	m_displ.resize(m_ndeg);
};

/*!It sets the number of degrees of freedom output of the manipulator object.
 * \param[in] ndeg #Degrees of freedom output.
 */
void
BaseManipulation::setNDegOut(uint32_t ndeg){
	m_ndegout = ndeg;
	m_displout.resize(m_ndegout);
};

/*!It sets the displacement of the degree of freedom currently stored in the object.
 * \param[in] displacements Displacements of the degrees of freedom.
 */
void
BaseManipulation::setDisplacements(dvecarr3E & displacements){
	m_displ = displacements;
	m_ndeg = m_displ.size();
};

/*!It sets the displacement of the degree of freedom ouput of the object.
 * \param[in] displacements Displacements of the degrees of freedom in output.
 */
void
BaseManipulation::setDisplacementsOut(dvecarr3E & displacements){
	m_displout = displacements;
	m_ndegout = m_displout.size();
};

/*!It sets the manipulator object linked by this object.
 * \param[in] manipulator Pointer to parent manipulator object.
 */
void
BaseManipulation::setParent(BaseManipulation* parent){
	m_parent = parent;
};

/*!It adds a child manipulator object to the children linked by this object.
 * \param[in] child Pointer to child manipulator object.
 */
void
BaseManipulation::addChild(BaseManipulation* child){
	m_child.push_back(child);
	m_child[m_child.size()-1]->setParent(this);
};

/*!It sets the geometry linked by the manipulator object.
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
BaseManipulation::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!It clears the manipulator object linked by this object.
 * It sets to NULL the pointer to parent manipulator object.
 */
void
BaseManipulation::unsetParent(){
	m_parent = NULL;
};

/*!It clears the children objects linked by this object.
 * It sets to NULL the pointers to children manipulator object and clear the vector.
 */
void
BaseManipulation::unsetChild(){
	for (int i=0; i<m_child.size(); i++) m_child[i] = NULL;
	m_child.clear();
};

/*!It clears the displacement of the degree of freedom currently stored in the object.
 */
void
BaseManipulation::clearDisplacements(){
	m_displ.clear();
};

/*!It clears the displacement of the degree of freedom output of the object.
 */
void
BaseManipulation::clearDisplacementsOut(){
	m_displout.clear();
};

/*!It clears the pointer to the geometry linked by the object.
 */
void
BaseManipulation::unsetGeometry(){
	m_geometry = NULL;
};

/*!It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
BaseManipulation::clear(){
	unsetParent();
	unsetChild();
	unsetGeometry();
	clearDisplacements();
	clearDisplacementsOut();
};


/*!It recovers the information on the number of the degrees of freedom in input and their
 * displacements from the parent manipulator object.
 */
void
BaseManipulation::recoverDisplacementsIn(){
	if (m_parent == NULL) return;
	if (m_parent->getNDegOut() > 0){
		setNDeg(m_parent->getNDegOut());
		setDisplacements(*(m_parent->getDisplacementsOut()));
	}
};

/*!It recovers the information on the number of the degrees of freedom in output and their
 * initial displacements from the children manipulator objects.
 */
void
BaseManipulation::recoverDisplacementsOut(){
	//TODO PAY ATTENTION!!! IN THIS MOMENT IT IS ONLY THE LAST CHILD THAT GIVES THE NDEGOUT AND DISPLACEMENTSOUT !!
	for (int i=0; i<m_child.size(); i++){
		if (m_child[i] == NULL) return;
		if (m_child[i]->getNDeg() > 0 && m_ndegout == 0){
			setNDegOut(m_child[i]->getNDeg());
			setDisplacementsOut(*(m_child[i]->getDisplacements()));
		}
	}
};

/*!It initializes the number of degrees of freedom of the children with the
 * info of current object.
 */
void
BaseManipulation::initChild(){
	for (int i=0; i<m_child.size(); i++){
		if (m_child[i] == NULL) return;
		if (m_child[i]->getNDeg() == 0) m_child[i]->setNDeg(m_ndegout);
	}
};

/*!It updates the displacements of the degrees of freedom of the children with the
 * output of current object.
 */
void
BaseManipulation::updateChild(){
	for (int i=0; i<m_child.size(); i++){
		if (m_child[i] == NULL) return;
		m_child[i]->setDisplacements(m_displout);
	}
};

/*!Execution command. It runs recoverDisplacements applyFilters and execute.
 * execute is pure virtual and it has to be implemented in a derived class.
 */
void
BaseManipulation::exec(){
	initChild();
	recoverDisplacementsOut();
	execute();
	updateChild();
}

