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

/*!Default constructor of BaseManipulation.
 * It sets to zero/null each member/pointer.
 */
BaseManipulation::BaseManipulation(){
	m_displ.clear();
	m_ndeg = 0;
	m_manipulator = NULL;
	m_geometry = NULL;
	m_gdispl.clear();
};

/*!Custom constructor of BaseManipulation.
 * \param[in] geometry Pointer to target geometry to be linked.
 * \param[in] parent Pointer to reference manipulator object to be linked (default value = NULL).
 */
BaseManipulation::BaseManipulation(MimmoObject* geometry, BaseManipulation* parent){
	m_displ.clear();
	m_ndeg = 0;
	m_manipulator = parent;
	m_geometry = geometry;
	m_gdispl.clear();
};

/*!Custom constructor of BaseManipulation.
 * \param[in] parent Pointer to reference manipulator object to be linked.
 */
BaseManipulation::BaseManipulation(BaseManipulation* parent){
	m_displ.clear();
	m_ndeg 			= 0;
	m_manipulator 	= parent;
	m_geometry 		= NULL;
	m_gdispl.clear();
};

/*!Default destructor of BaseManipulation.
 */
BaseManipulation::~BaseManipulation(){
	clear();
};

/*!Copy constructor of BaseManipulation.
 */
BaseManipulation::BaseManipulation(const BaseManipulation & other){
	m_displ 		= other.m_displ;
	m_ndeg 			= other.m_ndeg;
	m_manipulator 	= other.m_manipulator;
	m_geometry 		= other.m_geometry;
	m_gdispl 		= other.m_gdispl;
};

/*!Assignement operator of BaseManipulation.
 */
BaseManipulation & BaseManipulation::operator=(const BaseManipulation & other){
	m_displ 		= other.m_displ;
	m_ndeg 			= other.m_ndeg;
	m_manipulator 	= other.m_manipulator;
	m_geometry 		= other.m_geometry;
	m_gdispl 		= other.m_gdispl;
};

/*!It gets the number of degrees of freedom of the manipulator object.
 * \return #Degrees of freedom.
 */
uint32_t
BaseManipulation::getNDeg(){
	return m_ndeg;
};

/*!It gets the displacement of the degree of freedom currently stored in the object.
 * \return Displacements of the degrees of freedom.
 */
dvecarr3E*
BaseManipulation::getDisplacements(){
	return &m_displ;
};

/*!It gets the manipulator object linked by this object.
 * \return Pointer to parent manipulator object.
 */
BaseManipulation*
BaseManipulation::getManipulator(){
	return m_manipulator;
};

/*!It gets the geometry linked by the manipulator object.
 * \return Pointer to geometry to be deformed by the manipulator object.
 */
MimmoObject*
BaseManipulation::getGeometry(){
	return m_geometry;
};

/*!It gets the displacements of the geometry linked by the manipulator object.
 * \return Displacements to be applied to the geometry linked by the manipulator object.
 */
dvecarr3E*
BaseManipulation::getGeometryDisplacements(){
	return &m_gdispl;
};

/*!It sets the number of degrees of freedom of the manipulator object.
 * \param[in] ndeg #Degrees of freedom.
 */
void
BaseManipulation::setNDeg(uint32_t ndeg){
	m_ndeg = ndeg;
};

/*!It sets the displacement of the degree of freedom currently stored in the object.
 * \param[in] displacements Displacements of the degrees of freedom.
 */
void
BaseManipulation::setDisplacements(dvecarr3E & displacements){
	m_displ = displacements;
};

/*!It sets the manipulator object linked by this object.
 * \param[in] manipulator Pointer to parent manipulator object.
 */
void
BaseManipulation::setManipulator(BaseManipulation* manipulator){
	m_manipulator = manipulator;
};

/*!It sets the geometry linked by the manipulator object.
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
BaseManipulation::setGeometry(MimmoObject* geometry){
	m_geometry = geometry;
};

/*!It sets the displacements of the geometry linked by the manipulator object.
 * \param[in] gdisplacements Displacements to be applied to the geometry linked by the manipulator object.
 */
void
BaseManipulation::setGeometryDisplacements(dvecarr3E & gdisplacements){
	m_gdispl = gdisplacements;
};

/*!It clears the manipulator object linked by this object.
 * It sets to NULL the pointer to parent manipulator object.
 */
void
BaseManipulation::unsetManipulator(){
	m_manipulator = NULL;
};

/*!It clears the displacement of the degree of freedom currently stored in the object.
 */
void
BaseManipulation::clearDisplacements(){
	m_displ.clear();
};

/*!It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
BaseManipulation::clear(){
	m_displ.clear();
	m_manipulator = NULL;
	m_geometry = NULL;
	m_gdispl.clear();
};

/*!It recovers the information on the number of the degrees of freedom and their
 * displacements from the parent manipulator object.
 */
void
BaseManipulation::recoverDisplacements(){
	setNDeg(m_manipulator->getNDeg());
	setDisplacements(*(m_manipulator->getDisplacements()));
};
