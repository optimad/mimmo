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
#include "MaskFilter.hpp"

/*!Default constructor of MaskFilter
 */
MaskFilter::MaskFilter():BaseManipulation(){};

/*!Default destructor of MaskFilter
 */
MaskFilter::~MaskFilter(){};

/*!Copy constructor of MaskFilter.
 */
MaskFilter::MaskFilter(const MaskFilter & other):BaseManipulation(other){};

/*!Assignement operator of MaskFilter.
 */
MaskFilter & MaskFilter::operator=(const MaskFilter & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
};

void
MaskFilter::setCoords(dvecarr3E & coords){
	m_coords = coords;
};

void
MaskFilter::setThresholds(darray3E & thres){
	m_thres = thres;
};

void
MaskFilter::setForward(bool forward){
	m_forward = forward;
};

/*!It recovers the GEOMETRY displacements as result of the linked parent manipulator.
 * The recovered displacements are pushed in the geometry displacements member of the base manipulation class.
 */
void
MaskFilter::recoverDisplacements(){
	if (m_manipulator == NULL) return;
	setNDeg(m_manipulator->getNDeg());
	setDisplacements(*(m_manipulator->getDisplacements()));
	return;
};

/*!Execution command. It applies the deformation given by the parent manipulation
 * to the linked geometry. After exec() the original geometry will be permanently modified.
 */
void
MaskFilter::execute(){
	recoverDisplacements();
	for (int i=0; i<m_ndeg; i++){
		for (int j=0; j<3; j++){
			if (m_coords[i][j]>m_thres[j]){
				m_displ[i][j] = (1-m_forward)*m_displ[i][j];
			}
			else{
				m_displ[i][j] = (m_forward)*m_displ[i][j];
			}
		}
	}

	if (m_manipulator != NULL) m_manipulator->setDisplacements(m_displ);
	return;
};
