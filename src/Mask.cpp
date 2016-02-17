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
#include "Mask.hpp"

///*!Default constructor of Mask
// */
Mask::Mask(){};
//
///*!Default destructor of Mask
// */
Mask::~Mask(){};

///*!Copy constructor of Mask.
// */
Mask::Mask(const Mask & other):BaseManipulation(other){};

/*!Assignement operator of Mask.
 */
Mask & Mask::operator=(const Mask & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
};

void
Mask::setCoords(dvecarr3E & coords){
	m_coords = coords;
};

void
Mask::setThresholds(darray3E & thres){
	m_thres = thres;
};

void
Mask::setForward(int i, bool forward){
	if (i >= 0 && i < 3) m_forward[i] = forward;
};

/*!Execution command. It applies the deformation given by the parent manipulation
 * to the linked geometry. After exec() the original geometry will be permanently modified.
 */
void
Mask::execute(){
//	recoverDisplacementsOut();
	for (int i=0; i<m_ndegout; i++){
		if (m_coords[i][0]>m_thres[0] && m_coords[i][1]>m_thres[1] && m_coords[i][2]>m_thres[2]){
			for (int j=0; j<3; j++){
				m_displout[i][j] = (1-m_forward[j])*m_displout[i][j];
			}
		}
		else{
			for (int j=0; j<3; j++){
				m_displout[i][j] = (m_forward[j])*m_displout[i][j];
			}
		}
	}
	return;
};
