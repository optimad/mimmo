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

/*!Execution command. It modifies the displacements given by the child manipulation object
 * with the masking conditions. After exec() the original displacements will be permanently modified.
 */
void
Mask::execute(){
	for (int j=0; j<getNChild(); j++){
		dvecarr3E* displout = getDisplacementsOut(j);
		for (int i=0; i<getNDegOut(j); i++){
			if (m_coords[i][0]>m_thres[0] && m_coords[i][1]>m_thres[1] && m_coords[i][2]>m_thres[2]){
				for (int j=0; j<3; j++){
					(*displout)[i][j] = (1-m_forward[j])*(*displout)[i][j];
				}
			}
			else{
				for (int j=0; j<3; j++){
					(*displout)[i][j] = (m_forward[j])*(*displout)[i][j];
				}
			}
		}
	}
	return;
};
