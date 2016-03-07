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
#include "DeformBox.hpp"

///*!Default constructor of DeformBox
// */
DeformBox::DeformBox(){};
//
///*!Default destructor of DeformBox
// */
DeformBox::~DeformBox(){};

///*!Copy constructor of DeformBox.
// */
DeformBox::DeformBox(const DeformBox & other):BaseManipulation(other){};

/*!Assignement operator of DeformBox.
 */
DeformBox & DeformBox::operator=(const DeformBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
};

void
DeformBox::useInfo(){
	if (m_ndeg !=  m_info->m_coords.size() || m_info->m_naxes != 3){
		std::cout << "Incoherent Size ---> end of process " << std::endl;
		exit(1001);
	}
	m_coords.resize(m_info->m_coords.size());
	for (int i=0; i<m_ndeg; i++){
		for (int j=0; j<3; j++){
			m_coords[i][j] = m_info->m_coords[i][j];
		}
	}
}

/*!Execution command. It modifies the coordinates given by the child manipulation object
 * with the transform conditions. After exec() the original coordinates will be permanently modified.
 */
void
DeformBox::execute(){
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			m_coords[i][j] = dotProduct(m_displ[j], m_coords[i]);
		}
	}
	for (int i=0; i<m_child.size(); i++){
		if (m_child[i] != NULL){
//			*(static_cast<FFDLatticeBox*>(m_child[i]))->setCoords(m_coords);
		}
	}
	return;
};
