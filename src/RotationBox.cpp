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
#include "RotationBox.hpp"

///*!Default constructor of RotationBox
// */
RotationBox::RotationBox(darray3E origin, darray3E direction){
	m_ndeg = 1;
	m_origin = origin;
	m_direction = direction;
};
//
///*!Default destructor of RotationBox
// */
RotationBox::~RotationBox(){};

///*!Copy constructor of RotationBox.
// */
RotationBox::RotationBox(const RotationBox & other):BaseManipulation(other){
	m_origin = other.m_origin;
	m_direction = other.m_direction;
};

/*!Assignement operator of RotationBox.
 */
RotationBox & RotationBox::operator=(const RotationBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_origin = other.m_origin;
	m_direction = other.m_direction;
	return(*this);
};

void
RotationBox::setAxis(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
}

void
RotationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

void
RotationBox::setDirection(darray3E direction){
	m_direction = direction;
	double L = sqrt(m_direction[0]*m_direction[0] + m_direction[1]*m_direction[1] + m_direction[2]*m_direction[2]);
	for (int i=0; i<3; i++)
		m_direction[i] /= L;
//		m_direction[i] /= norm2(m_direction);
}

void
RotationBox::setRotation(double alpha){
	m_displ.resize(1);
	m_displ[0] = { {alpha, 0 , 0} };
}

void
RotationBox::useInfo(){
	// 3 axes
	m_axes.resize(m_info->m_naxes);
	for (int i=0; i<m_info->m_naxes; i++){
		for (int j=0; j<m_info->m_naxes; j++){
			m_axes[i][j] = m_info->m_axes[i][j];
		}
	}
	for (int i=0; i<3; i++)
		m_axes_origin[i] = m_info->m_origin[i];
}

/*!Execution command. It modifies the coordinates of the origin given by the child manipulation object
 * with the rotation conditions. After exec() the original origin will be permanently modified.
 * Set the translated origin only for one child (the first one) and it has to be a FFDLattice
 * (static cast to use setOrigin method of basic shape).
 */
void
RotationBox::execute(){

	//Rotation of origin
	dvecarr3E rotated(1, {{0,0,0}});
	m_axes_origin -= m_origin;
	//rodrigues formula
	rotated[0] = m_axes_origin * cos(m_displ[0][0]) +
			dotProduct(m_direction, m_axes_origin) * (1 - cos(m_displ[0][0])) * m_direction +
			crossProduct(m_direction, m_axes_origin) * sin(m_displ[0][0]);

	rotated[0] += m_origin;
	if (m_child[0] != NULL){
		static_cast<FFDLattice*>(m_child[0])->changeOrigin(rotated[0]);
	}

	//rotation of axes
	rotated.clear();
	rotated.resize(3, {{0,0,0}});

	for (int i=0; i<3; i++){
		rotated[i] = m_axes[i] * cos(m_displ[0][0]) +
				dotProduct(m_direction, m_axes[i]) * (1 - cos(m_displ[0][0])) * m_direction +
				crossProduct(m_direction, m_axes[i]) * sin(m_displ[0][0]);
	}
	if (m_child[0] != NULL){
		static_cast<FFDLattice*>(m_child[0])->setRefSystem(rotated[0], rotated[1], rotated[2]);
	}

	return;
};

