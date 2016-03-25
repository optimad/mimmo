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

/*!Default constructor of RotationBox
 */
RotationBox::RotationBox(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
	m_name = "MiMMO.RotationBox";
};

/*!Default destructor of RotationBox
 */
RotationBox::~RotationBox(){};

/*!Copy constructor of RotationBox.
 */
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
	double L = norm2(m_direction);
	for (int i=0; i<3; i++)
		m_direction[i] /= L;
}

void
RotationBox::setRotation(double alpha){
	setInput(alpha);
}

void
RotationBox::setAxes(dmatrix33E axes){
	m_axes = axes;
}

void
RotationBox::setAxesOrigin(darray3E axes_origin){
	m_axes_origin = axes_origin;
}

dmatrix33E
RotationBox::getRotatedAxes(){
	return(m_rotax);
}

darray3E
RotationBox::getRotatedOrigin(){
	return(m_rotax_origin);
}

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions, to be furnished by a pin to the child object.
 */
void
RotationBox::execute(){

	//Rotation of origin
	m_rotax_origin = {{0,0,0}};
	m_axes_origin -= m_origin;
	double alpha = *getInput<double>();
	//rodrigues formula
	m_rotax_origin = m_axes_origin * cos(alpha) +
			dotProduct(m_direction, m_axes_origin) * (1 - cos(alpha)) * m_direction +
			crossProduct(m_direction, m_axes_origin) * sin(alpha);

	m_rotax_origin += m_origin;
	m_axes_origin += m_origin;

	//rotation of axes
	m_rotax.fill(darray3E{{0,0,0}});
	//rodrigues formula
	for (int i=0; i<3; i++){
		m_rotax[i] = m_axes[i] * cos(alpha) +
				dotProduct(m_direction, m_axes[i]) * (1 - cos(alpha)) * m_direction +
				crossProduct(m_direction, m_axes[i]) * sin(alpha);
	}

	return;
};

