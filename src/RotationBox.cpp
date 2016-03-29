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

using namespace mimmo;

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

/*!It sets the origin and direction of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationBox::setAxis(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
}

/*!It sets the origin of the rotation axis.
 * \param[in] origin Origin of rotation axis.
 */
void
RotationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!It sets the direction of the rotation axis.
 * \param[in] direction Direction of rotation axis.
 */
void
RotationBox::setDirection(darray3E direction){
	m_direction = direction;
	double L = norm2(m_direction);
	for (int i=0; i<3; i++)
		m_direction[i] /= L;
}

/*!It sets the value of the rotation.
 * \param[in] alpha Value of rotation axis.
 */
void
RotationBox::setRotation(double alpha){
	setInput(alpha);
}

/*!It sets the reference system to be rotated.
 * \param[in] axes Original reference system.
 */
void
RotationBox::setAxes(dmatrix33E axes){
	m_axes = axes;
}

/*!It sets the origin of the reference system to be rotated.
 * \param[in] axes_origin Origin of reference system.
 */
void
RotationBox::setAxesOrigin(darray3E axes_origin){
	m_axes_origin = axes_origin;
}

/*!It gets the rotated reference system.
 * \return Rotated reference system.
 */
dmatrix33E
RotationBox::getRotatedAxes(){
	return(m_rotax);
}

/*!It gets the rotated origin of the reference system.
 * \return Rotated origin of reference system.
 */
darray3E
RotationBox::getRotatedOrigin(){
	return(m_rotax_origin);
}

/*!Execution command. It saves in "rot"-terms the modified axes and origin, by the
 * rotation conditions. This terms can be recovered and passed by a pin to a child object
 * by the related get-methods.
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

