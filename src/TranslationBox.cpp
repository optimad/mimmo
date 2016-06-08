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
#include "TranslationBox.hpp"

using namespace mimmo;

/*!Default constructor of TranslationBox
 */
TranslationBox::TranslationBox(darray3E direction){
	m_direction = direction;
	m_name = "MiMMO.TranslationBox";
};

/*!Default destructor of TranslationBox
 */
TranslationBox::~TranslationBox(){};

/*!Copy constructor of TranslationBox.
 */
TranslationBox::TranslationBox(const TranslationBox & other):BaseManipulation(other){
	m_direction = other.m_direction;
};

/*!Assignement operator of TranslationBox.
 */
TranslationBox & TranslationBox::operator=(const TranslationBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_direction = other.m_direction;
	return(*this);
};

/*! It builds the input/output ports of the object
 */
void TranslationBox::buildPorts(){
	bool built = true;
	built = (built && createPortIn<darray3E, TranslationBox>(&m_origin, POINT, 20 , {POINT2}));
	built = (built && createPortIn<darray3E, TranslationBox>(&m_direction, AXIS, 21));
	built = (built && createPortIn<double, TranslationBox>(&m_alpha, VALUED, 30, {VALUED2}));
	built = (built && createPortOut<darray3E, TranslationBox>(this, &mimmo::TranslationBox::getOrigin, POINT, 20));
	built = (built && createPortOut<darray3E, TranslationBox>(this, &mimmo::TranslationBox::getDirection, AXIS, 21));
	built = (built && createPortOut<double, TranslationBox>(this, &mimmo::TranslationBox::getTranslation, VALUED, 30));
	m_arePortsBuilt = built;
};

/*!It gets the direction of the translation.
 * \return Direction of translation.
 */
darray3E
TranslationBox::getDirection(){
	return(m_direction);
}

/*!It gets the value of the translation.
 * \return Value of translation.
 */
double
TranslationBox::getTranslation(){
	return(m_alpha);
}

/*!It gets the original position of the point to be translated (before the execution)
 * or the position of the translated point (after the execution of the object).
 * \return Position of the point.
 */
darray3E
TranslationBox::getOrigin(){
	return(m_origin);
}

/*!It sets the direction of the translation.
 * \param[in] direction Direction of translation.
 */
void
TranslationBox::setDirection(darray3E direction){
	m_direction = direction;
}

/*!It sets the value of the translation.
 * \param[in] alpha Value of translation.
 */
void
TranslationBox::setTranslation(double alpha){
	m_alpha = alpha;
}

/*!It sets the original coordinates of the point to be translated.
 * \param[in] origin Position of the point.
 */
void
TranslationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!Execution command. It modifies the coordinates of the origin
 * with the translation conditions.
 * The result of the translation is stored in member result of base class
 * and in the member m_origin.
 * After exec() the original point coordinates will be permanently modified.
 */
void
TranslationBox::execute(){
	for (int i=0; i<3; i++){
			m_origin[i] += m_alpha * m_direction[i];
	}
	return;
};
