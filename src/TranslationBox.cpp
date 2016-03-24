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

///*!Default constructor of TranslationBox
// */
TranslationBox::TranslationBox(darray3E direction){
	m_direction = direction;
	m_name = "MiMMO.TranslationBox";
};
//
///*!Default destructor of TranslationBox
// */
TranslationBox::~TranslationBox(){};

///*!Copy constructor of TranslationBox.
// */
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

darray3E
TranslationBox::getDirection(){
	return(m_direction);
}

double
TranslationBox::getTranslation(){
	return(m_alpha);
}

darray3E
TranslationBox::getOrigin(){
	return(m_origin);
}

void
TranslationBox::setDirection(darray3E direction){
	m_direction = direction;
}

void
TranslationBox::setTranslation(double alpha){
	m_alpha = alpha;
}

void
TranslationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

/*!Execution command. It modifies the coordinates of the origin given by the child manipulation object
 * with the translation conditions. After exec() the original origin will be permanently modified.
 * Set the translated origin only for one child (the first one) and it has to be a FFDLattice
 * (static cast to use setOrigin method of basic shape).
 */
void
TranslationBox::execute(){
	for (int i=0; i<3; i++){
			m_origin[i] += m_alpha * m_direction[i];
	}
	setResult(m_origin);
	return;
};
