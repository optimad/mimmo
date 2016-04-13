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
#include "DofOut.hpp"

using namespace std;
using namespace mimmo;

//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION								//
//==============================================================//

/*!Default constructor of DofOut
*/
DofOut::DofOut(){};

/*!Default destructor of DofOut
*/
DofOut::~DofOut(){
	m_obj	= NULL;
};

/*!Copy constructor of DofOut.
*/
DofOut::DofOut(const DofOut & other){
	m_obj 	= other.m_obj;
};

/*!Assignement operator of DofOut.
 */
DofOut & DofOut::operator=(const DofOut & other){
	m_obj 	= other.m_obj;
	return (*this);
};

/*!Compare operator of DofOut.
 */
bool DofOut::operator==(const DofOut & other){
	return(m_obj == other.m_obj);
};

/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
BaseManipulation*
DofOut::getLink(){
	return(m_obj);
}

int
DofOut::getNuse(){
	return(m_nuse);
}

int
DofOut::getNgdof(){
	return(m_nglob);
}


bool
DofOut::isActive(int j){
	if (j>=getNgdof()) return false;
	return m_actives[j];
}






