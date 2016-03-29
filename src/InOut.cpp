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
#include "InOut.hpp"
#include "BaseManipulation.hpp"

using namespace std;
using namespace mimmo;

//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION								//
//==============================================================//

/*!Default constructor of InOut
*/
InOut::InOut(){};

/*!Default destructor of InOut
*/
InOut::~InOut(){
	m_objLink	= NULL;
};

/*!Copy constructor of InOut.
*/
InOut::InOut(const InOut & other){
	m_objLink 	= other.m_objLink;
};

/*!Assignement operator of InOut.
 */
InOut & InOut::operator=(const InOut & other){
	m_objLink 	= other.m_objLink;
	return (*this);
};

/*!Compare operator of InOut.
 */
bool InOut::operator==(const InOut & other){
	return(m_objLink == other.m_objLink);
};

/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
BaseManipulation*
InOut::getLink(){
	return(m_objLink);
}

