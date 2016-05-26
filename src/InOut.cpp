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

/*!Default constructor of PinOut
*/
PinOut::PinOut(){};

/*!Default destructor of PinOut
*/
PinOut::~PinOut(){};

/*!Copy constructor of PinOut.
*/
PinOut::PinOut(const PinOut & other){
	m_objLink 	= other.m_objLink;
	m_obuffer	= other.m_obuffer;
	m_portLink	= other.m_portLink;
	return;
};

/*!Assignement operator of PinOut.
 */
PinOut & PinOut::operator=(const PinOut & other){
	m_objLink 	= other.m_objLink;
	m_obuffer	= other.m_obuffer;
	m_portLink	= other.m_portLink;
	return (*this);
};

/*!Compare operator of PinOut.
 */
bool PinOut::operator==(const PinOut & other){
	bool check = true;
	check = check && (m_portLink == other.m_portLink);
	check = check && (m_objLink == other.m_objLink);
	return(check);
};

/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
std::vector<mimmo::BaseManipulation*>
PinOut::getLink(){
	return(m_objLink);
}

std::vector<int>
PinOut::getPortLink(){
	return(m_portLink);
}

/*! Execution of the pin.
 * All the pins of an object are called in execute of the owner after its own execution.
 */
void
mimmo::PinOut::exec(){
	if (m_objLink.size() > 0){
		writeBuffer();
		bitpit::IBinaryStream input(m_obuffer.data());
		for (int j=0; j<m_objLink.size(); j++){
			if (m_objLink[j] != NULL){
				m_objLink[j]->getPinsIn()[m_portLink[j]]->m_ibuffer = input;
				m_objLink[j]->readBuffer(m_portLink[j]);
			}
		}
	}
};

/*!Default constructor of PinIn
*/
PinIn::PinIn(){};

/*!Default destructor of PinIn
*/
PinIn::~PinIn(){
	m_objLink	= NULL;
};

/*!Copy constructor of PinIn.
*/
PinIn::PinIn(const PinIn & other){
	m_objLink 	= other.m_objLink;
	m_ibuffer	= other.m_ibuffer;
	return;
};

/*!Assignement operator of PinIn.
 */
PinIn & PinIn::operator=(const PinIn & other){
	m_objLink 	= other.m_objLink;
	m_ibuffer	= other.m_ibuffer;
	return (*this);
};

/*!Compare operator of PinIn.
 */
bool PinIn::operator==(const PinIn & other){
	bool check = true;
	check = check && (m_objLink == other.m_objLink);
	return(check);
};

/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
BaseManipulation*
PinIn::getLink(){
	return(m_objLink);
}
