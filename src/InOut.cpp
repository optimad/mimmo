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


/*!
	Output stream operator for dvecarr3E

	\param[in] buffer is the output stream
	\param[in] var is the element to be streamed
	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const dvecarr3E &var)
{
	int nP = var.size();
	buffer << nP;
	for (int i = 0; i < nP; ++i) {
		for (int j = 0; j < 3; ++j) {
			buffer << var[i][j];
		}
	}

	return buffer;
}


/*!
	Input stream operator for dvecarr3E

	\param[in] buffer is the input stream
	\param[in] var is the element to be streamed
	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, dvecarr3E &var)
{
	int nP;
	buffer >> nP;
	var.resize(nP);
	for (int i = 0; i < nP; ++i) {
		for (int j = 0; j < 3; ++j) {
			buffer >> var[i][j];
		}
	}

	return buffer;
}

//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION								//
//==============================================================//

/*!Default constructor of PortOut
*/
PortOut::PortOut(){};

/*!Default destructor of PortOut
*/
PortOut::~PortOut(){};

/*!Copy constructor of PortOut.
*/
PortOut::PortOut(const PortOut & other){
	m_objLink 	= other.m_objLink;
	m_obuffer	= other.m_obuffer;
	m_portLink	= other.m_portLink;
	return;
};

/*!Assignement operator of PortOut.
 */
PortOut & PortOut::operator=(const PortOut & other){
	m_objLink 	= other.m_objLink;
	m_obuffer	= other.m_obuffer;
	m_portLink	= other.m_portLink;
	return (*this);
};

/*!Compare operator of PortOut.
 */
bool PortOut::operator==(const PortOut & other){
	bool check = true;
	check = check && (m_portLink == other.m_portLink);
	check = check && (m_objLink == other.m_objLink);
	return(check);
};

/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
std::vector<mimmo::BaseManipulation*>
PortOut::getLink(){
	return(m_objLink);
}

std::vector<int>
PortOut::getPortLink(){
	return(m_portLink);
}

void
mimmo::PortOut::cleanBuffer(){
	m_obuffer.setCapacity(0);
	m_obuffer.eof();
}

void
mimmo::PortOut::clear(){
	m_objLink.clear();
	m_portLink.clear();
}

void
mimmo::PortOut::clear(int j){
	if (j < m_objLink.size() && j >= 0){
		m_objLink.erase(m_objLink.begin() + j);
		m_portLink.erase(m_portLink.begin() + j);
	}
}

/*! Execution of the pin.
 * All the pins of an object are called in execute of the owner after its own execution.
 */
void
mimmo::PortOut::exec(){
	if (m_objLink.size() > 0){
		writeBuffer();
		bitpit::IBinaryStream input(m_obuffer.data());
		cleanBuffer();
		for (int j=0; j<m_objLink.size(); j++){
			if (m_objLink[j] != NULL){
				m_objLink[j]->setBufferIn(m_portLink[j], input);
				m_objLink[j]->readBufferIn(m_portLink[j]);
				m_objLink[j]->cleanBufferIn(m_portLink[j]);
			}
		}
	}
};

/*!Default constructor of PortIn
*/
PortIn::PortIn(){};

/*!Default destructor of PortIn
*/
PortIn::~PortIn(){
	m_objLink	= NULL;
};

/*!Copy constructor of PortIn.
*/
PortIn::PortIn(const PortIn & other){
	m_objLink 	= other.m_objLink;
	m_ibuffer	= other.m_ibuffer;
	return;
};

/*!Assignement operator of PortIn.
 */
PortIn & PortIn::operator=(const PortIn & other){
	m_objLink 	= other.m_objLink;
	m_ibuffer	= other.m_ibuffer;
	return (*this);
};

/*!Compare operator of PortIn.
 */
bool PortIn::operator==(const PortIn & other){
	bool check = true;
	check = check && (m_objLink == other.m_objLink);
	return(check);
};

void
PortIn::addCompatibility(PortType label){
	m_labelOK.push_back(label);
}

const std::vector<PortType>&
PortIn::getCompatibility(){
	return m_labelOK;
}


/*!It gets the linked object by this pin.
 * \return Pointer to linked object.
 */
BaseManipulation*
PortIn::getLink(){
	return(m_objLink);
}

void
mimmo::PortIn::clear(){
	m_objLink = NULL;
}

void
mimmo::PortIn::cleanBuffer(){
	m_ibuffer.setCapacity(0);
	m_ibuffer.eof();
}

