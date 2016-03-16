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

///*!Default constructor of InOut
// */
InOut::InOut(){};
//
///*!Default destructor of InOut
// */
InOut::~InOut(){
	m_obj 		= NULL;
	m_objIn		= NULL;
	m_objOut	= NULL;
	m_getVal 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOut.
// */
InOut::InOut(const InOut & other){
	m_obj 		= other.m_obj;
	m_objIn 	= other.m_objIn;
	m_objOut 	= other.m_objOut;
	m_getVal 	= other.m_getVal;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOut.
 */
InOut & InOut::operator=(const InOut & other){
	m_obj 		= other.m_obj;
	m_objIn 	= other.m_objIn;
	m_objOut 	= other.m_objOut;
	m_getVal 	= other.m_getVal;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOut::setInput(BaseManipulation* objIn, function<dvecarr3E*(void)> getVal, function<void(dvecarr3E&)> setVal){
	m_objOut	= NULL;
	m_objIn		= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOut::setOutput(BaseManipulation* objOut, function<void(dvecarr3E&)> setVal, function<dvecarr3E*(void)> getVal){
	m_objOut	= objOut;
	m_objIn		= NULL;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

