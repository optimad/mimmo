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

//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION								//
//==============================================================//

///*!Default constructor of InOut
// */
InOut::InOut(){};
//
///*!Default destructor of InOut
// */
InOut::~InOut(){
	m_type 		= false;
	m_obj 		= NULL;
	m_objLink	= NULL;
};

///*!Copy constructor of InOut.
// */
InOut::InOut(const InOut & other){
	m_type 		= other.m_type;
	m_obj 		= other.m_obj;
	m_objLink 	= other.m_objLink;
};

/*!Assignement operator of InOut.
 */
InOut & InOut::operator=(const InOut & other){
	m_type 		= other.m_type;
	m_obj 		= other.m_obj;
	m_objLink 	= other.m_objLink;
	return (*this);
};

bool
InOut::getType(){
	return(m_type);
}

BaseManipulation*
InOut::getLink(){
	return(m_objLink);
}


//==============================================================//
// DOUBLE DERIVED INOUT CLASS									//
//==============================================================//

///*!Default constructor of InOutDVA3
// */
InOutD::InOutD(){
	m_getValR 	= NULL;
	m_getVal 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutDVA3
// */
InOutD::~InOutD(){
	m_getValR 	= NULL;
	m_getVal 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutD.
// */
InOutD::InOutD(const InOutD & other){
	m_getValR 	= other.m_getValR;
	m_getVal 	= other.m_getVal;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutD.
 */
InOutD & InOutD::operator=(const InOutD & other){
	m_getValR 	= other.m_getValR;
	m_getVal 	= other.m_getVal;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutD::setInput(BaseManipulation* objIn, function<double(void)> getVal, function<void(double)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutD::setOutput(BaseManipulation* objOut, function<void(double)> setVal, function<double(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutD::setInput(BaseManipulation* objIn, function<double&(void)> getValR, function<void(double)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutD::setOutput(BaseManipulation* objOut, function<void(double)> setVal, function<double&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutD::exec(){
	if (m_getVal != NULL){
		double val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		double& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// INT DERIVED INOUT CLASS										//
//==============================================================//

///*!Default constructor of InOutIVA3
// */
InOutI::InOutI(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutIVA3
// */
InOutI::~InOutI(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutI.
// */
InOutI::InOutI(const InOutI & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutI.
 */
InOutI & InOutI::operator=(const InOutI & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutI::setInput(BaseManipulation* objIn, function<int(void)> getVal, function<void(int)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutI::setOutput(BaseManipulation* objOut, function<void(int)> setVal, function<int(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutI::setInput(BaseManipulation* objIn, function<int&(void)> getValR, function<void(int)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutI::setOutput(BaseManipulation* objOut, function<void(int)> setVal, function<int&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutI::exec(){
	if (m_getVal != NULL){
		int val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		int& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// VECTOR OF DOUBLE DERIVED INOUT CLASS							//
//==============================================================//

///*!Default constructor of InOutDV1VA3
// */
InOutDV1::InOutDV1(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutDV1VA3
// */
InOutDV1::~InOutDV1(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutDV1.
// */
InOutDV1::InOutDV1(const InOutDV1 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutDV1.
 */
InOutDV1 & InOutDV1::operator=(const InOutDV1 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutDV1::setInput(BaseManipulation* objIn, function<dvector1D(void)> getVal, function<void(dvector1D)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDV1::setOutput(BaseManipulation* objOut, function<void(dvector1D)> setVal, function<dvector1D(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDV1::setInput(BaseManipulation* objIn, function<dvector1D&(void)> getValR, function<void(dvector1D)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDV1::setOutput(BaseManipulation* objOut, function<void(dvector1D)> setVal, function<dvector1D&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDV1::exec(){
	if (m_getVal != NULL){
		dvector1D val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		dvector1D& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// VECTOR OF INT DERIVED INOUT CLASS							//
//==============================================================//

///*!Default constructor of InOutIV1VA3
// */
InOutIV1::InOutIV1(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutIV1VA3
// */
InOutIV1::~InOutIV1(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutIV1.
// */
InOutIV1::InOutIV1(const InOutIV1 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutIV1.
 */
InOutIV1 & InOutIV1::operator=(const InOutIV1 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutIV1::setInput(BaseManipulation* objIn, function<ivector1D(void)> getVal, function<void(ivector1D)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIV1::setOutput(BaseManipulation* objOut, function<void(ivector1D)> setVal, function<ivector1D(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIV1::setInput(BaseManipulation* objIn, function<ivector1D&(void)> getValR, function<void(ivector1D)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutIV1::setOutput(BaseManipulation* objOut, function<void(ivector1D)> setVal, function<ivector1D&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};


void
InOutIV1::exec(){
	if (m_getVal != NULL){
		ivector1D val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		ivector1D& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// ARRAY OF DOUBLE OF SIZE 3 DERIVED INOUT CLASS				//
//==============================================================//

///*!Default constructor of InOutDA3VA3
// */
InOutDA3::InOutDA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutDA3VA3
// */
InOutDA3::~InOutDA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutDA3.
// */
InOutDA3::InOutDA3(const InOutDA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutDA3.
 */
InOutDA3 & InOutDA3::operator=(const InOutDA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutDA3::setInput(BaseManipulation* objIn, function<darray3E(void)> getVal, function<void(darray3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDA3::setOutput(BaseManipulation* objOut, function<void(darray3E)> setVal, function<darray3E(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDA3::setInput(BaseManipulation* objIn, function<darray3E&(void)> getValR, function<void(darray3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDA3::setOutput(BaseManipulation* objOut, function<void(darray3E)> setVal, function<darray3E&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDA3::exec(){
	if (m_getVal != NULL){
		darray3E val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		darray3E& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// ARRAY OF INT OF SIZE 3 DERIVED INOUT CLASS					//
//==============================================================//

///*!Default constructor of InOutIA3VA3
// */
InOutIA3::InOutIA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutIA3VA3
// */
InOutIA3::~InOutIA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutIA3.
// */
InOutIA3::InOutIA3(const InOutIA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutIA3.
 */
InOutIA3 & InOutIA3::operator=(const InOutIA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutIA3::setInput(BaseManipulation* objIn, function<iarray3E(void)> getVal, function<void(iarray3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIA3::setOutput(BaseManipulation* objOut, function<void(iarray3E)> setVal, function<iarray3E(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIA3::setInput(BaseManipulation* objIn, function<iarray3E&(void)> getValR, function<void(iarray3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutIA3::setOutput(BaseManipulation* objOut, function<void(iarray3E)> setVal, function<iarray3E&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutIA3::exec(){
	if (m_getVal != NULL){
		iarray3E val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		iarray3E& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// VECTOR OF ARRAY OF SIZE 3 OF INT DERIVED INOUT CLASS			//
//==============================================================//

///*!Default constructor of InOutIVA3
// */
InOutIVA3::InOutIVA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutIVA3
// */
InOutIVA3::~InOutIVA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutIVA3.
// */
InOutIVA3::InOutIVA3(const InOutIVA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutIVA3.
 */
InOutIVA3 & InOutIVA3::operator=(const InOutIVA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutIVA3::setInput(BaseManipulation* objIn, function<ivecarr3E(void)> getVal, function<void(ivecarr3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIVA3::setOutput(BaseManipulation* objOut, function<void(ivecarr3E)> setVal, function<ivecarr3E(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutIVA3::setInput(BaseManipulation* objIn, function<ivecarr3E&(void)> getValR, function<void(ivecarr3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutIVA3::setOutput(BaseManipulation* objOut, function<void(ivecarr3E)> setVal, function<ivecarr3E&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutIVA3::exec(){
	if (m_getVal != NULL){
		ivecarr3E val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		ivecarr3E& val = m_getValR();
		m_setVal(val);
	}
};


//==============================================================//
// VECTOR OF ARRAY OF SIZE 3 OF DOUBLE DERIVED INOUT CLASS		//
//==============================================================//

///*!Default constructor of InOutDVA3
// */
InOutDVA3::InOutDVA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

//
///*!Default destructor of InOutDVA3
// */
InOutDVA3::~InOutDVA3(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_setVal 	= NULL;
};

///*!Copy constructor of InOutDVA3.
// */
InOutDVA3::InOutDVA3(const InOutDVA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
};

/*!Assignement operator of InOutDVA3.
 */
InOutDVA3 & InOutDVA3::operator=(const InOutDVA3 & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_setVal 	= other.m_setVal;
	return (*this);
};

void
InOutDVA3::setInput(BaseManipulation* objIn, function<dvecarr3E(void)> getVal, function<void(dvecarr3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDVA3::setOutput(BaseManipulation* objOut, function<void(dvecarr3E)> setVal, function<dvecarr3E(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

void
InOutDVA3::setInput(BaseManipulation* objIn, function<dvecarr3E&(void)> getValR, function<void(dvecarr3E)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDVA3::setOutput(BaseManipulation* objOut, function<void(dvecarr3E)> setVal, function<dvecarr3E&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

void
InOutDVA3::exec(){
	if (m_getVal != NULL){
		dvecarr3E val = m_getVal();
		m_setVal(val);
	}else if (m_getValR != NULL){
		dvecarr3E& val = m_getValR();
		m_setVal(val);
	}
};






