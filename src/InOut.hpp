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
#ifndef __INOUT_HPP__
#define __INOUT_HPP__

#include "MiMMO_TypeDef.hpp"
#include <functional>

class BaseManipulation;
/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief InOut is the input-output PIN class.
 *
 *	OR input OR output.
 *
 */
class InOut{
public:
	//members
	bool				m_type;		/**<Type of pin 0/1 = input/output pin.*/
	BaseManipulation*	m_obj;		/**<Owner object of this input-output PIN. */
	BaseManipulation*	m_objLink;	/**<Input/Output object from/to which recover/give the target variable. */


public:
	InOut();
	~InOut();

	InOut(const InOut & other);
	InOut & operator=(const InOut & other);

	bool				getType();
	BaseManipulation*	getLink();

	virtual void exec() = 0;

};

//==============================================================//
// DOUBLE DERIVED INOUT CLASS									//
//==============================================================//

class InOutD: public InOut {

public:
	std::function<double(void)>		m_getVal;
	std::function<double&(void)>	m_getValR;
	std::function<void(double)>		m_setVal;

public:
	InOutD();
	~InOutD();

	InOutD(const InOutD & other);
	InOutD & operator=(const InOutD & other);

	void setInput(BaseManipulation* objIn, std::function<double(void)> getVal, std::function<void(double)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(double)> setVal, std::function<double(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<double&(void)> getValR, std::function<void(double)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(double)> setVal, std::function<double&(void)> getValR);

	void exec();

};


//==============================================================//
// INT DERIVED INOUT CLASS										//
//==============================================================//

class InOutI: public InOut {

public:
	std::function<int(void)>	m_getVal;
	std::function<int&(void)>	m_getValR;
	std::function<void(int)>	m_setVal;

public:
	InOutI();
	~InOutI();

	InOutI(const InOutI & other);
	InOutI & operator=(const InOutI & other);

	void setInput(BaseManipulation* objIn, std::function<int(void)> getVal, std::function<void(int)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(int)> setVal, std::function<int(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<int&(void)> getVal, std::function<void(int)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(int)> setVal, std::function<int&(void)> getVal);

	void exec();

};


//==============================================================//
// VECTOR OF DOUBLE DERIVED INOUT CLASS							//
//==============================================================//

class InOutDV1: public InOut {

public:
	std::function<dvector1D(void)>	m_getVal;
	std::function<dvector1D&(void)>	m_getValR;
	std::function<void(dvector1D)>	m_setVal;

public:
	InOutDV1();
	~InOutDV1();

	InOutDV1(const InOutDV1 & other);
	InOutDV1 & operator=(const InOutDV1 & other);

	void setInput(BaseManipulation* objIn, std::function<dvector1D(void)> getVal, std::function<void(dvector1D)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(dvector1D)> setVal, std::function<dvector1D(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<dvector1D&(void)> getVal, std::function<void(dvector1D)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(dvector1D)> setVal, std::function<dvector1D&(void)> getVal);

	void exec();

};


//==============================================================//
// VECTOR OF INT DERIVED INOUT CLASS							//
//==============================================================//

class InOutIV1: public InOut {

public:
	std::function<ivector1D(void)>	m_getVal;
	std::function<ivector1D&(void)>	m_getValR;
	std::function<void(ivector1D)>	m_setVal;

public:
	InOutIV1();
	~InOutIV1();

	InOutIV1(const InOutIV1 & other);
	InOutIV1 & operator=(const InOutIV1 & other);

	void setInput(BaseManipulation* objIn, std::function<ivector1D(void)> getVal, std::function<void(ivector1D)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(ivector1D)> setVal, std::function<ivector1D(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<ivector1D&(void)> getVal, std::function<void(ivector1D)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(ivector1D)> setVal, std::function<ivector1D&(void)> getVal);

	void exec();

};


//==============================================================//
// ARRAY OF SIZE 3 OF DOUBLE DERIVED INOUT CLASS				//
//==============================================================//

class InOutDA3: public InOut {

public:
	std::function<darray3E(void)>	m_getVal;
	std::function<darray3E&(void)>	m_getValR;
	std::function<void(darray3E)>	m_setVal;

public:
	InOutDA3();
	~InOutDA3();

	InOutDA3(const InOutDA3 & other);
	InOutDA3 & operator=(const InOutDA3 & other);

	void setInput(BaseManipulation* objIn, std::function<darray3E(void)> getVal, std::function<void(darray3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(darray3E)> setVal, std::function<darray3E(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<darray3E&(void)> getVal, std::function<void(darray3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(darray3E)> setVal, std::function<darray3E&(void)> getVal);

	void exec();

};


//==============================================================//
// ARRAY OF SIZE 3 OF INT DERIVED INOUT CLASS					//
//==============================================================//

class InOutIA3: public InOut {

public:
	std::function<iarray3E(void)>	m_getVal;
	std::function<iarray3E&(void)>	m_getValR;
	std::function<void(iarray3E)>	m_setVal;

public:
	InOutIA3();
	~InOutIA3();

	InOutIA3(const InOutIA3 & other);
	InOutIA3 & operator=(const InOutIA3 & other);

	void setInput(BaseManipulation* objIn, std::function<iarray3E(void)> getVal, std::function<void(iarray3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(iarray3E)> setVal, std::function<iarray3E(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<iarray3E&(void)> getVal, std::function<void(iarray3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(iarray3E)> setVal, std::function<iarray3E&(void)> getVal);

	void exec();

};

//==============================================================//
// VECTOR OF ARRAY OF SIZE 3 OF DOUBLE DERIVED INOUT CLASS		//
//==============================================================//

class InOutDVA3: public InOut {

public:
	std::function<dvecarr3E(void)>	m_getVal;
	std::function<dvecarr3E&(void)>	m_getValR;
	std::function<void(dvecarr3E)>	m_setVal;

public:
	InOutDVA3();
	~InOutDVA3();

	InOutDVA3(const InOutDVA3 & other);
	InOutDVA3 & operator=(const InOutDVA3 & other);

	void setInput(BaseManipulation* objIn, std::function<dvecarr3E(void)> getVal, std::function<void(dvecarr3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(dvecarr3E)> setVal, std::function<dvecarr3E(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<dvecarr3E&(void)> getVal, std::function<void(dvecarr3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(dvecarr3E)> setVal, std::function<dvecarr3E&(void)> getVal);

	void exec();

};

//==============================================================//
// VECTOR OF ARRAY OF SIZE 3 OF INT DERIVED INOUT CLASS			//
//==============================================================//

class InOutIVA3: public InOut {

public:
	std::function<ivecarr3E(void)>	m_getVal;
	std::function<ivecarr3E&(void)>	m_getValR;
	std::function<void(ivecarr3E)>	m_setVal;

public:
	InOutIVA3();
	~InOutIVA3();

	InOutIVA3(const InOutIVA3 & other);
	InOutIVA3 & operator=(const InOutIVA3 & other);

	void setInput(BaseManipulation* objIn, std::function<ivecarr3E(void)> getVal, std::function<void(ivecarr3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(ivecarr3E)> setVal, std::function<ivecarr3E(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<ivecarr3E&(void)> getValR, std::function<void(ivecarr3E)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(ivecarr3E)> setVal, std::function<ivecarr3E&(void)> getValR);

	void exec();

};

#endif /* __INOUT_HPP__ */
