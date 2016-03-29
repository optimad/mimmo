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
 *	\brief InOut is the input-output PIN base class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 */
class InOut{
public:
	//members
	BaseManipulation*	m_objLink;	/**<Input/Output object from/to which
										recover/give the target variable. */

public:
	InOut();
	~InOut();

	InOut(const InOut & other);
	InOut & operator=(const InOut & other);
	bool operator==(const InOut & other);

	BaseManipulation*	getLink();

	virtual void exec() = 0;

};

//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief InOutT is the input-output templated PIN class derived from base InOut pin class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 *	The parent object gives the output value by a function set by the user
 *	during the creation of the pin (getVal functions), while the child object
 *	use this value by a set function (setVal function).
 *	The functions are given to the pin by get/set function objects with any
 *	type of value (the same for the two functions of the linked object) by using
 *	the template formulation. The functions have to be function objects of the standard library
 *	(functional include) created by the bind method.
 *	The input/output value can be returned by copy, reference or pointer by the get function
 *	of the parent and it can be passed by copy or pointer to the set value of the child object.
 *	The output pins of an object are executed after its own execution.
 *
 */
template<typename T>
class InOutT: public InOut {

public:
	std::function<T(void)>	m_getVal;	/**<Pointer to get function with copy return.*/
	std::function<T&(void)>	m_getValR;	/**<Pointer to get function with reference return.*/
	std::function<T*(void)>	m_getValP;	/**<Pointer to get function with pointer return.*/
	std::function<void(T)>	m_setVal;	/**<Pointer to set function with copy argument.*/
	std::function<void(T*)>	m_setValP;	/**<Pointer to set function with pointer argument.*/

public:
	InOutT();
	~InOutT();

	InOutT(const InOutT & other);
	InOutT & operator=(const InOutT & other);
	bool operator==(const InOutT & other);

	void setInput(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<T&(void)> getValR, std::function<void(T)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getValR);
	void setInput(BaseManipulation* objIn, std::function<T*(void)> getValP, std::function<void(T)> setVal);
	void setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getValP);

	void setInput(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setValP);
	void setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T(void)> getVal);
	void setInput(BaseManipulation* objIn, std::function<T&(void)> getValR, std::function<void(T*)> setValP);
	void setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T&(void)> getValR);
	void setInput(BaseManipulation* objIn, std::function<T*(void)> getValP, std::function<void(T*)> setValP);
	void setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T*(void)> getValP);

	void exec();

};

#include "InOut.tpp"

#endif /* __INOUT_HPP__ */
