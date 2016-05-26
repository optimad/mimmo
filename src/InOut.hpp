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
#include "binary_stream.hpp"
#include <functional>

namespace mimmo{

class BaseManipulation;

/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief PinOut is the input-output PIN base class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 */
class PinOut{
public:
	//members
	bitpit::OBinaryStream				m_obuffer;
	std::vector<BaseManipulation*>		m_objLink;	/**<Input/Output object from/to which
										recover/give the target variable. */
	std::vector<int>					m_portLink;


public:
	PinOut();
	virtual ~PinOut();

	PinOut(const PinOut & other);
	PinOut & operator=(const PinOut & other);
	bool operator==(const PinOut & other);

	std::vector<BaseManipulation*>	getLink();
	std::vector<int>				getPortLink();

	virtual void	writeBuffer() = 0;

	/*!Execution method.
	 */
	void exec();

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief PinOutT is the input-output templated PIN class derived from base PinOut pin class.
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
class PinOutT: public PinOut {

public:

	T		*m_var_;

public:
	PinOutT();
	PinOutT(T *var_);
	virtual ~PinOutT();

	PinOutT(const PinOutT & other);
	PinOutT & operator=(const PinOutT & other);
	bool operator==(const PinOutT & other);

	void writeBuffer();
	void readBuffer();

};



/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief PinOut is the input-output PIN base class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 */
class PinIn{
public:
	//members
	bitpit::IBinaryStream				m_ibuffer;
	BaseManipulation*					m_objLink;	/**<Input/Output object from/to which
										recover/give the target variable. */

public:
	PinIn();
	virtual ~PinIn();

	PinIn(const PinIn & other);
	PinIn & operator=(const PinIn & other);
	bool operator==(const PinIn & other);

	BaseManipulation*	getLink();

	virtual void	readBuffer() = 0;

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
 *	\date			14/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief PinInT is the input-output templated PIN class derived from base PinIn pin class.
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
class PinInT: public PinIn {

public:

	T		*m_var_;

public:
	PinInT();
	PinInT(T *var_);
	virtual ~PinInT();

	PinInT(const PinInT & other);
	PinInT & operator=(const PinInT & other);
	bool operator==(const PinInT & other);

	void writeBuffer();
	void readBuffer();

};


}

#include "InOut.tpp"

#endif /* __INOUT_HPP__ */
