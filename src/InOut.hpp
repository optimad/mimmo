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
	BaseManipulation*	m_obj;		/**<Owner object of this input-output PIN. */
	BaseManipulation*	m_objIn;	/**<Input object from which recover the target variable. */
	BaseManipulation*	m_objOut;	/**<Output object to which give the target variable. */

	std::function<dvecarr3E*(void)>	m_getVal;
	std::function<void(dvecarr3E&)>	m_setVal;


public:
	InOut();
	~InOut();

	InOut(const InOut & other);
	InOut & operator=(const InOut & other);

	void	setInput(BaseManipulation* objIn, std::function<dvecarr3E*(void)> getVal, std::function<void(dvecarr3E&)> setVal);
	void	setOutput(BaseManipulation* objOut, std::function<void(dvecarr3E&)> setVal, std::function<dvecarr3E*(void)> getVal);


private:


};

#endif /* __INOUT_HPP__ */
