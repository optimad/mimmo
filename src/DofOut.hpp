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
#ifndef __DOFOUT_HPP__
#define __DOFOUT_HPP__

#include "MiMMO_TypeDef.hpp"
#include <functional>

namespace mimmo{

class BaseManipulation;

/*!
 *	\date			31/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief DofOut is the input-output base class for chain.
 *
 *
 */
class DofOut{

	//friendship declaration
	friend class Chain;

private:
	//members
	BaseManipulation*	m_obj;		/**<InputDOF/Output object from/to which
										recover/give the target variables.*/
	int					m_nglob;	/**<Global number of variables given by linked object.*/
	int					m_nuse;		/**<Number of variables to be used as dof/output.*/
	std::vector<bool>	m_actives;	/**<List of true/false to map the variables to be used as dof/output.*/

public:
	DofOut();
	~DofOut();

	DofOut(const DofOut & other);
	DofOut & operator=(const DofOut & other);
	bool operator==(const DofOut & other);

private:
	BaseManipulation*	getLink();
	int					getNuse();

};

//==============================================================//
// TEMPLATE DERIVED DOFOUT CLASS									//
//==============================================================//

/*!
 *	\date			31/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief DofOutT is the input-output templated class derived from base DofOut class.
 *
 *
 */
template<typename T, typename C>
class DofOutT: public DofOut {
public:

	std::function<T(void)>	m_getVal;	/**<Pointer to get function with copy return.*/
	std::function<T*(void)>	m_getValP;	/**<Pointer to get function with pointer return.*/

	DofOutT();
	DofOutT(BaseManipulation* obj, T (C::*fget) ());
	DofOutT(BaseManipulation* obj, T* (C::*fget) ());
	~DofOutT();

	DofOutT(const DofOutT & other);
	DofOutT & operator=(const DofOutT & other);
	bool operator==(const DofOutT & other);

	void setGetFunction(T (C::*fget) ());
	void setGetFunction(T* (C::*fget) ());

	T get();
	T* getP();

};


}

#include "DofOut.tpp"

#endif /* __DOFOUT_HPP__ */
