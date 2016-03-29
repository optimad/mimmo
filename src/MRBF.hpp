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
#ifndef __MRBF_HPP__
#define __MRBF_HPP__

#include "BaseManipulation.hpp"
#include "rbf.hpp"

namespace mimmo{

/*!
 *	\date			25/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Radial Basis Function evaluation from clouds of control points.
 *
 *	This class is derived from BaseManipulation class of MiMMO and from RBF class
 *	of bitpit library.
 *	It evaluates the result of RBF built over a set of control point given by the user
 *	or stored in a MimmoObject (geometry container).
 *
 */
class MRBF: public BaseManipulation, public bitpit::RBF {

protected:

public:

private:

	double m_tol; 	/**< Tolerance value for greedy method in node insert.*/


public:
	MRBF();
	virtual ~MRBF();

	//copy operators/constructors
	MRBF(const MRBF & other);
	MRBF & operator=(const MRBF & other);

	//internal methods
	void setTol(double tol);
	void addNodes(dvecarr3E nodes);
	void addNodes(MimmoObject* geometry);
	void addField(dvecarr3E field);

	//execute deformation methods
	void 		execute();

};

}

#endif /* __MRBF_HPP__ */
