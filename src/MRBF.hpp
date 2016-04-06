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
 *	or stored in a MimmoObject (geometry container). Class is built as default in 
 *  bitpit::RBFType::PARAM mode. See bitpit::RBF docs for further information. 
 * 
 *  \\TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF
 */
class MRBF: public BaseManipulation, public bitpit::RBF {
private:
	double m_tol;		/**< Tolerance for greedy algorithm.*/

public:
	MRBF();
	virtual ~MRBF();

	//copy operators/constructors
	MRBF(const MRBF & other);
	MRBF & operator=(const MRBF & other);

	int 			addNode(darray3E &);
	ivector1D		addNode(dvecarr3E &);
	ivector1D	 	addNode(MimmoObject* geometry);
	
	ivector1D		checkDuplicatedNodes(double tol=1.0E-12);
	bool 			removeDuplicatedNodes(ivector1D * list=NULL);
	
	void 			setTol(double tol);
	void 			setDisplacements(dvecarr3E);
	void			setActiveDisplacements(dvecarr3E);
	//execute deformation methods
	void 		execute();
	
private:
	dvector1D convertActiveToTotal(dvector1D &);
};

}

#endif /* __MRBF_HPP__ */
