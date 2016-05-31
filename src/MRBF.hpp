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
 * @enum MRBFSol
 * @brief Solver enum for your RBF data fields interpolation/ direct parameterization
 */
enum class MRBFSol{
	NONE = 0, 	/**< activate class as pure parameterizator. Set freely your RBF coefficients/weights */
	WHOLE = 1,	/**< activate class as pure interpolator, with RBF coefficients evaluated solving a full linear system for all active nodes.*/ 
	GREEDY= 2   /**< activate class as pure interpolator, with RBF coefficients evaluated using a greedy algorithm on active nodes.*/
};
	
/*!
 *	\date			25/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *  @class MRBF
 *	\brief Radial Basis Function evaluation from clouds of control points.
 *
 *	This class is derived from BaseManipulation class of MiMMO and from RBF class
 *	of bitpit library.
 *	It evaluates the result of RBF built over a set of control point given by the user
 *	or stored in a MimmoObject (geometry container). Default solver in execution is
 *	MRBFSol::NONE for direct parameterization. Use MRBFSol::GREEDY or MRBFSol::SOLVE to activate
 *  interpolation features.
 *  See bitpit::RBF docs for further information.
 *
 * 
 */

//TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF
//TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF

class MRBF: public BaseManipulation, public bitpit::RBF {

private:
	double 		m_tol;		/**< Tolerance for greedy algorithm.*/
	MRBFSol		m_solver; 	/**<Type of solver specified for the class as default in execution*/
	dvector1D	m_filter;	/**<Filter field for displacements modulation */
	bool		m_bfilter;		/**<boolean to recognize if a filter field is applied */
	
	bool		m_srset;	/**<True if the support radius is set by the user.*/

	dvecarr3E	m_displ;

public:
	MRBF();
	virtual ~MRBF();

	//copy operators/constructors
	MRBF(const MRBF & other);
	MRBF & operator=(const MRBF & other);

	void 			setGeometry(MimmoObject* geometry);
	
	dvecarr3E*		getNodes();

	MRBFSol			getMode();
	void			setMode(MRBFSol);
	void			setMode(int);
	dvector1D		getFilter();
	
	std::pair<MimmoObject * , dvecarr3E * >	getDeformedField();
	
	int 			addNode(darray3E);
	ivector1D		addNode(dvecarr3E);
	ivector1D	 	addNode(MimmoObject* geometry);

	void 			setNode(darray3E);
	void			setNode(dvecarr3E);
	void		 	setNode(MimmoObject* geometry);
	void			setFilter(dvector1D );
	
	ivector1D		checkDuplicatedNodes(double tol=1.0E-12);
	bool 			removeDuplicatedNodes(ivector1D * list=NULL);
	
	void 			setSupportRadius(const double & suppR);
	void 			setTol(double tol);
	void 			setDisplacements(dvecarr3E displ);

	void 		clear();
	void 		clearFilter();
	
	//execute deformation methods
	void 			execute();

private:
	void			setWeight(dvector2D value);

};

}

#endif /* __MRBF_HPP__ */
