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
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------------|
 *	|                  Port Input                                   |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType       		|
 *	|-------|-----------|-------------------|-----------------------|
 *	| 0     | M_COORDS  | setNode           | (VECARR3, FLOAT)		|
 *	| 10    | M_DISPLS  | setDisplacements  | (VECARR3, FLOAT)		|
 *	| 12    | M_FILTER  | setFilter         | (VECTOR, FLOAT) 		|
 *	| 30    | M_VALUED  | setSupportRadius  | (SCALAR, FLOAT)	  	|
 *	| 99    | M_GEOM    | m_geometry        | (SCALAR, MIMMO_)  		|
 *	| 130   | M_VALUED2 | setTol            | (SCALAR, FLOAT)		|
 *	|-------|-----------|-------------------|-----------------------|
 *
 *
 *	|--------------------------------------------|------------------------------|
 *	|             Port Output	                 |                      		|
 *	|-------|----------------|-------------------|------------------------------|
 *	|PortID | PortType  	 | variable/function | DataType    		   			|
 *	|-------|----------------|-------------------|------------------------------|
 *	| 11    | M_GDISPLS 	 | getDisplacements  | (VECARR3, FLOAT)				|
 *	| 80    | M_PAIRVECFIELD | getDeformedField  | (PAIR, MIMMO_VECARR3FLOAT_)	|
 *	|-------|----------------|-------------------|------------------------------|
 * ~~~
 *	=========================================================
 *
 */
//TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF
class MRBF: public BaseManipulation, public bitpit::RBF {

private:
	double 		m_tol;		/**< Tolerance for greedy algorithm.*/
	MRBFSol		m_solver; 	/**<Type of solver specified for the class as default in execution*/
	dvector1D	m_filter;	/**<Filter field for displacements modulation */
	bool		m_bfilter;	/**<boolean to recognize if a filter field is applied */
	double		m_SRRatio;	/**<support Radius ratio */

	dvecarr3E	m_displ;	/**<Rsulting displacements of geometry vertex.*/

public:
	MRBF();
	virtual ~MRBF();

	//copy operators/constructors
	MRBF(const MRBF & other);
	MRBF & operator=(const MRBF & other);

	void buildPorts();

	void 			setGeometry(MimmoObject* geometry);
	
	dvecarr3E*		getNodes();

	MRBFSol			getMode();
	void			setMode(MRBFSol);
	void			setMode(int);
	dvector1D		getFilter();
	double			getSupportRadius();
	double			getSupportRadiusValue();
	
	std::pair<MimmoObject * , dvecarr3E * >	getDeformedField();
	dvecarr3E		getDisplacements();
	
	int 			addNode(darray3E);
	ivector1D		addNode(dvecarr3E);
	ivector1D	 	addNode(MimmoObject* geometry);

	void 			setNode(darray3E);
	void			setNode(dvecarr3E);
	void		 	setNode(MimmoObject* geometry);
	void			setFilter(dvector1D );
	
	ivector1D		checkDuplicatedNodes(double tol=1.0E-12);
	bool 			removeDuplicatedNodes(ivector1D * list=NULL);
	
	void 			setSupportRadius(double suppR_);
	void 			setTol(double tol);
	void 			setDisplacements(dvecarr3E displ);

	void 		clear();
	void 		clearFilter();
	
	//execute deformation methods
	void 			execute();

	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	void			setWeight(dvector2D value);

};

}

#endif /* __MRBF_HPP__ */
