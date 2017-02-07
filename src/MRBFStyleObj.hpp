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
#ifndef __MRBFSTYLEOBJ_HPP__
#define __MRBFSTYLEOBJ_HPP__

#include "MRBF.hpp"

namespace mimmo{
	
/*!
 *	\date			25/mar/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *  @class MRBFStyleObj
 *	\brief Radial Basis Function evaluation from discrete meshes as control nodes .
 *
 *	This class is derived from BaseManipulation class of MiMMO and from RBFAbstract class
 *	of bitpit library.
 *	It evaluates the result of RBF built over a set of control node geometries (3D curves or surfaces) given by the User
 *	stored in a MimmoObject (geometry container). Default solver in execution is
 *	MRBFSol::NONE for direct parameterization. Use MRBFSol::GREEDY or MRBFSol::SOLVE to activate
 *  interpolation features.
 *  See bitpit::RBFAbstract docs for further information.
 *
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------|
 *	|                  Port Input                                       |
 *	|-------|-----------|-----------------------|-----------------------|
 *	|PortID | PortType  | variable/function     | DataType       		|
 *	|-------|-----------|-----------------------|-----------------------|
 *  | 98    | M_GEOM2   | setAddNode            | (SCALAR, MIMMO_)		|
 *	| 100   | M_VECGEOM | setAddNode            | (VECTOR, MIMMO_)		|
 *	| 10    | M_DISPLS  | setDisplacements      | (VECARR3, FLOAT)		|
 *	| 12    | M_FILTER  | setFilter             | (VECTOR, FLOAT) 		|
 *  | 30    | M_VALUED  | setSupportRadius      | (SCALAR, FLOAT)       |
 *  | 130   | M_VALUED2 | setSupportRadiusValue | (SCALAR, FLOAT)       |
 *	| 99    | M_GEOM    | m_geometry            | (SCALAR, MIMMO_)  	|
 *	|-------|-----------|-----------------------|-----------------------|
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
class MRBFStyleObj: public BaseManipulation, public bitpit::RBFAbstract {

private:
	double 		m_tol;	    	/**< Tolerance for greedy algorithm.*/
	MRBFSol		m_solver;   	/**<Type of solver specified for the class as default in execution*/
	dvector1D	m_filter;   	/**<Filter field for displacements modulation */
	bool		m_bfilter;  	/**<boolean to recognize if a filter field is applied */
	double		m_SRRatio;	    /**<support Radius ratio */
	dvecarr3E	m_displ;	    /**<Resulting displacements of geometry vertex.*/
	bool        m_supRIsValue;  /**<True if support radius is defined as absolute value, false if is ratio of bounding box diagonal.*/

	std::vector<MimmoObject * > m_node; /**< list of nodes for your RBF parameterization */
public:
	MRBFStyleObj();
	virtual ~MRBFStyleObj();

	//copy operators/constructors
	MRBFStyleObj(const MRBFStyleObj & other);
	MRBFStyleObj & operator=(const MRBFStyleObj & other);

	void buildPorts();

	void 			setGeometry(MimmoObject* geometry);
	
	std::vector<MimmoObject*>		getNodes();
	int 							getTotalNodesCount();
	
	MRBFSol			getMode();
	void			setMode(MRBFSol);
	void			setMode(int);
	dvector1D		getFilter();
	double			getSupportRadius();
	double			getSupportRadiusValue();
    bool            getIsSupportRadiusValue();
	
	std::pair<MimmoObject * , dvecarr3E * >	getDeformedField();
	dvecarr3E		getDisplacements();
	
	int 			setAddNode(MimmoObject*);
	ivector1D		setAddNode(std::vector<MimmoObject*>);

	void			setNode(std::vector<MimmoObject*>)
	void			setFilter(dvector1D );
	
	ivector1D		checkDuplicatedNodes();
	bool 			removeDuplicatedNodes(ivector1D * list=NULL);
	bool			removeNode(int id);
	bool			removeNode(ivector1D & list);
	bool			removeAllNodes();
	
    void            setSupportRadius(double suppR_);
    void            setSupportRadiusValue(double suppR_);
	void 			setTol(double tol);
	void 			setDisplacements(dvecarr3E displ);

	void 			clear();
	void 			clearFilter();
	
	//execute deformation methods
	void 			execute();

	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	void			setWeight(dvector2D value);
	double			calcDist(int i,int j);
	double			calcDist(darray3E & point, int j)

};

}

#endif /* __MRBFSTYLEOBJ_HPP__ */
