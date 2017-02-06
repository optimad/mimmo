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
#ifndef __PROJPRIMITIVESONSURFACES_HPP__
#define __PROJPRIMITIVESONSURFACES_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief ProjPrimitivesOnSurfaces is an executable block class capable of  
 *         projecting elemental 1D or 2D primitives such as segments, circumferences or circles, triangles etc...
 *         on a 3D surface mesh defined by a MimmoObject. 
 *
 *	ProjPrimitivesOnSurfaces is a pure virtual class 
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|PortID | PortType | variable/function 					| DataType		        |
 *	|-------|----------|------------------------------------|-----------------------|
 *  | 99    | M_GEOM   | setGeometry						| (SCALAR, MIMMO_)      |  
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 99    | M_GEOM    | getProjectedElement 				 | (SCALAR, MIMMO_)	     |
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */

class ProjPrimitivesOnSurfaces: public BaseManipulation{
	
protected:
	int 				m_topo;			/**<Mark topology of your primitive element 1-one dimensional, 2- bi-dimensional*/
	int 				m_nC;           /**<Number of target elements of your 3D curve discrete projection */
	bool				m_buildBvTree;	/**<If true build BvTree of of the projected element mesh. */
	bool				m_buildKdTree;	/**<If true build KdTree of the projected element mesh. */
	std::unique_ptr<MimmoObject> m_patch;			/**< resulting projected elements stored as MimmoObject */	

	
public:
	ProjPrimitivesOnSurfaces();
	virtual ~ProjPrimitivesOnSurfaces();

	ProjPrimitivesOnSurfaces(const ProjPrimitivesOnSurfaces & other);
	ProjPrimitivesOnSurfaces & operator=(const ProjPrimitivesOnSurfaces & other);

	virtual void buildPorts();
	
	int 							getTopology();
	int 							getProjElementTargetNCells();
	MimmoObject * 					getProjectedElement();
	
	void		setGeometry(MimmoObject * geo);
	
	void		setBuildBvTree(bool build);
	void		setBuildKdTree(bool build);
	void 		setProjElementTargetNCells(int nC);
	
	bool 		isEmpty();
	
	void 		clear();
	void 		execute();

	void		plotOptionalResults();
	
protected:
	virtual void projection() = 0;
};

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief ProjSegmentOnSurface is an executable block class capable of  
 *         projecting elemental 1D segments on a 3D surface mesh defined by a MimmoObject. 
 *
 *	ProjSegmentOnSurface project a given segment on a given 3D surface mesh and return a discrete 3D curve mesh 
 *  in MimmoObject container
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|PortID | PortType | variable/function 					| DataType		        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */

class ProjSegmentOnSurface: public ProjPrimitivesOnSurfaces{
	
private:
	darray3E				 m_pointA;		/**<origin of your segment*/
	darray3E				 m_pointB;		/**<end of your segment*/
	
public:
	ProjSegmentOnSurface();
	virtual ~ProjSegmentOnSurface();
	
	ProjSegmentOnSurface(const ProjSegmentOnSurface & other);
	ProjSegmentOnSurface & operator=(const ProjSegmentOnSurface & other);
	
	void 		clear();
	
	void 		setSegment(darray3E pointA, darray3E pointB);
	void		setSegment(darray3E origin, darray3E dir, double length);
	
	void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
protected:
	void projection();
};

}

#endif /* __PROJPRIMITIVESONSURFACES_HPP__ */
