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
#ifndef __STITCHGEOMETRY_HPP__
#define __STITCHGEOMETRY_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\date			30/nov/2016
 *	\authors		Rocco Arpa
 *	\authors		
 *
 *	\brief StitchGeometry is an executable block class capable of  
 *         stitch multiple MimmoObject geometries of the same topology 
 *
 *	StitchGeometry is the object to append two or multiple MimmoObject of the same topology
 *  in a unique MimmoObject container. 
 * 
 *	=========================================================
 * ~~~
 *	|-------------------------------------------------------------------------------|
 *	|                 Port Input                           					        |
 *	|-------|----------|------------------------------------|-----------------------|
 *	|PortID | PortType | variable/function 					| DataType		        |
 *	|-------|----------|------------------------------------|-----------------------|
 *  | 99    | M_GEOM   | setAddGeometry						| (SCALAR, MIMMO_)      |  
 *	| 100   | M_VECGEOM| setGeometry      					| (VECTOR, MIMMO_)	    |
 *	|-------|----------|------------------------------------|-----------------------|
 *
 *
 *	|--------------------------------------------------------|-----------------------|
 *	|            Port Output                				 |                       |
 *	|-------|-----------|------------------------------------|-----------------------|
 *	|PortID | PortType  | variable/function 				 | DataType		      	 | 
 *	|-------|-----------|------------------------------------|-----------------------|
 *	| 99    | M_GEOM    | getGeometry      					 | (SCALAR, MIMMO_)	     |
 *  | 100   | M_VECGEOM | getOriginalGeometries				 | (VECTOR, MIMMO_)	     |
 *  | 104   | M_MAPDCELL| getCellDivisionMap   				 | (UN_MAP, LONGPAIRINTLONG) |
 * 	| 105	| M_MAPDVERT| getVertDivisionMap				 | (UN_MAP, LONGPAIRINTLONG) |	
 *	|-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */

class StitchGeometry: public BaseManipulation{
	
private:
	int 									m_topo;		/**<Mark topology of your stitcher 1-surface, 2-volume, 3-pointcloud*/
	std::unordered_map<MimmoObject*,int>	m_extgeo;	/**< pointers to external geometries*/
	
	std::unique_ptr<MimmoObject> m_patch;	/**< resulting patch geometry */	
	
	std::unordered_map<long, std::pair<int, long> > m_mapCellDivision; /**< division map of actual ID-cell, part Id, original ID-cell*/
	std::unordered_map<long, std::pair<int, long> > m_mapVertDivision; /**< division map of actual ID-vertex, part Id, original ID-vertex*/
	
	bool		m_buildBvTree;				/**<If true build BvTree of stitched geometry. */
	bool		m_buildKdTree;				/**<If true build KdTree of stitched geometry. */

	int m_geocount;							/**<Internal geometry counter */
	
public:
	StitchGeometry(int topo);
	virtual ~StitchGeometry();

	StitchGeometry(const StitchGeometry & other);
	StitchGeometry & operator=(const StitchGeometry & other);

	void buildPorts();
	
	int 							getTopology();
	std::vector< MimmoObject *> 	getOriginalGeometries();
	MimmoObject * 					getGeometry();
	
	std::unordered_map<long, std::pair<int, long> >	getCellDivisionMap();
	std::unordered_map<long, std::pair<int, long> >	getVertDivisionMap();
	
	void		setAddGeometry(MimmoObject * geo);
	void		setGeometry( std::vector<MimmoObject *> list);
	
	void		setBuildBvTree(bool build);
	void		setBuildKdTree(bool build);
		
	bool 		isEmpty();
	
	void 		clear();
	void 		execute();
	
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
	
	void 	plotOptionalResults();	
};


}

#endif /* __STITCHGEOMETRY_HPP__ */
