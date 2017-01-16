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
#ifndef __CLIPGEOMETRY_HPP__
#define __CLIPGEOMETRY_HPP__

#include "BaseManipulation.hpp"
#include "MimmoObject.hpp"

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief ClipGeometry is a class that clip a 3D geometry according to a plane intersecting it.
 *
 *	ClipGeometry is derived from BaseManipulation class.
 *	It needs a target MimmoObject geometry, alongside a clipping plane definition [a,b,c,d], in its
 *  implicit form a*x + b*y + c*z + d =0; 
 *	Returns geometry clipped in an independent MimmoObject. Controls clipping direction with an "insideout" boolean.
 *  Class plot as optional result the clipped portion of geoemetry, and absorb/flush its parameters from xml 
 *  dictionaries
 *
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------|
 *	|                 Port Input                       		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType 		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *	| 29    | M_PLANE  | setClipPlane      | (ARRAY4, FLOAT)		 |
 *	| 32    | M_VALUEB | setInsideOut      | (SCALAR, BOOL)			 |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)		 |
 *	|-------|----------|-------------------|-------------------------|
 *
 *
 *	|----------------------------------------------------------------|
 *	|            Port Output                          		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType 		 		 |
 *	|-------|----------|-------------------|-------------------------|
 *	| 99    | M_GEOM   | getClippedPatch   | (SCALAR, MIMMO_)		 |
 *	|-------|----------|-------------------|-------------------------|
 * ~~~
 *	=========================================================
 *
 */
class ClipGeometry: public BaseManipulation{
private:
	darray4E						m_plane;		/**<Coefficients of implicit plane a*x +b*y+c*z + d =0.*/
	bool						m_insideout;	/**<set direction of clipping, false along current plane normal, true the opposite*/
	std::unique_ptr<MimmoObject>	m_patch;		/**<Resulting Clipped Patch.*/
	
public:
	ClipGeometry();
	~ClipGeometry();

	ClipGeometry(const ClipGeometry & other);
	ClipGeometry & operator=(const ClipGeometry & other);

	void	buildPorts();

	bool 			isInsideOut();
	MimmoObject * 	getClippedPatch();
	darray4E 		getClipPlane();	

	void	setClipPlane(darray4E plane);
	void	setClipPlane(darray3E origin, darray3E normal);
	void	setInsideOut(bool flag);

	void 	execute();
	
	virtual void absorbSectionXML(bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	
protected:
	virtual void plotOptionalResults();
	
private:	
	livector1D	clipPlane();

};

}

#endif /* __CLIPGEOMETRY_HPP__ */
