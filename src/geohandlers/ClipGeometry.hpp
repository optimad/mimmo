/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __CLIPGEOMETRY_HPP__
#define __CLIPGEOMETRY_HPP__

#include "BaseManipulation.hpp"
#include "MimmoObject.hpp"

namespace mimmo{

/*!
 * \class ClipGeometry
 * \ingroup geohandlers
 * \brief ClipGeometry is a class that clip a 3D geometry according to a plane intersecting it.
 *
 * ClipGeometry is derived from BaseManipulation class.
 * It needs a target MimmoObject geometry, alongside a clipping plane definition through its origin
 * and normal or by a set af value [a,b,c,d], describing its
 * implicit form a*x + b*y + c*z + d = 0;
 * Returns geometry clipped in an independent MimmoObject.
 * Controls clipping direction (based on normal direction) with an "insideout" boolean.
 * Class plot as optional result the clipped portion of geoemetry, and absorb/flush its parameters from xml 
 * dictionaries
 *
 *	=========================================================
 * ~~~
 *	|----------------------------------------------------------------|
 *	|                 Port Input                                     |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType                |
 *	|-------|----------|-------------------|-------------------------|
 *  | 29    | M_PLANE  | setClipPlane      | (ARRAY4, FLOAT)         |
 *  | 20    | M_POINT  | setOrigin         | (ARRAY3, FLOAT)         |
 *  | 21    | M_AXIS   | setNormal         | (ARRAY3, FLOAT)         |
 *	| 32    | M_VALUEB | setInsideOut      | (SCALAR, BOOL)          |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)        |
 *	|-------|----------|-------------------|-------------------------|
 *
 *
 *	|----------------------------------------------------------------|
 *	|            Port Output                                         |
 *	|-------|----------|-------------------|-------------------------|
 *	|PortID | PortType | variable/function | DataType                |
 *	|-------|----------|-------------------|-------------------------|
 *	| 99    | M_GEOM   | getClippedPatch   | (SCALAR, MIMMO_)        |
 *	|-------|----------|-------------------|-------------------------|
 * ~~~
 *	=========================================================
 *
 */
class ClipGeometry: public BaseManipulation{
private:
    darray4E                        m_plane;    /**<Coefficients of implicit plane a*x +b*y+c*z + d =0.*/
	bool						m_insideout;	/**<set direction of clipping, false along current plane normal, true the opposite*/
	std::unique_ptr<MimmoObject>	m_patch;	/**<Resulting Clipped Patch.*/
    darray3E                        m_origin;   /**<Origin of plane. */
    darray3E                        m_normal;   /**<Normal of plane. */
    bool                            m_implicit; /**<True if an implicit definition of plane is set. */
	
public:
	ClipGeometry();
	ClipGeometry(const bitpit::Config::Section & rootXML);
	~ClipGeometry();

	ClipGeometry(const ClipGeometry & other);
	ClipGeometry & operator=(const ClipGeometry & other);

	void	buildPorts();

	bool 			isInsideOut();
	MimmoObject * 	getClippedPatch();
	darray4E 		getClipPlane();	

	void	setClipPlane(darray4E plane);
    void    setClipPlane(darray3E origin, darray3E normal);
    void    setOrigin(darray3E origin);
    void    setNormal(darray3E normal);
	void	setInsideOut(bool flag);

	void 	execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	
protected:
	virtual void plotOptionalResults();
	
private:	
	livector1D	clipPlane();

};

REGISTER(BaseManipulation, ClipGeometry, "mimmo.ClipGeometry")
}

#endif /* __CLIPGEOMETRY_HPP__ */
