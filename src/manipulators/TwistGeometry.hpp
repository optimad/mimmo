/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
#ifndef __TWISTGEOMETRY_HPP__
#define __TWISTGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\class TwistGeometry
 *	\brief TwistGeometry is the class that applies a twist to a given geometry patch.
 *
 *	The used parameters are a twist max reference value in radians, the direction, the origin
 *	of the twist axis and the distance from the origin where the max value is fixed and maintained
 *	for higher distances from the origin.
 *	The twist is applied following the right-hand rule by using as direction the distance
 *	vector from the origin.
 *
 *	=========================================================
 * ~~~
 *	|--------------------------------------------------------------|
 *	|                 Port Input                                   |
 *	|-------|----------|-------------------|-----------------------|
 *	|PortID | PortType | variable/function | DataType              |
 *	|-------|----------|-------------------|-----------------------|
 *  | 20    | M_POINT  | setOrigin         | (ARRAY3, FLOAT)       |
 *	| 21    | M_AXIS   | setDirection      | (ARRAY3, FLOAT)       |
 *	| 30    | M_VALUED | setTwist          | (SCALAR, FLOAT)       |
 *  | 130   | M_VALUED2| setMaxDistance    | (SCALAR, FLOAT)       |
 *	| 12    | M_FILTER | setFilter         | (VECTOR, FLOAT)       |
 *	| 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |                       |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType              |
 *	|-------|-----------|-------------------|-----------------------|
 *	| 11    | M_GDISPLS | getDisplacements  | (VECARR3, FLOAT)      |
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class TwistGeometry: public BaseManipulation{
private:
	//members
	darray3E	m_origin;		/**<Origin of the twist axis.*/
	darray3E	m_direction;	/**<Components of the twist axis.*/
    double      m_alpha;        /**<Angle of twist in radiant. */
    double      m_distance;     /**<Maximum distance where the input angle of twist is reached. */
	dvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
	TwistGeometry(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
	TwistGeometry(const bitpit::Config::Section & rootXML);
	~TwistGeometry();

	TwistGeometry(const TwistGeometry & other);
	TwistGeometry & operator=(const TwistGeometry & other);

	void        buildPorts();

	void        setAxis(darray3E origin, darray3E direction);
	void        setOrigin(darray3E origin);
	void        setDirection(darray3E direction);
    void        setTwist(double alpha);
    void        setMaxDistance(double distance);
    void        setFilter(dvector1D filter);

    dvecarr3E   getDisplacements();

	void 	    execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
};

REGISTER(BaseManipulation, TwistGeometry, "mimmo.TwistGeometry")

};

#endif /* __TWISTGEOMETRY_HPP__ */
