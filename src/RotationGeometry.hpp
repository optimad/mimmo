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
#ifndef __ROTATIONGEOMETRY_HPP__
#define __ROTATIONGEOMETRY_HPP__

#include "BaseManipulation.hpp"
#include "FFDLattice.hpp"

namespace mimmo{

/*!
 *	\date			10/oct/2016
 *	\authors		Edoardo Lombardi
 *
 *	\brief RotationGeometry is the class that applies a rotation to a given geometry patch.
 *
 *	The used parameters are the rotation value and the direction and the origin
 *	of the rotation axis.
 *
 *	=========================================================
 * ~~~
 *	|--------------------------------------------------------------|
 *	|                 Port Input                                   |
 *	|-------|----------|-------------------|-----------------------|
 *	|PortID | PortType | variable/function | DataType		       |
 *	|-------|----------|-------------------|-----------------------|
 *	| 20    | M_POINT  | setOrigin         | (ARRAY3, FLOAT)	   |
 *	| 21    | M_AXIS   | setDirection      | (ARRAY3, FLOAT)	   |
 *	| 30    | M_VALUED | setRotation       | (SCALAR, FLOAT)	   |
 *  | 12    | M_FILTER | setFilter         | (VECTOR, FLOAT)       |
 *  | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |                       |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType		 		|
 *	|-------|-----------|-------------------|-----------------------|
 *  | 11    | M_GDISPLS | getDisplacements  | (VECARR3, FLOAT)      |
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class RotationGeometry: public BaseManipulation{
private:
	//members
	darray3E	m_origin;		/**<Origin of the rotation axis.*/
	darray3E	m_direction;	/**<Components of the rotation axis.*/
	double		m_alpha;        /**<Angle of rotation in radiant. */
	dvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
	RotationGeometry(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
	RotationGeometry(const bitpit::Config::Section & rootXML);
	~RotationGeometry();

	RotationGeometry(const RotationGeometry & other);
	RotationGeometry & operator=(const RotationGeometry & other);

	void        buildPorts();

	void        setAxis(darray3E origin, darray3E direction);
	void        setOrigin(darray3E origin);
	void        setDirection(darray3E direction);
    void        setRotation(double alpha);
    void        setFilter(dvector1D filter);

    dvecarr3E   getDisplacements();

	void 	    execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
};

}

#endif /* __ROTATIONGEOMETRY_HPP__ */
