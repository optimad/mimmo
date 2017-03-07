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
#ifndef __ROTATIONBOX_HPP__
#define __ROTATIONBOX_HPP__

#include "BaseManipulation.hpp"
#include "FFDLattice.hpp"

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief RotationBox is the class that applies a rotation to a given reference system.
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
 *	| 20    | M_POINT  | m_origin          | (ARRAY3, FLOAT)	   |
 *	| 21    | M_AXIS   | m_direction       | (ARRAY3, FLOAT)	   |
 *	| 30    | M_VALUED | m_alpha           | (SCALAR, FLOAT)	   |
 *	| 120   | M_POINT2 | m_axes_origin     | (ARRAY3, FLOAT)	   |
 *	| 22    | M_AXES   | m_axes            | (ARR3ARR3, FLOAT)	   |
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|---------------------------------------|-----------------------|
 *	|            Port Output                |                       |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType		 		|
 *	|-------|-----------|-------------------|-----------------------|
 *	| 20    | M_POINT   | getRotatedOrigin  | (ARRAY3, FLOAT)	  	|
 *	| 22    | M_AXES    | getRotatedAxes    | (ARR3ARR3, FLOAT)	  	|
 *	|-------|-----------|-------------------|-----------------------|
 * ~~~
 *	=========================================================
 *
 */
class RotationBox: public BaseManipulation{
private:
	//members
	darray3E	m_origin;		/**<Origin of the rotation axis.*/
	darray3E	m_direction;	/**<Components of the rotation axis.*/
	dmatrix33E	m_axes;			/**<Axes of box to be deformed.*/
	darray3E	m_axes_origin;	/**<Origin of the axes to be rotated. */
	dmatrix33E	m_rotax;		/**<Axes of box deformed.*/
	darray3E	m_rotax_origin;	/**<Origin of the axes rotated. */
	double		m_alpha;

public:
	RotationBox(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
	RotationBox(const bitpit::Config::Section & rootXML);
	~RotationBox();

	RotationBox(const RotationBox & other);
	RotationBox & operator=(const RotationBox & other);

	void buildPorts();

	void setAxis(darray3E origin, darray3E direction);
	void setOrigin(darray3E origin);
	void setDirection(darray3E direction);
	void setRotation(double alpha);
	void setAxes(dmatrix33E axes);
	void setAxesOrigin(darray3E axes_origin);

	dmatrix33E getRotatedAxes();
	darray3E getRotatedOrigin();

	void 	execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	

};

REGISTER(BaseManipulation, RotationBox, "MiMMO.RotationBox")

};

#endif /* __ROTATIONBOX_HPP__ */
