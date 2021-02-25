/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#ifndef __ROTATIONAXES_HPP__
#define __ROTATIONAXES_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\class RotationAxes
 *	\ingroup utils
 *	\brief RotationAxes is the class that applies a rotation to a given reference system.
 *
 *	The used parameters are the rotation value (in randians), the rotation axis
    direction and origin.
 *
 * Ports available in GenericInput Class :
 *
 *	=========================================================
 *

    |Port Input | | |
    |-|-|-|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_POINT           | setOrigin                 | (MC_ARRAY3, MD_FLOAT)       |
 	| M_AXIS            | setDirection              | (MC_ARRAY3, MD_FLOAT)       |
 	| M_VALUED          | setRotation               | (MC_SCALAR, MD_FLOAT)       |
 	| M_POINT2          | setAxesOrigin             | (MC_ARRAY3, MD_FLOAT)       |
 	| M_AXES            | setAxes                   | (MC_ARR3ARR3, MD_FLOAT)     |


    |Port Output | | |
 	|-|-|-|
    | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
    | M_POINT   | getRotatedOrigin  | (MC_ARRAY3, MD_FLOAT)       |
    | M_AXES    | getRotatedAxes    | (MC_ARR3ARR3, MD_FLOAT)     |

 *	=========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.RotationAxes</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>Origin</B>: rotation axis origin coordinates (space separated);
 * - <B>Direction</B>: axis direction coordinates (space separated);
 * - <B>Rotation</B>: rotation angle in radians. Positive on counterclockwise rotations around reference axis;
 * - <B>RefSystem</B>: axes of current shape reference system. written in XML as: \n
 *   <tt><B>\<RefSystem\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<axis0\></B> 1.0 0.0 0.0 <B>\</axis0\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<axis1\></B> 0.0 1.0 0.0 <B>\</axis1\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<axis2\></B> 0.0 0.0 1.0 <B>\</axis2\></B> \n
 *   <B>\</RefSystem\></B> </tt> \n\n

 * - <B>OriginRS</B>: origin coordinates of the target reference system to be rotated (space separated);
 *
 */
class RotationAxes: public BaseManipulation{

protected:
	//members
	darray3E	m_origin;		/**<Origin of the rotation axis.*/
	darray3E	m_direction;	/**<Components of the rotation axis.*/
	dmatrix33E	m_axes;			/**<Axes of box to be deformed.*/
	darray3E	m_axes_origin;	/**<Origin of the axes to be rotated. */
	dmatrix33E	m_rotax;		/**<Axes of box deformed.*/
	darray3E	m_rotax_origin;	/**<Origin of the axes rotated. */
	double		m_alpha;        /**<Rotation angle. */

public:
	RotationAxes(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
	RotationAxes(const bitpit::Config::Section & rootXML);
	~RotationAxes();

	RotationAxes(const RotationAxes & other);
	RotationAxes & operator=(RotationAxes other);

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

protected:
    void swap(RotationAxes & x) noexcept;

};

REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__ROTATIONAXES_HPP__)
REGISTER_PORT(M_AXIS, MC_ARRAY3, MD_FLOAT,__ROTATIONAXES_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT,__ROTATIONAXES_HPP__)
REGISTER_PORT(M_POINT2, MC_ARRAY3, MD_FLOAT,__ROTATIONAXES_HPP__)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT,__ROTATIONAXES_HPP__)

REGISTER(BaseManipulation, RotationAxes, "mimmo.RotationAxes")

};

#endif /* __ROTATIONAXES_HPP__ */
