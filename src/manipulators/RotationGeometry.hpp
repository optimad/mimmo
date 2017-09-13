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
#ifndef __ROTATIONGEOMETRY_HPP__
#define __ROTATIONGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class RotationGeometry
 *    \ingroup manipulators
 *    \brief RotationGeometry is the class that applies a rotation to a given geometry patch.
 *
 *    The used parameters are the rotation value and the direction and the origin
 *    of the rotation axis.
 *
 * \n
 * Ports available in RotationGeometry Class :
 *
 *    =========================================================

     |Port Input | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 20    | M_POINT  | setOrigin         | (ARRAY3, FLOAT)       |
     | 21    | M_AXIS   | setDirection      | (ARRAY3, FLOAT)       |
     | 30    | M_VALUED | setRotation       | (SCALAR, FLOAT)       |
     | 12    | M_FILTER | setFilter         | (MPVECTOR, FLOAT)       |
     | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 
     |Port Output | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | 11    | M_GDISPLS | getDisplacements  | (MPVECARR3, FLOAT)      |
     | 99    | M_GEOM   | getGeometry       | (SCALAR,MIMMO_) |
 
 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.RotationGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 *
 * Proper of the class:
 * - <B>Origin</B>: rotation axis origin;
 * - <B>Direction</B>: axis direction coordinates;
 * - <B>Rotation</B>: rotation angle in radians. Positive on counterclockwise rotations around reference axis.
 *
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class RotationGeometry: public BaseManipulation{
private:
    //members
    darray3E    m_origin;        /**<Origin of the rotation axis.*/
    darray3E    m_direction;    /**<Components of the rotation axis.*/
    double        m_alpha;        /**<Angle of rotation in radiant. */
    dmpvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dmpvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
    RotationGeometry(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
    RotationGeometry(const bitpit::Config::Section & rootXML);
    ~RotationGeometry();

    RotationGeometry(const RotationGeometry & other);
    RotationGeometry & operator=(RotationGeometry other);
    
    void        buildPorts();

    void        setAxis(darray3E origin, darray3E direction);
    void        setOrigin(darray3E origin);
    void        setDirection(darray3E direction);
    void        setRotation(double alpha);
    void        setFilter(dmpvector1D filter);

    dmpvecarr3E   getDisplacements();

    void         execute();
    void         apply();
    void         checkFilter();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(RotationGeometry & x) noexcept;
};

REGISTER(BaseManipulation, RotationGeometry, "mimmo.RotationGeometry")

};

#endif /* __ROTATIONGEOMETRY_HPP__ */
