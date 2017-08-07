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
#ifndef __TWISTGEOMETRY_HPP__
#define __TWISTGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class TwistGeometry
 *    \ingroup manipulators
 *    \brief TwistGeometry is the class that applies a twist to a given geometry patch.
 *
 *    The used parameters are a twist max reference value in radians, the direction, the origin
 *    of the twist axis and the distance from the origin where the max value is fixed and maintained
 *    for higher distances from the origin.
 *    The twist is applied following the right-hand rule by using as direction the distance
 *    vector from the origin.
 *
 * \n
 * Ports available in TwistGeometry Class :
 *
 *    =========================================================
 
     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT  | setOrigin         | (MC_ARRAY3, MD_FLOAT)       |
     | M_AXIS   | setDirection      | (MC_ARRAY3, MD_FLOAT)       |
     | M_VALUED | setTwist          | (MC_SCALAR, MD_FLOAT)       |
     | M_VALUED2| setMaxDistance    | (MC_SCALAR, MD_FLOAT)       |
     | M_FILTER | setFilter         | (MC_VECTOR, MD_FLOAT)       |
     | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)      |
 
     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GDISPLS       | getDisplacements  | (MC_VECARR3, MD_FLOAT) |
     | M_GEOM          | getGeometry       | (MC_SCALAR,MD_MIMMO_) |
 
 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.TwistGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply result directly in execution;
 *
 * Proper of the class:
 * - <B>Origin</B>: twist axis origin;
 * - <B>Direction</B>: axis direction coordinates;
 * - <B>Twist</B>: twist angle in radians. Positive on counterclockwise twists around reference axis;
 * - <B>Distance</B>: distance from the origin (on twist axis) where the twist angle is reached.
 * - <B>Symmetric</B>: boolean 0/1 activate symmetric twisting 
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class TwistGeometry: public BaseManipulation{
private:
    //members
    darray3E    m_origin;        /**<Origin of the twist axis.*/
    darray3E    m_direction;    /**<Components of the twist axis.*/
    double      m_alpha;        /**<Angle of twist in radiant. */
    double      m_distance;     /**<Maximum distance where the input angle of twist is reached. */
    dmpvector1D m_filter;       /**<Filter field for displacements modulation. */
    dmpvecarr3E m_displ;        /**<Resulting displacements of geometry vertex.*/
    bool        m_sym;          /**<Propagate twist for negative local coordinate.*/

public:
    TwistGeometry(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
    TwistGeometry(const bitpit::Config::Section & rootXML);
    ~TwistGeometry();

    TwistGeometry(const TwistGeometry & other);
    TwistGeometry & operator=(TwistGeometry other);

    void        buildPorts();

    void        setAxis(darray3E origin, darray3E direction);
    void        setOrigin(darray3E origin);
    void        setDirection(darray3E direction);
    void        setTwist(double alpha);
    void        setSym(bool sym);
    void        setMaxDistance(double distance);
    void        setFilter(dmpvector1D filter);

    dmpvecarr3E getDisplacements();

    void         execute();
    void         apply();
    void         checkFilter();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(TwistGeometry & x) noexcept;
};

REGISTER(BaseManipulation, TwistGeometry, "mimmo.TwistGeometry")

};

#endif /* __TWISTGEOMETRY_HPP__ */
