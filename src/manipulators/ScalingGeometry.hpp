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
#ifndef __SCALINGGEOMETRY_HPP__
#define __SCALINGGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class ScalingGeometry
 *    \ingroup manipulators
 *    \brief ScalingGeometry is the class that applies a scaling to a given geometry
 *    patch in respect to the mean point of the vertices.
 *
 *    The used parameters are the scaling factor values for each direction of the cartesian
 *    reference system.
 *
 * \n
 * Ports available in ScalingGeometry Class :
 *
 *    =========================================================

     |Port Input | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 23    | M_SPAN   | setScaling        | (ARRAY3, FLOAT)       |
     | 12    | M_FILTER | setFilter         | (VECTOR, FLOAT)       |
     | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 
     |Port Output | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | 11    | M_GDISPLS | getDisplacements  | (VECARR3, FLOAT)      |
     | 80    | M_PAIRVECFIELD | getDeformedField  | (PAIR, MIMMO_VECARR3FLOAT_)  |
 
 *    =========================================================
 *
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ScalingGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Scaling</B>: scaling factor values for each cartesian axis.
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class ScalingGeometry: public BaseManipulation{
private:
    darray3E    m_scaling;        /**<Values of the three fundamental scaling factors (1 = original geometry).*/
    dvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
    ScalingGeometry(darray3E scaling = { {1.0, 1.0, 1.0} });
    ScalingGeometry(const bitpit::Config::Section & rootXML);
    ~ScalingGeometry();

    ScalingGeometry(const ScalingGeometry & other);
    ScalingGeometry & operator=(const ScalingGeometry & other);

    void        buildPorts();

    void        setScaling(darray3E scaling);
    void        setFilter(dvector1D filter);

    dvecarr3E   getDisplacements();
    std::pair<MimmoObject * , dvecarr3E * >    getDeformedField();

    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
};

REGISTER(BaseManipulation, ScalingGeometry, "mimmo.ScalingGeometry")

};

#endif /* __SCALINGGEOMETRY_HPP__ */
