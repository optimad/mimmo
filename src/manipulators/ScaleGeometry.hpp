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
#ifndef __SCALEGEOMETRY_HPP__
#define __SCALEGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class ScaleGeometry
 *    \ingroup manipulators
 *    \brief ScaleGeometry is the class that applies a scaling to a given geometry
 *    patch in respect to the mean point of the vertices.
 *
 *    The used parameters are the scaling factor values for each direction of the cartesian
 *    reference system.
 *
 * \n
 * Ports available in ScaleGeometry Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT  | setOrigin         | (MC_ARRAY3, MD_FLOAT)       |
     | M_SPAN   | setScaling        | (MC_ARRAY3, MD_FLOAT)       |
     | M_FILTER | setFilter         | (MC_VECTOR, MD_FLOAT)       |
     | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)      |
 
     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GDISPLS | getDisplacements  | (MC_VECARR3, MD_FLOAT)      |
     | M_GEOM   | getGeometry       | (MC_SCALAR,MD_MIMMO_) |
 
 *    =========================================================
 *
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ScaleGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 *
 * Proper of the class:
 * - <B>Origin</B>: scaling center point;
 * - <B>MeanPoint</B>: use mean point as center of scaling;
 * - <B>Scaling</B>: scaling factor values for each cartesian axis.
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class ScaleGeometry: public BaseManipulation{
private:
    darray3E    m_origin;        /**<Center point of scaling.*/
    darray3E    m_scaling;        /**<Values of the three fundamental scaling factors (1 = original geometry).*/
    dmpvector1D   m_filter;       /**<Filter field for displacements modulation. */
    dmpvecarr3E   m_displ;        /**<Resulting displacements of geometry vertex.*/
    bool        m_meanP;       /**<Use mean point as center of scaling.*/

public:
    ScaleGeometry(darray3E scaling = { {1.0, 1.0, 1.0} });
    ScaleGeometry(const bitpit::Config::Section & rootXML);
    ~ScaleGeometry();

    ScaleGeometry(const ScaleGeometry & other);
    ScaleGeometry& operator=(ScaleGeometry other);
    
    void        buildPorts();

    void        setScaling(darray3E scaling);
    void        setFilter(dmpvector1D filter);
    void        setOrigin(darray3E origin);
    void        setMeanPoint(bool meanP);

    dmpvecarr3E   getDisplacements();

    void         execute();
    void         apply();
    void         checkFilter();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(ScaleGeometry & x) noexcept;
};

REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__SCALEGEOMETRY_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT,__SCALEGEOMETRY_HPP__)
REGISTER_PORT(M_FILTER, MC_VECTOR, MD_FLOAT,__SCALEGEOMETRY_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__SCALEGEOMETRY_HPP__)
REGISTER_PORT(M_GDISPLS, MC_VECARR3, MD_FLOAT,__SCALEGEOMETRY_HPP__)
REGISTER_PORT(M_PAIRVECFIELD, MC_PAIR, MD_MIMMO_VECARR3FLOAT_,__SCALEGEOMETRY_HPP__)


REGISTER(BaseManipulation, ScaleGeometry, "mimmo.ScaleGeometry")

};

#endif /* __SCALEGEOMETRY_HPP__ */
