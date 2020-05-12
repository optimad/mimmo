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
 * It returns geometry clipped in an independent MimmoObject.
 * The class controls clipping direction (based on normal direction) with an "insideout" boolean.
 * ClipGeometry plots as optional result the clipped portion of input geometry.
 *
 * Ports available in ClipGeometry Class :
 *
 *    =========================================================
 *
     |                 Port Input    ||                              |
     |----------|-------------------|-------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_PLANE  | setClipPlane      | (MC_ARRAY4, MD_FLOAT)         |
     | M_POINT  | setOrigin         | (MC_ARRAY3, MD_FLOAT)         |
     | M_AXIS   | setNormal         | (MC_ARRAY3, MD_FLOAT)         |
     | M_VALUEB | setInsideOut      | (MC_SCALAR, MD_BOOL)          |
     | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)        |


     |            Port Output           ||                           |
     |----------|-------------------|-------------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM   | getClippedPatch   | (MC_SCALAR, MD_MIMMO_)        |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ClipGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>InsideOut</B>: boolean 0/1 to get direction of clipping according to given plane normal;
 * - <B>ClipPlane</B>: section defining the plane's normal and a point belonging to it: \n\n
 *         <tt> <B>\<ClipPlane\></B> \n
 *              &nbsp;&nbsp;&nbsp;<B>\<Point\></B> 0.0 0.0 0.0 <B>\</Point\></B> \n
 *              &nbsp;&nbsp;&nbsp;<B>\<Normal\></B> 0.0 1.0 0.0 <B>\</Normal\></B> \n
 *              <B>\</ClipPlane\></B> </tt>;
 *
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class ClipGeometry: public BaseManipulation{
private:
    darray4E                        m_plane;        /**<Coefficients of implicit plane a*x +b*y+c*z + d =0.*/
    bool                            m_insideout;    /**<set direction of clipping, false along current plane normal, true the opposite*/
    MimmoSharedPointer<MimmoObject> m_patch;        /**<Resulting Clipped Patch.*/
    darray3E                        m_origin;       /**<Origin of plane. */
    darray3E                        m_normal;       /**<Normal of plane. */
    bool                            m_implicit;     /**<True if an implicit definition of plane is set. */

public:
    ClipGeometry();
    ClipGeometry(const bitpit::Config::Section & rootXML);
    ~ClipGeometry();

    ClipGeometry(const ClipGeometry & other);

    void    buildPorts();

    bool             isInsideOut();
    MimmoSharedPointer<MimmoObject>    getClippedPatch();
    darray4E         getClipPlane();

    void    setClipPlane(darray4E plane);
    void    setClipPlane(darray3E origin, darray3E normal);
    void    setOrigin(darray3E origin);
    void    setNormal(darray3E normal);
    void    setInsideOut(bool flag);
    using BaseManipulation::setGeometry;

    void     execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    virtual void plotOptionalResults();

private:
    livector1D    clipPlane();

};


REGISTER_PORT(M_PLANE, MC_ARRAY4, MD_FLOAT,__CLIPGEOMETRY_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__CLIPGEOMETRY_HPP__)
REGISTER_PORT(M_AXIS, MC_ARRAY3, MD_FLOAT,__CLIPGEOMETRY_HPP__)
REGISTER_PORT(M_VALUEB, MC_SCALAR, MD_BOOL,__CLIPGEOMETRY_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__CLIPGEOMETRY_HPP__)

REGISTER(BaseManipulation, ClipGeometry, "mimmo.ClipGeometry")
}

#endif /* __CLIPGEOMETRY_HPP__ */
