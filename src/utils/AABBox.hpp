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
#ifndef __AABBox_HPP__
#define __AABBox_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class AABBox
 *    \ingroup utils
 *    \brief Axis Aligned Bounding Box calculator.
 *
 *    Builds the axis aligned bounding box of a 3D object
      (Point Clouds or superficial tessellations), passed as MimmoObject.

      A Custom Reference System of axes can be provided by the User, instead of the
      Elemental s.d.r.: in this case the AABB will be specialized onto the new sdr of axes.
 *
 * \n
 * Ports available in AABBox Class :
 *
 *    =========================================================


     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM      | setGeometry                           |(MC_SCALAR, MD_MIMMO_)       |
     | M_AXES      | setAxes                               | (MC_ARR3ARR3, MD_FLOAT)   |
     | M_VECGEOM   | setGeometries                         |(MC_VECTOR, MD_MIMMO_)       |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | M_POINT     | getOrigin         | (MC_ARRAY3, MD_FLOAT)     |
     | M_AXES      | getAxes           | (MC_ARR3ARR3, MD_FLOAT)   |
     | M_SPAN      | getSpan           | (MC_ARRAY3, MD_FLOAT)     |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.AABBox</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>WriteInfo</B>: boolean(0/1) if true write info of AABB on file, in plotOptionalResults directory, false do nothing.
 * - <B>Axes</B>: reference system to bind box calculation: \n\n
 *                <tt><B>\<RefSystem\></B> \n
 *                     &nbsp;&nbsp;&nbsp;<B>\<axis0\></B> 1.0 0.0 0.0 <B>\</axis0\></B> \n
 *                     &nbsp;&nbsp;&nbsp;<B>\<axis1\></B> 0.0 1.0 0.0 <B>\</axis1\></B> \n
 *                     &nbsp;&nbsp;&nbsp;<B>\<axis2\></B> 0.0 0.0 1.0 <B>\</axis2\></B> \n
 *                  <B>\</Axes\></B> </tt> \n\n

 * Geometries have to be mandatorily added/passed through ports.
 *
 */
class AABBox: public BaseManipulation {

protected:
    darray3E    m_origin;       /**< Origin of the AABB.*/
    darray3E    m_span;         /**< Span of the AABB. */
    dmatrix33E    m_axes;       /**< reference system of the bbox, ordered aaccording maximum shape variance */

    std::unordered_map<MimmoSharedPointer<MimmoObject>, int> m_listgeo; /**< list of geometries linked in input, according to type */
    bool m_writeInfo; /**< write AABB info on file */

public:
    AABBox();
    AABBox(const bitpit::Config::Section & rootXML);
    virtual ~AABBox();

    //copy operators/constructors
    AABBox(const AABBox & other);

    void buildPorts();

    //clean structure;
    void         clearAABBox();

    //internal methods
    std::vector<MimmoSharedPointer<MimmoObject> >   getGeometries();
    darray3E                         getOrigin();
    darray3E                         getSpan();
    dmatrix33E                       getAxes();

    void        setGeometry(MimmoSharedPointer<MimmoObject> geo);
    void        setGeometries(std::vector<MimmoSharedPointer<MimmoObject> > listgeo);
    void        setAxes(dmatrix33E axes);
    void        setWriteInfo(bool flag);

    //plotting wrappers
    void        plot(std::string directory, std::string filename, int counter, bool binary);

    //building method
    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");


protected:
    virtual void plotOptionalResults();
    void swap(AABBox & x) noexcept;
    dmatrix33E transpose(const dmatrix33E & mat);
    dmatrix33E inverse (const dmatrix33E & mat);

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__AABBox_HPP__)
REGISTER_PORT(M_VECGEOM, MC_VECTOR, MD_MIMMO_,__AABBox_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__AABBox_HPP__)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT,__AABBox_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT,__AABBox_HPP__)


REGISTER(BaseManipulation, AABBox, "mimmo.AABBox")

};

#endif /* __LATTICE_HPP__ */
