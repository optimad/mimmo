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
#ifndef __OBBox_HPP__
#define __OBBox_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class OBBox
 *    \ingroup utils
 *    \brief Oriented Bounding Box calculator.
 *
 *    Builds the oriented bounding box of a 3D object (Point Clouds or superficial tessellations), passed as MimmoObject.
 *
 * \n
 * Ports available in OBBox Class :
 *
 *    =========================================================
 
 
     |Port Input | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> | 
     | 99    | M_GEOM      | m_geometry                            |(SCALAR, MIMMO_)       |
     | 100   | M_VECGEOM   | setGeometries                         |(VECTOR, MIMMO_)       |

     |Port Output | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | 20    | M_POINT     | getOrigin         | (ARRAY3, FLOAT)     |
     | 22    | M_AXES      | getAxes           | (ARR3ARR3, FLOAT)   |
     | 23    | M_SPAN      | getSpan           | (ARRAY3, FLOAT)     |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.OBBox</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>ForceAABB</B>: if true calculate the simple AABB of the union of target geometries linked
 * Geometries have to be mandatorily added/passed through ports.
 *
 */
class OBBox: public BaseManipulation {

protected:
    darray3E    m_origin;       /**< Origin of the OBB.*/
    darray3E    m_span;         /**< Span of the OBB. */
    dmatrix33E    m_axes;       /**< reference system of the bbox, ordered aaccording maximum shape variance */

    std::unordered_map<MimmoObject*, int> m_listgeo; /**< list of geometries linked in input, according to type */
    bool m_forceAABB; /**< force class to evaluate a simple AABB, not oriented */
public:
    OBBox();
    OBBox(const bitpit::Config::Section & rootXML);
    virtual ~OBBox();

    //copy operators/constructors
    OBBox(const OBBox & other);

    void buildPorts();

    //clean structure;
    void         clearOBBox();

    //internal methods
    std::vector<MimmoObject*>        getGeometries();
    darray3E                         getOrigin();
    darray3E                         getSpan();
    dmatrix33E                       getAxes();
    bool                             isForcedAABB();
    
    void        setGeometry(MimmoObject* geo);
    void        setGeometries(std::vector<MimmoObject*> listgeo);
    void        setForceAABB(bool flag);
    //plotting wrappers
    void        plot(std::string directory, std::string filename, int counter, bool binary);

    //building method
    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");


protected:
    virtual void plotOptionalResults();

private:
    dmatrix33E      evaluateCovarianceMatrix(dmatrix33E &, darray3E &);
    void            assemblyCovContributes(std::vector<MimmoObject *> list, bool flag, dmatrix33E &);
    darray3E        evaluateMassCenter(std::vector<MimmoObject *> list, bool flag );
    dmatrix33E      eigenVectors( dmatrix33E &, darray3E & eigenValues);
    void            adjustBasis( dmatrix33E &, darray3E & eigenValues);
};

REGISTER(BaseManipulation, OBBox, "mimmo.OBBox")

};

#endif /* __LATTICE_HPP__ */
