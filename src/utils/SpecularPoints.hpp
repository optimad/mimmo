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
#ifndef __SPECULARPOINTS_HPP__
#define __SPECULARPOINTS_HPP__

#include "ProjectCloud.hpp"

namespace mimmo{

/*!
 *    \class SpecularPoints
 *    \ingroup utils
 *    \brief SpecularPoints is a class that mirrors a point cloud w.r.t. a reference plane, on a target geometry if any
 *
 *  SpecularPoints is derived from ProjectCloud class. Given a certain number of points and a reference plane,
 *  the class mirrors such points with respect to this plane. If a geometry is linked, project mirrored points on 
 *  this target geometry. \n
 *  Any data attached, as scalar/vector float data format, are mirrored as well.
 *
 * \n
 * Ports available in GenericInput Class :
 *
 *    =========================================================

   |Port Input | | | |
   |-|-|-|-|
   |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 10    | M_DISPLS      | setVectorData     | (VECARR3, FLOAT)          |
   | 19    | M_SCALARFIELD | setScalarData     | (VECTOR, FLOAT)           |
   | 29    | M_PLANE       | setPlane          | (ARRAY4, FLOAT)           |
   | 20    | M_POINT       | setOrigin         | (ARRAY3, FLOAT)           |
   | 21    | M_AXIS        | setNormal         | (ARRAY3, FLOAT)           |
   | 32    | M_VALUEB      | setInsideOut      | (SCALAR, BOOL)            |
   | 140   | M_VALUEB      | setForce          | (SCALAR, BOOL)            |

   |Port Output | | | |
   |-|-|-|-|
   |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
   | 10    | M_DISPLS      | getCloudVectorData| (VECARR3, FLOAT)       |
   | 19    | M_SCALARFIELD | getCloudScalarData| (VECTOR, FLOAT)        |


  Inherited from ProjectCloud

   |Port Input | | | |
   |-|-|-|-|
   |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 0     | M_COORDS | setCoords         | (VECARR3, FLOAT)    |
   | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)    |

   |Port Output | | | |
   |-|-|-|-|
   |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
   | 0     | M_COORDS | getCloudResult    | (VECARR3, FLOAT)    |


 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SpecularPoints</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 * 
 * Proper of the class:
 * - <B>Force</B>: boolean 0/1. If 1, force mirroring of points that lies on the plane;
 * - <B>InsideOut</B>: boolean 0/1 to get direction of clipping according to given plane;
 * - <B>Plane</B>: section defining the plane's normal and a point belonging to it : \n
 *              <tt> \<Plane\> \n
 *                  \<Point\> 0.0 0.0 0.0 \</Point\> \n
 *                  \<Normal\> 0.0 1.0 0.0 \</Normal\> \n
 *              \</Plane\> </tt> \n

 *
 * Points list and data have to be mandatorily passed through port.
 *
 */
class SpecularPoints: public ProjectCloud{
private:
    bool        m_insideout;        /**< plane direction for mirroring */
    bool        m_force;            /**< if true force the mirroring of points that belong to the symmetry plane. */
    darray4E    m_plane;            /**< reference plane */
    dvector1D   m_scalar;           /**< float scalar data attached to original points*/
    dvecarr3E   m_vector;           /**< float vector data attached to original points*/
    dvector1D   m_scalarMirrored;   /**< resulting float scalar data after mirroring*/
    dvecarr3E   m_vectorMirrored;   /**< resulting float scalar data after mirroring*/
    darray3E    m_origin;           /**<Origin of plane. */
    darray3E    m_normal;           /**<Normal of plane. */
    bool        m_implicit;         /**<True if an implicit definition of plane is set. */

public:
    SpecularPoints();
    SpecularPoints(const bitpit::Config::Section & rootXML);
    ~SpecularPoints();

    SpecularPoints(const SpecularPoints & other);
    SpecularPoints & operator=(const SpecularPoints & other);

    void    buildPorts();

    dvector1D getOriginalScalarData();
    dvecarr3E getOriginalVectorData();

    dvector1D getCloudScalarData();
    dvecarr3E getCloudVectorData();

    darray4E  getPlane();
    bool      isInsideOut();
    bool      isForce();

    void    setVectorData(dvecarr3E data);
    void    setScalarData(dvector1D data);
    void    setPlane(darray4E plane);
    void    setPlane(darray3E origin, darray3E normal);
    void    setOrigin(darray3E origin);
    void    setNormal(darray3E normal);
    void    setInsideOut(bool flag);
    void    setForce(bool flag);

    void     execute();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name= "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    virtual void plotOptionalResults();
};

REGISTER(BaseManipulation, SpecularPoints, "mimmo.SpecularPoints")
};

#endif /* __SPECULARPOINTS_HPP__ */
