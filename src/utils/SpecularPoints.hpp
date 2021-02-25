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
#ifndef __SPECULARPOINTS_HPP__
#define __SPECULARPOINTS_HPP__

#include "ProjPatchOnSurface.hpp"

namespace mimmo{

/*!
 *  \class SpecularPoints
 *  \ingroup utils
 *  \brief SpecularPoints is a class that mirrors a point cloud w.r.t. a
             reference plane, on a target surface geometry if any
 *
 *  SpecularPoints is a custom derivation of ProjPatchOnSurface class, specialized for Point Clouds.
    Given a point cloud and a reference plane, the class mirrors such points with
    respect to this plane. If a surface geometry is linked,
    it projects the final mirrored points on it. \n
 *  Any data attached, as scalar/vector float data format, are mirrored as well.
 *
 * Ports available in SpecularPoints Class :
 *
 *    =========================================================

 |Port Input | | |
 |-|-|-|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM2       | setPointCloud     | (MC_SCALAR, MD_MIMMO_)     |
 | M_GEOM        | setGeometry       | (MC_SCALAR, MD_MIMMO_)     |
 | M_VECTORFIELD | setVectorData     | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
 | M_SCALARFIELD | setScalarData     | (MC_SCALAR, MD_MPVECFLOAT_)     |
 | M_PLANE       | setPlane          | (MC_ARRAY4, MD_FLOAT)           |
 | M_POINT       | setOrigin         | (MC_ARRAY3, MD_FLOAT)           |
 | M_AXIS        | setNormal         | (MC_ARRAY3, MD_FLOAT)           |

 |Port Output | | |
 |-|-|-|
 | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>           |
 | M_VECTORFIELD | getMirroredVectorData  | (MC_SCALAR, MD_MPVECARR3FLOAT_) |
 | M_SCALARFIELD | getMirroredScalarData | (MC_SCALAR, MD_MPVECFLOAT_)     |
 | M_GEOM        | getMirroredPointCloud  | (MC_SCALAR, MD_MIMMO_)      |


 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SpecularPoints</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from ProjPatchOnSurface:
 * - <B>KdTree</B> : evaluate kdTree true 1/false 0 for final Point Cloud;
 *
 * Proper of the class:
 * - <B>Force</B>: boolean 0/1. If true, force mirroring of points that lies on the plane;
 * - <B>InsideOut</B>: boolean 0/1 to align/reverse mirroring with/against the plane normal;
 * - <B>Plane</B>: section defining the plane's normal and a point belonging to it : \n\n
 *   <tt> <B>\<Plane\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<Point\></B> 0.0 0.0 0.0 <B>\</Point\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<Normal\></B> 0.0 1.0 0.0 <B>\</Normal\></B> \n
 *   <B>\</Plane\></B> </tt> \n\n
 *
 * Points list and data have to be mandatorily passed through port.
 *
 */
class SpecularPoints: public ProjPatchOnSurface{

private:
    bool        m_insideout;        /**< plane direction for mirroring */
    bool        m_force;            /**< if true force the mirroring of points that belong to the symmetry plane. */
    darray4E    m_plane;            /**< reference plane */
    darray3E    m_origin;           /**<Origin of plane. */
    darray3E    m_normal;           /**<Normal of plane. */
    bool        m_implicit;         /**<True if an implicit definition of plane is set. */
    dmpvector1D * m_scalar;           /**< pointer to original float scalar data attached to original PointCloud*/
    dmpvecarr3E * m_vector;           /**< pointer to original float vector data attached to original PointCloud*/
    dmpvector1D   m_scalarMirrored;   /**< resulting mirrored Scalar Data*/
    dmpvecarr3E   m_vectorMirrored;   /**< resulting mirrored Vector Data*/
    MimmoSharedPointer<MimmoObject> m_pc; /**<original point cloud */

public:
    SpecularPoints();
    SpecularPoints(const bitpit::Config::Section & rootXML);
    ~SpecularPoints();

    SpecularPoints(const SpecularPoints & other);
    SpecularPoints& operator=(SpecularPoints other);
    void    buildPorts();

    dmpvector1D * getOriginalScalarData();
    dmpvecarr3E * getOriginalVectorData();

    dmpvector1D * getMirroredScalarData();
    dmpvecarr3E * getMirroredVectorData();
    dvecarr3E  getMirroredRawCoords();
    livector1D getMirroredLabels();
    MimmoSharedPointer<MimmoObject> getMirroredPointCloud();

    darray4E  getPlane();
    bool      isInsideOut();
    bool      isForced();

    void    setVectorData(dmpvecarr3E * vdata);
    void    setScalarData(dmpvector1D * data);
    void    setPointCloud(MimmoSharedPointer<MimmoObject> targetpatch);
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
    void swap(SpecularPoints & x) noexcept;

private:
    //interface method disabling
    void setPatch(MimmoSharedPointer<MimmoObject>geo){ProjPatchOnSurface::setPatch(geo);};
    MimmoSharedPointer<MimmoObject> getProjectedElement(){return ProjPatchOnSurface::getProjectedElement();};
    void setBuildSkdTree(bool build){ProjPrimitivesOnSurfaces::setBuildSkdTree(false);};
    void setWorkingOnTarget(bool flag){ProjPatchOnSurface::setWorkingOnTarget(flag);}
    bool isWorkingOnTarget(){return ProjPatchOnSurface::isWorkingOnTarget();}

};

REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_PLANE, MC_ARRAY4, MD_FLOAT,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__SPECULARPOINTS_HPP__)
REGISTER_PORT(M_AXIS, MC_ARRAY3, MD_FLOAT,__SPECULARPOINTS_HPP__)

REGISTER(BaseManipulation, SpecularPoints, "mimmo.SpecularPoints")
};

#endif /* __SPECULARPOINTS_HPP__ */
