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
#ifndef __PROJECTCLOUD_HPP__
#define __PROJECTCLOUD_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *  \class ProjectCloud
 *  \ingroup utils
 *  \brief ProjectCloud is the class that project a 3D cloud of points on a target geometry
 *
 * \n
 * Ports available in ProjectCloud Class :
 *
 *  =========================================================

   |Port Input | | |
   |-|-|-|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_COORDS | setCoords         | (MC_VECARR3, MD_FLOAT)    |
   | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)    |
   | M_GEOM2  | setCoords         | (MC_SCALAR, MD_MIMMO_)    |


   |Port Output  | | |
   |-|-|-|
   | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>   |
   | M_COORDS | getProjectedCoords         | (MC_VECARR3, MD_FLOAT)    |
   | M_GEOM   | getProjectedCloud          | (MC_SCALAR, MD_MIMMO_)    |


 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ProjectCloud</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Geometry and points have to be mandatorily passed through port.
 *
 */

class ProjectCloud: public BaseManipulation{

protected:
    dvecarr3E            m_points;    /**<Coordinates of 3D points in the cloud.*/
    dvecarr3E            m_proj;     /**<Projected points coordinates.*/
    livector1D           m_labels;   /**<labels associated to 3D points in the cloud.*/

    MimmoSharedPointer<MimmoObject> m_internalPC; /**< internal PointCloud MimmoObject for I/O purposes*/

public:
    ProjectCloud();
    ProjectCloud(const bitpit::Config::Section & rootXML);
    ~ProjectCloud();

    ProjectCloud(const ProjectCloud & other);
    ProjectCloud& operator=(ProjectCloud other);

    void         buildPorts();

    dvecarr3E    getOriginalCoords(livector1D * labels = nullptr);
    dvecarr3E    getProjectedCoords();
    dvecarr3E    getProjectedCoords(livector1D * labels);
    MimmoSharedPointer<MimmoObject> getProjectedCloud();

    void     setCoords(dvecarr3E coords);
    void     setCoords(MimmoSharedPointer<MimmoObject> coordsPC);
    void     setGeometry(MimmoSharedPointer<MimmoObject> refSurface);
    void     execute();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name= "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    virtual void plotOptionalResults();
    void swap(ProjectCloud & x) noexcept;
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__PROJECTCLOUD_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__PROJECTCLOUD_HPP__)
REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__PROJECTCLOUD_HPP__)

REGISTER(BaseManipulation, ProjectCloud,"mimmo.ProjectCloud")

};

#endif // __PROJECTCLOUD_HPP__ 
