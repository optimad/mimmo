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

#ifndef __SURFACETRIANGULATOR_HPP__
#define __SURFACETRIANGULATOR_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *
 *  \class SurfaceTriangulator
 *  \ingroup geohandlers
 *  \brief Triangulate a target MimmoObject non-homogeneous and/or non-triangular surface mesh.
 *
 * Class force a given 3D surface tessellation to be a full homogeneous triangular one.
 * By default class create a new indipendent copy of the target geometry and triangulate it.
 * A workOnTarget option can be activated if the User want to apply modifications directly
 * on the target mesh.
 *
 * Ports available in SurfaceTriangulator Class :
 *
 *    =========================================================
 *
    |                   Port Input       ||                               |
    |----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry        | (MC_SCALAR, MD_MIMMO_)      |


    |             Port Output     ||                                      |
    |----------------|---------------------|--------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | getGeometry        | (MC_SCALAR, MD_MIMMO_)      |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.SurfaceTriangulator</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : plot optional results in execution;
 * - <B>OutputPlot</B> : path to store optional results.
 *
 * Proper of the class:
 * - <B>WorkOnTarget</B>  : 0 -create a new homogeneously triangulated and independent copy of target mesh, 1 - modify directly the target mesh;
 *
 * Geometry has to be mandatorily passed through ports.
 *
 */
class SurfaceTriangulator: public mimmo::BaseManipulation {

private:

    bool m_workOnTarget;                      /**< true apply on target, false create an indipendent copy */
    mimmo::MimmoSharedPointer<MimmoObject> m_intPatch;     /**< internal patch */

public:

    SurfaceTriangulator();
    SurfaceTriangulator(const bitpit::Config::Section & rootXML);
    virtual ~SurfaceTriangulator();
    SurfaceTriangulator(const SurfaceTriangulator & other);
    SurfaceTriangulator & operator=(SurfaceTriangulator other);

    void    buildPorts();

    //get-set methods
    mimmo::MimmoSharedPointer<MimmoObject>    getGeometry();
    bool                            isWorkingOnTarget();

    void            setGeometry(mimmo::MimmoSharedPointer<MimmoObject> geo);
    void            setWorkOnTarget(bool flag = false);

    //execute
    void        execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
    virtual void plotOptionalResults();

protected:
    void    swap(SurfaceTriangulator & x) noexcept;
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__SURFACETRIANGULATOR_HPP__)
REGISTER(BaseManipulation, SurfaceTriangulator, "mimmo.SurfaceTriangulator")

};

#endif /* __SURFACETRIANGULATOR_HPP__ */
