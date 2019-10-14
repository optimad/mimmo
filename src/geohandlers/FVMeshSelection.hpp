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
 \ *---------------------------------------------------------------------------*/

#ifndef __FVMESHSELECTION_HPP__
#define __FVMESHSELECTION_HPP__

#include "BaseManipulation.hpp"
#include "MimmoObject.hpp"
#include "mimmo_common.hpp"

namespace mimmo{

/*!
 * \enum FVSelectionType
 * \ingroup geohandlers
 * \brief Enum class for choiche of method to select sub-patch of a volume mesh
 */
enum class FVSelectionType{
    UNDEFINED    = 0,
    BOX          = 1,
    CYLINDER     = 2,
    SPHERE       = 3
};

/*!
 * \class FVGenericSelection
 * \ingroup geohandlers
 * \brief Abstract Interface for selection classes related to Volume Unstructured Mesh
 *
 * Class/BaseManipulation Object managing selection of sub-patches of MimmoFvMesh data structures.
 *
 * Ports available in FVGenericSelection Class :
 *
 *    =========================================================
 *
     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VALUEB       | setDual              | (MC_SCALAR, MD_BOOL)    |
     | M_GEOM         | setGeometry          | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | setBoundary          | (MC_SCALAR, MD_MIMMO_)  |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM         | getVolumePatch        | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | getBoundaryPatch      | (MC_SCALAR, MD_MIMMO_)  |

 *    =========================================================
 *
 */
class FVGenericSelection: public mimmo::BaseManipulation {

protected:

    FVSelectionType                 m_type;      /**< Type of enum class SelectionType for selection method */
    std::unique_ptr<MimmoObject>    m_volpatch;  /**< Pointer to result volume sub-patch */
    std::unique_ptr<MimmoObject>    m_bndpatch;  /**< Pointer to result boundary sub-patch */
    int                             m_topo;      /**< 1 = volume (default value), 2 = surface */
    bool                            m_dual;      /**< False selects w/ current set up, true gets its "negative". False is default. */
    MimmoObject *                   m_bndgeometry; /**<target boundary geometry */

public:

    FVGenericSelection(int topo = 1);
    virtual ~FVGenericSelection();
    FVGenericSelection(const FVGenericSelection & other);
    FVGenericSelection & operator=(const FVGenericSelection & other);

    void    buildPorts();

    FVSelectionType    whichMethod();
    void               setGeometry(MimmoObject *);
    void               setBoundaryGeometry(MimmoObject *);
    void               setDual(bool flag=false);

    const MimmoObject*    getVolumePatch()const;
    const MimmoObject*    getBoundaryPatch()const;
    MimmoObject    *        getVolumePatch();
    MimmoObject    *        getBoundaryPatch();

    bool                isDual();

    void        execute();

    virtual void plotOptionalResults();
    bool checkCoherenceBulkBoundary();

protected:
    /*!
     * Extract selection from target geometry
     */
    virtual void extractSelection(livector1D &, livector1D &) = 0;
    void swap(FVGenericSelection &x) noexcept;
};


/*!
 * \class FVSelectionByBox
 * \ingroup geohandlers
 * \brief Selection through volume box primitive.
 *
 * Selection Object managing selection of sub-patches of a MimmoFvMesh Data structure.
 * Select all elements contained in a box.
 *
 * Ports available in FVSelectionByBox Class :
 *
 *    =========================================================

     |                   Port Input  ||                                       |
     |----------------|---------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT        | setOrigin           | (MC_ARRAY3, MD_FLOAT)       |
     | M_AXES         | setRefSystem        | (MC_ARR3ARR3, MD_FLOAT)     |
     | M_SPAN         | setSpan             | (MC_ARRAY3, MD_FLOAT)       |


     |             Port Output        ||                                      |
     |----------------|---------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from FVGenericSelection

     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VALUEB       | setDual              | (MC_SCALAR, MD_BOOL)    |
     | M_GEOM         | setGeometry          | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | setBoundary          | (MC_SCALAR, MD_MIMMO_)  |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM         | getVolumePatch        | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | getBoundaryPatch      | (MC_SCALAR, MD_MIMMO_)  |

 *    ===============================================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.FVSelectionByBox</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Topology</B> : 1 for volume mesh+surface boundary, 2 for surface mesh + 3Dcurve boundary. No other value allowed.
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual;
 * - <B>Origin</B>: array of 3 doubles identifying origin;
 * - <B>Span</B>: span of the box (width, height, depth);
 * - <B>RefSystem</B>: reference system of the box: \n
 *                <tt>\<RefSystem\> \n
 *                      \<axis0\> 1.0 0.0 0.0 \</axis0\> \n
 *                      \<axis1\> 0.0 1.0 0.0 \</axis1\> \n
 *                      \<axis2\> 0.0 0.0 1.0 \</axis2\> \n
 *                  \</RefSystem\> </tt> \n
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class FVSelectionByBox: public FVGenericSelection, public mimmo::Cube {

public:
    FVSelectionByBox(int topo = 1);
    FVSelectionByBox(const bitpit::Config::Section & rootXML);
    FVSelectionByBox(int topo, darray3E origin, darray3E span);
    virtual ~FVSelectionByBox();
    FVSelectionByBox(const FVSelectionByBox & other);
    FVSelectionByBox & operator=(FVSelectionByBox other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    void extractSelection(livector1D &, livector1D &);
    void swap(FVSelectionByBox &) noexcept;
};

/*!
 * \class FVSelectionByCylinder
 * \ingroup geohandlers
 * \brief Selection through cylinder primitive.
 *
 * Selection Object managing selection of sub-patches oof a MimmoFvMesh Data structure.
 * Select all elements contained in a cylinder.
 *
 * Ports available in FVSelectionByCylinder Class :
 *
 *    =========================================================

     |                   Port Input  ||                                       |
     |----------------|----------------------|----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT        | setOrigin            | (MC_ARRAY3, MD_FLOAT)      |
     | M_AXES         | setRefSystem         | (MC_ARR3ARR3, MD_FLOAT)    |
     | M_SPAN         | setSpan              | (MC_ARRAY3, MD_FLOAT)      |
     | M_INFLIMITS    | setInfLimits         | (MC_ARRAY3, MD_FLOAT)      |

     |             Port Output        ||                                      |
     |----------------|---------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from FVGenericSelection

     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VALUEB       | setDual              | (MC_SCALAR, MD_BOOL)    |
     | M_GEOM         | setGeometry          | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | setBoundary          | (MC_SCALAR, MD_MIMMO_)  |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM         | getVolumePatch        | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | getBoundaryPatch      | (MC_SCALAR, MD_MIMMO_)  |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :

 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.FVSelectionByCylinder</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Topology</B> : 1 for volume mesh+surface boundary, 2 for surface mesh + 3Dcurve boundary. No other value allowed.
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual;
 * - <B>Origin</B>: array of 3 doubles identifying origin of cylinder;
 * - <B>Span</B>: span of the cylinder (base radius, angular azimuthal width, height);
 * - <B>RefSystem</B>: reference system of the cylinder(axis2 along the cylinder's height): \n
 *                <tt>\<RefSystem\> \n
 *                      \<axis0\> 1.0 0.0 0.0 \</axis0\> \n
 *                      \<axis1\> 0.0 1.0 0.0 \</axis1\> \n
 *                      \<axis2\> 0.0 0.0 1.0 \</axis2\> \n
 *                  \</RefSystem\> </tt> \n
 * - <B>InfLimits</B>: set starting point for each cylindrical coordinate. Useful to assign different starting angular coordinate to azimuthal width.
 *
 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class FVSelectionByCylinder: public FVGenericSelection, public mimmo::Cylinder {

public:
    FVSelectionByCylinder(int topo = 1);
    FVSelectionByCylinder(const bitpit::Config::Section & rootXML);
    FVSelectionByCylinder(int topo, darray3E origin, darray3E span, double infLimTheta, darray3E mainAxis);
    virtual ~FVSelectionByCylinder();
    FVSelectionByCylinder(const FVSelectionByCylinder & other);
    FVSelectionByCylinder & operator=(FVSelectionByCylinder other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    void extractSelection(livector1D &, livector1D &);
    void swap(FVSelectionByCylinder &) noexcept;
};

/*!
 * \class FVSelectionBySphere
 * \ingroup geohandlers
 * \brief Selection through sphere primitive.
 *
 * Selection Object managing selection of sub-patches of a MimmoFvMesh Data structure.
 * Select all elements contained in a sphere.
 *
 * Ports available in FVSelectionBySphere Class :
 *
 *    =========================================================
 *
     |                   Port Input  ||                                       |
     |----------------|----------------------|----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_POINT        | setOrigin            | (MC_ARRAY3, MD_FLOAT)      |
     | M_AXES         | setRefSystem         | (MC_ARR3ARR3, MD_FLOAT)    |
     | M_SPAN         | setSpan              | (MC_ARRAY3, MD_FLOAT)      |
     | M_INFLIMITS    | setInfLimits         | (MC_ARRAY3, MD_FLOAT)      |

     |             Port Output    ||                                          |
     |----------------|---------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from FVGenericSelection

     |                   Port Input       ||                               |
     |----------------|----------------------|-------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VALUEB       | setDual              | (MC_SCALAR, MD_BOOL)    |
     | M_GEOM         | setGeometry          | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | setBoundary          | (MC_SCALAR, MD_MIMMO_)  |


     |             Port Output     ||                                      |
     |----------------|---------------------|--------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM         | getVolumePatch        | (MC_SCALAR, MD_MIMMO_)  |
     | M_GEOM2        | getBoundaryPatch      | (MC_SCALAR, MD_MIMMO_)  |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:

 * - <B>ClassName</B>: name of the class as <tt>mimmo.FVSelectionBySphere</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Topology</B> : 1 for volume mesh+surface boundary, 2 for surface mesh + 3Dcurve boundary. No other value allowed.
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual;
 * - <B>Origin</B>: array of 3 doubles identifying origin of sphere;
 * - <B>Span</B>: span of the sphere (radius, angular azimuthal width, angular polar width);
 * - <B>RefSystem</B>: reference system of the sphere: \n
 *                <tt>\<RefSystem\> \n
 *                      \<axis0\> 1.0 0.0 0.0 \</axis0\> \n
 *                      \<axis1\> 0.0 1.0 0.0 \</axis1\> \n
 *                      \<axis2\> 0.0 0.0 1.0 \</axis2\> \n
 *                  \</RefSystem\> </tt> \n
 * - <B>InfLimits</B>: set starting point for each spherical coordinate. Useful to assign different starting angular coordinate to azimuthal width/ polar width.
 *
 *
 * Geometry has to be mandatorily passed through port.

 *
 */
class FVSelectionBySphere: public FVGenericSelection, public mimmo::Sphere {

public:
    FVSelectionBySphere(int topo = 1);
    FVSelectionBySphere(const bitpit::Config::Section & rootXML);
    FVSelectionBySphere(int topo, darray3E origin, darray3E span, double infLimTheta, double infLimPhi);
    virtual ~FVSelectionBySphere();
    FVSelectionBySphere(const FVSelectionBySphere & other);
    FVSelectionBySphere & operator=(FVSelectionBySphere other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    void extractSelection(livector1D &, livector1D &);
    void swap(FVSelectionBySphere &) noexcept;
};


REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __FVMESHSELECTION_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_, __FVMESHSELECTION_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT, __FVMESHSELECTION_HPP__)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT, __FVMESHSELECTION_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT, __FVMESHSELECTION_HPP__)
REGISTER_PORT(M_INFLIMITS, MC_ARRAY3, MD_FLOAT, __FVMESHSELECTION_HPP__)


REGISTER(BaseManipulation, FVSelectionByBox,"mimmo.FVSelectionByBox")
REGISTER(BaseManipulation, FVSelectionByCylinder, "mimmo.FVSelectionByCylinder")
REGISTER(BaseManipulation, FVSelectionBySphere,"mimmo.FVSelectionBySphere")
};


#endif /* __FVMESHSELECTION_HPP__ */
