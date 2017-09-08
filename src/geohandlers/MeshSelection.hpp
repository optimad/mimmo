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

#ifndef __MESHSELECTION_HPP__
#define __MESHSELECTION_HPP__

#include "BaseManipulation.hpp"
#include "MimmoObject.hpp"
#include "BasicShapes.hpp"
#include "MimmoGeometry.hpp"
#include "SkdTreeUtils.hpp"

#include <memory>
#include <utility>
namespace mimmo{

/*!
 * \enum SelectionType
 * \ingroup geohandlers
 * \brief Enum class for choiche of method to select sub-patch
 * of a tessellated mesh.
 */
enum class SelectionType{
    UNDEFINED    = 0,
            BOX            = 1,
            CYLINDER    = 2,
            SPHERE        = 3,
            MAPPING        = 4,
            PID            = 5,
            BOXwSCALAR  = 11
};

/*!
 * \class GenericSelection
 * \ingroup geohandlers
 * \brief Abstract Interface for selection classes
 * 
 * Class/BaseManipulation Object managing selection of sub-patches of a 3D open 
 * unstructured surface/volume mesh.
 *
 * Ports available in GenericSelection Class :
 *
 *    =========================================================
 *
     |                   Port Input       |||                               |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output     |||                                      |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *    =========================================================
 *
 */
class GenericSelection: public mimmo::BaseManipulation {

protected:

    SelectionType                   m_type;      /**< Type of enum class SelectionType for selection method */
    std::unique_ptr<MimmoObject>    m_subpatch;  /**< Pointer to result sub-patch */
    int                             m_topo;      /**< 1 = surface (default value), 2 = volume, 3 = points cloud */
    bool                            m_dual;      /**< False selects w/ current set up, true gets its "negative". False is default. */
public:

    GenericSelection();
    virtual ~GenericSelection();
    GenericSelection(const GenericSelection & other);
    GenericSelection & operator=(const GenericSelection & other);

    void    buildPorts();

    SelectionType    whichMethod();
    virtual void             setGeometry(MimmoObject *);
    void             setDual(bool flag=false);

    const MimmoObject*    getPatch()const;
    MimmoObject    *        getPatch();
    bool                isDual();

    livector1D    constrainedBoundary();

    void        execute();

    virtual void plotOptionalResults();

protected:
    /*!
     * Extract selection from target geometry
     */
    virtual livector1D extractSelection() = 0;
    void swap(GenericSelection &x) noexcept;
};


/*!
 * \class SelectionByBox
 * \ingroup geohandlers
 * \brief Selection through volume box primitive.
 * 
 * Selection Object managing selection of sub-patches of a 3D open 
 * unstructured surface mesh. Select all simplex contained in a box.
 * 
 * Ports available in SelectionByBox Class :
 *
 *    =========================================================

     |                   Port Input  |||                                       |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 20    | M_POINT        | setOrigin           | (ARRAY3, FLOAT)       |
     | 22    | M_AXES         | setRefSystem        | (ARR3ARR3, FLOAT)     |
     | 23    | M_SPAN         | setSpan             | (ARRAY3, FLOAT)       |


     |             Port Output        |||                                      |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *    ===============================================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionByBox</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 
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
class SelectionByBox: public GenericSelection, public mimmo::Cube {

public:
    SelectionByBox();
    SelectionByBox(const bitpit::Config::Section & rootXML);
    SelectionByBox(darray3E origin, darray3E span, MimmoObject * target);
    virtual ~SelectionByBox();
    SelectionByBox(const SelectionByBox & other);
    SelectionByBox & operator=(SelectionByBox other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    livector1D extractSelection();
    void swap(SelectionByBox &) noexcept;
};

/*!
 * \class SelectionByCylinder
 * \ingroup geohandlers
 * \brief Selection through cylinder primitive.
 * 
 * Selection Object managing selection of sub-patches of a 3D open 
 * unstructured surface mesh. Select all simplex contained in a cylinder.
 *
 * Ports available in SelectionByCylinder Class :
 *
 *    =========================================================

     |                   Port Input  |||                                       |
     |-------|----------------|----------------------|----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 20    | M_POINT        | setOrigin            | (ARRAY3, FLOAT)      |
     | 22    | M_AXES         | setRefSystem         | (ARR3ARR3, FLOAT)    |
     | 23    | M_SPAN         | setSpan              | (ARRAY3, FLOAT)      |
     | 25    | M_INFLIMITS    | setInfLimits         | (ARRAY3, FLOAT)      |

     |             Port Output        |||                                      |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionByCylinder</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
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
class SelectionByCylinder: public GenericSelection, public mimmo::Cylinder {

public:
    SelectionByCylinder();
    SelectionByCylinder(const bitpit::Config::Section & rootXML);
    SelectionByCylinder(darray3E origin, darray3E span, double infLimTheta, darray3E mainAxis, MimmoObject * target);
    virtual ~SelectionByCylinder();
    SelectionByCylinder(const SelectionByCylinder & other);
    SelectionByCylinder & operator=(SelectionByCylinder other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    livector1D extractSelection();
    void swap(SelectionByCylinder &) noexcept;
};

/*!
 * \class SelectionBySphere
 * \ingroup geohandlers
 * \brief Selection through sphere primitive.
 * 
 * Selection Object managing selection of sub-patches of a 3D open 
 * unstructured surface mesh. Select all simplex contained in a sphere.
 *
 * Ports available in SelectionBySphere Class :
 *
 *    =========================================================
 *
     |                   Port Input  |||                                       |
     |-------|----------------|----------------------|----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 20    | M_POINT        | setOrigin            | (ARRAY3, FLOAT)      |
     | 22    | M_AXES         | setRefSystem         | (ARR3ARR3, FLOAT)    |
     | 23    | M_SPAN         | setSpan              | (ARRAY3, FLOAT)      |
     | 25    | M_INFLIMITS    | setInfLimits         | (ARRAY3, FLOAT)      |

     |             Port Output    |||                                          |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionBySphere</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
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
class SelectionBySphere: public GenericSelection, public mimmo::Sphere {

public:
    SelectionBySphere();
    SelectionBySphere(const bitpit::Config::Section & rootXML);
    SelectionBySphere(darray3E origin, darray3E span, double infLimTheta, double infLimPhi, MimmoObject * target);
    virtual ~SelectionBySphere();
    SelectionBySphere(const SelectionBySphere & other);
    SelectionBySphere & operator=(SelectionBySphere other);

    void buildPorts();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    livector1D extractSelection();
    void swap(SelectionBySphere &) noexcept;
};

/*!
 * 
 * \class SelectionByMapping
 * \ingroup geohandlers
 * \brief Selection mapping external surfaces on the target mesh.
 * 
 * Selection Object managing selection of sub-patches of a 3D open 
 * unstructured surface mesh. Extract portion of mesh in common between 
 * a target geometry and a second one, provided externally. Extraction criterium
 * is based on euclidean nearness, within a prescribed tolerance.
 *
 * Ports available in SelectionByMapping Class :
 *
 *    =========================================================
 *
     |                   Port Input     |||                                                     |
     |-------|----------------|--------------------------|----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 98    | M_GEOM2        | addMappingGeometry       | (SCALAR, MIMMO_)     |

     |             Port Output       |||                                       |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *    =========================================================
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionByMapping</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Topology</B>: number indentifying topology of tesselated mesh. 1-surfaces, 2-voume. no other types are supported;
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual;
 * - <B>Tolerance</B>: proximity threshold to activate mapping;
 * - <B>Files</B>: list of external files to map on the target surface: \n
 *                <tt>\<Files\> \n
 *                      \<file0\> \n
 *                          \<fullpath> absolute path to your file with extension \<fullpath\> \n
 *                          \<tag\> tag extension as supported format STL,NAS,STVTU,etc...\<tag\> \n
 *                      \<file0\> \n
 *                      \<file1\> \n
 *                          \<fullpath> absolute path to your file with extension \<fullpath\> \n
 *                          \<tag\> tag extension as supported format STL,NAS,STVTU,etc...\<tag\> \n
 *                      \<file1\> \n
 *                      ... 
 *                  \</Files\> </tt> \n
 *
 * Geometry has to be mandatorily passed through port.
 */
class SelectionByMapping: public GenericSelection {

private:
    double                                  m_tolerance;     /**< tolerance for proximity detection*/
    std::unordered_map<std::string, int>    m_geolist;      /**< list of file for geometrical proximity check*/
    std::unordered_set<MimmoObject*>        m_mimmolist;    /**< list of mimmo objects for geometrical proximity check*/
    std::vector<std::set< int > >           m_allowedType;  /**< list of FileType actually allowed for the target geometry type*/
public:
    SelectionByMapping(int topo = 1);
    SelectionByMapping(const bitpit::Config::Section & rootXML);
    SelectionByMapping(std::unordered_map<std::string, int> & geolist, MimmoObject * target, double tolerance);
    virtual ~SelectionByMapping();
    SelectionByMapping(const SelectionByMapping & other);
    SelectionByMapping & operator=(SelectionByMapping other);

    void    buildPorts();

    double    getTolerance();
    void     setTolerance(double tol=1.e-8);

    void     setGeometry(MimmoObject * geometry);

    const std::unordered_map<std::string, int> &     getFiles() const;
    void    setFiles(std::unordered_map<std::string,int> );
    void    addFile(std::pair<std::string,int> );
    void     removeFile(std::string);
    void     removeFiles();
    void    addMappingGeometry(MimmoObject * obj);
    void    removeMappingGeometries();

    void clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    livector1D extractSelection();
    void swap(SelectionByMapping &) noexcept;

private:
    livector1D getProximity(std::pair<std::string, int> val);
    livector1D getProximity(MimmoObject * obj);
    svector1D extractInfo(std::string);
};

/*!
 * \class SelectionByPID
 * \ingroup geohandlers
 * \brief Selection using target mesh Part Identifiers.
 * 
 * Selection Object managing selection of sub-patches of a 3D open 
 * unstructured surface mesh. Extract portion of mesh getting its PID in 
 * target geometry.
 *
 * Ports available in SelectionByPID Class :
 *
 *    =========================================================
 *
     |                   Port Input   |||                                      |
     |-------|----------------|----------------------|----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORSI     | setPID               | (VECTOR, SHORT)      |
     | 35    | M_VALUESI      | setPID               | (SCALAR, SHORT)      |

     |             Port Output     |||                                         |
     |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *    =========================================================
 *
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionByMapping</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Dual</B>: boolean to get straight what given by selection method or its exact dual;
 * - <B>nPID</B>: number of PID to be selected relative to target geometry; 
 * - <B>PID</B>: list of PID (separated by blank spaces) to be selected relative to target geometry; 
 *
 * Geometry has to be mandatorily passed through port.
 */
class SelectionByPID: public GenericSelection {

private:
    std::unordered_set<short>       m_setPID;    /**< list of pid given by the user */
    std::unordered_map<short,bool> m_activePID; /**< list of pid available in geometry, to be flagged as active or not*/

public:
    SelectionByPID();
    SelectionByPID(const bitpit::Config::Section & rootXML);
    SelectionByPID(shivector1D & pidlist, MimmoObject * target);
    virtual ~SelectionByPID();
    SelectionByPID(const SelectionByPID & other);
    SelectionByPID & operator=(SelectionByPID other);

    void    buildPorts();

    shivector1D getPID();
    shivector1D    getActivePID(bool active= true);

    void     setGeometry(MimmoObject * );
    void     setPID(short i = -1);
    void     setPID(shivector1D);

    void     removePID(short i = -1);
    void     removePID(shivector1D);

    void     syncPIDList();
    void     clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

protected:
    livector1D extractSelection();
    void swap(SelectionByPID &) noexcept;

};

/*!
 * \class SelectionByBoxWithScalar
 * \ingroup geohandlers
 * \brief Selection through volume box primitive.
 * 
 * Selection Object managing selection of sub-patches
 * of a 3D unstructured surface mesh.
 * Select all simplex contained in a box and extract
 * a scalar field if it is present.
 *
 * Ports available in SelectionByBoxWithScalar Class :
 *
 *  =========================================================
 *
   |                   Port Input      |||                                   |
   |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 18    | M_SCALARFIELD  | setField            | (MPVECTOR, FLOAT)       |


   |             Port Output    |||                                          |
   |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 18    | M_SCALARFIELD  | getField            | (MPVECTOR, FLOAT)       |

   Inherited from SelectionByBox

   |                   Port Input  |||                                       |
   |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 20    | M_POINT        | setOrigin           | (ARRAY3, FLOAT)       |
   | 22    | M_AXES         | setRefSystem        | (ARR3ARR3, FLOAT)     |
   | 23    | M_SPAN         | setSpan             | (ARRAY3, FLOAT)       |


   |             Port Output       |||                                       |
   |-------|----------------|---------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     Inherited from GenericSelection

     |                   Port Input   |||                                   |
     |-------|----------------|----------------------|-------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 32    | M_VALUEB       | setDual              | (SCALAR, BOOL)    |
     | 99    | M_GEOM         | setGeometry          | (SCALAR, MIMMO_)  |


     |             Port Output          |||                                 |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 17    | M_VECTORLI     | constrainedBoundary | (VECTOR, LONG)     |
     | 99    | M_GEOM         | getPatch            | (SCALAR, MIMMO_)   |

 *  ===============================================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectionByBoxWithScalar</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 *Inherited from SelectionByBox:
 
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
 * Geometry and fields have to be mandatorily passed/recovered through ports.
 */
class SelectionByBoxWithScalar: public SelectionByBox{
protected:
    dmpvector1D     m_field;          /**<Scalar field attached to the patch
                                        (related to the whole patch before execution,
                                        related to the selected patch after execution).*/

public:
    SelectionByBoxWithScalar();
    SelectionByBoxWithScalar(const bitpit::Config::Section & rootXML);
    SelectionByBoxWithScalar(darray3E origin, darray3E span, MimmoObject * target);
    virtual ~SelectionByBoxWithScalar();
    SelectionByBoxWithScalar(const SelectionByBoxWithScalar & other);
    SelectionByBoxWithScalar & operator=(SelectionByBoxWithScalar other);

    void buildPorts();

    void clear();

    void        setField(dmpvector1D);
    dmpvector1D   getField();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="" );
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="" );

    void plotOptionalResults();
protected:
    void swap(SelectionByBoxWithScalar &) noexcept;
    
};


REGISTER(BaseManipulation, SelectionByBox,"mimmo.SelectionByBox")
REGISTER(BaseManipulation, SelectionByCylinder, "mimmo.SelectionByCylinder")
REGISTER(BaseManipulation, SelectionBySphere,"mimmo.SelectionBySphere")
REGISTER(BaseManipulation, SelectionByMapping, "mimmo.SelectionByMapping")
REGISTER(BaseManipulation, SelectionByPID, "mimmo.SelectionByPID")
REGISTER(BaseManipulation, SelectionByBoxWithScalar, "mimmo.SelectionByBoxWithScalar")
};

/*!
 * \}
 */


#endif /* __MESHSELECTION_HPP__ */
