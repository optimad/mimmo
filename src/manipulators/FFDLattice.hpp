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
#ifndef __FFDLATTICE_HPP__
#define __FFDLATTICE_HPP__

#include "Lattice.hpp"

namespace mimmo{

/*!
 * \class FFDLattice
 * \ingroup manipulators
 * \brief Free Form Deformation of a 3D surface and point clouds, with structured lattice.
 *
 * Free Form deformation tool for 3D geometries (surface and point clouds). Basically, it builds an elemental 3D shape
 *  (box, sphere, cylinder or part of them) around the geometry and set a structured cartesian mesh of control
 *  points on it (lattice). Displacements of each control point is linked to the geometry inside
 *  the shape by means of a NURBS volumetric parameterization. Deformation will be applied only to
 *  those portion of geometry encased into the 3D shape.
 *
 * \n
 * Ports available in FFDLattice Class :
 *
 *    =========================================================


     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_DISPLS         | setDisplacements              | (MC_VECARR3, MD_FLOAT)    |
     | M_FILTER         | setFilter                     | (MC_SCALAR, MD_MPVECFLOAT_)     |
     | M_DEG            | setDegrees                    | (MC_ARRAY3, MD_INT)       |
     | M_NURBSWEIGHTS   | setNodalWeight                | (MC_VECTOR, MD_FLOAT)     |
     | M_NURBSCOORDTYPE | setCoordType                  | (MC_ARRAY3, MD_COORDT)    |


     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GDISPLS         | getDeformation    | (MC_SCALAR, MD_MPVECARR3FLOAT_)              |
     | M_DEG             | getDegrees        | (MC_ARRAY3, MD_INT)                 |
     | M_FILTER          | getFilter         | (MC_SCALAR, MD_MPVECFLOAT_)               |
     | M_NURBSWEIGHTS    | getWeights        | (MC_VECTOR, MD_FLOAT)               |
     | M_NURBSCOORDTYPE  | getCoordType      | (MC_ARRAY3, MD_COORDT)              |


      Inherited from Lattice :

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM      | m_geometry                            | (MC_SCALAR, MD_MIMMO_)      |
     | M_DIMENSION | setDimension                          | (MC_ARRAY3, MD_INT)         |
     | M_INFLIMITS | setInfLimits                          | (MC_ARRAY3, MD_FLOAT)       |
     | M_AXES      | setRefSystem                          | (MC_ARR3ARR3, MD_FLOAT)     |
     | M_SPAN      | setSpan                               | (MC_ARRAY3, MD_FLOAT)       |
     | M_POINT     | setOrigin                             | (MC_ARRAY3, MD_FLOAT)       |
     | M_SHAPE     | setShape(mimmo::ShapeType)            | (MC_SCALAR, MD_SHAPET)      |
     | M_COPYSHAPE | setShape(const BasicShape  )          | (MC_SCALAR, MD_SHAPE_       |
     | M_SHAPEI    | setShape(int)                         | (MC_SCALAR, MD_INT)         |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GLOBAL    | getGlobalCoords   | (MC_VECARR3, MD_FLOAT)    |
     | M_LOCAL     | getLocalCoords    | (MC_VECARR3, MD_FLOAT)    |
     | M_POINT     | getOrigin         | (MC_ARRAY3, MD_FLOAT)     |
     | M_AXES      | getRefSystem      | (MC_ARR3ARR3, MD_FLOAT)   |
     | M_INFLIMITS | getInfLimits      | (MC_ARRAY3, MD_FLOAT)     |
     | M_SPAN      | getSpan           | (MC_ARRAY3, MD_FLOAT)     |
     | M_DIMENSION | getDimension      | (MC_ARRAY3, MD_INT)       |
     | M_COPYSHAPE | getShape          | (MC_SCALAR, MD_SHAPE_)    |
     | M_GEOM      | getGeometry       | (MC_SCALAR, MD_MIMMO_)    |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.Lattice"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 *  Inherited from Lattice:
 * - <B>Shape</B>: type of basic shape for your lattice. Available choice are CUBE, CYLINDER,SPHERE
 * - <B>Origin</B>: 3D point marking the shape barycenter
 * - <B>Span</B>: span dimensions of your shape (width-height-depth for CUBE, baseRadius-azimuthalspan-height for CYLINDER, radius-azimuthalspan-polarspan for SPHERE)
 * - <B>RefSystem</B>: axes of current shape reference system. written in XML as:
 *                  <tt> \n \<RefSystem\> \n
 *                       \<axis0\> 1.0 0.0 0.0 \</axis0\> \n
 *                       \<axis1\> 0.0 1.0 0.0 \</axis1\> \n
 *                       \<axis2\> 0.0 0.0 1.0 \</axis2\> \n
 *                   \</RefSystem\> </tt>
 * - <B>InfLimits</B>: inferior limits for shape coordinates (meaningful only for CYLINDER AND SPHERE curvilinear coordinates)
 * - <B>Dimension</B>: number of nodes in each coordinate direction to get the structured lattice mesh
 *
 * Proper of the class:
 * - <B>CoordType</B>: Set Boundary conditions for each NURBS interpolant on their extrema. Available choice are <B>CLAMPED,SYMMETRIC,UNCLAMPED, PERIODIC</B>;
 * - <B>Degrees</B>: degrees for NURBS interpolant in each spatial direction;
 * - <B>DisplGlobal</B>:0/1 use local-shape/global x,y,z reference system to define displacements of lattice node;
 *

 *
 * Geometry, displacements field and filter field have to be mandatorily passed through port.
 */
class FFDLattice: public Lattice {

protected:
    iarray3E    m_deg;            /**< Nurbs curve degree for each of the possible 3 direction in space*/
    dvector2D    m_knots;        /**< Nurbs curve knots for each of the possible 3 direction in space*/
    ivector2D    m_mapEff;        /**< Nurbs map of theoretical node distribution */
    dvector1D    m_weights;        /**< Weights of each control node*/
    dvecarr3E    m_displ;        /**< Displacements of control nodes.*/
    ivector2D     m_mapNodes;        /**< Internal map to access node index w/ knots structure theoretical indexing */
    dmpvecarr3E    m_gdispl;        /**< Displacements of geometry vertex.*/
private:
    iarray3E    m_mapdeg;        /**< Map of curves degrees. Increasing order of curves degrees. */
    bool        m_globalDispl;     /**< Choose type of displacements passed to lattice TRUE/Global XYZ displacement, False/local shape ref sys*/
    std::unordered_map<int, double> m_collect_wg; /**< temporary collector of nodal weights passed as parameter. Nodal weight can be applied by build() method */
    dmpvector1D   m_filter;        /**< Filter scalar field defined on geometry nodes for displacements modulation*/
    bool         m_bfilter;        /**< Boolean to recognize if a filter field for for displacements modulation is set or not */

public:
    FFDLattice();
    FFDLattice(const bitpit::Config::Section & rootXML);
    virtual ~FFDLattice();

    //copy operators/constructors
    FFDLattice(const FFDLattice & other);
    FFDLattice & operator=(FFDLattice other);

    void buildPorts();

    //clean structure
    void         clearLattice();
    void         clearFilter();

    //internal methods
    ivector1D     getKnotsDimension();
    dvector1D   getWeights();
    void         returnKnotsStructure(dvector2D &, ivector2D &);
    void         returnKnotsStructure( int, dvector1D &, ivector1D &);
    dvecarr3E*     getDisplacements();
    dmpvector1D*   getFilter();
    dmpvecarr3E*  getDeformation();
    bool         isDisplGlobal();
    iarray3E    getDegrees();

    void        setDegrees(iarray3E curveDegrees);
    void         setDisplacements(dvecarr3E displacements);
    void         setDisplGlobal(bool flag);
    void         setLattice(darray3E & origin, darray3E & span, ShapeType, iarray3E & dimensions, iarray3E & degrees);
    void         setLattice(darray3E & origin, darray3E & span, ShapeType, dvector1D & spacing, iarray3E & degrees);
    void         setLattice(BasicShape *, iarray3E & dimensions,  iarray3E & degrees);
    void         setLattice(BasicShape *, dvector1D & spacing,iarray3E & degrees);

    void         setNodalWeight(double , int );
    void         setNodalWeight(double , int, int, int);
    void         setNodalWeight(dvector1D );

    void        setFilter(dmpvector1D * );

    //plotting wrappers
    void        plotGrid(std::string directory, std::string filename, int counter, bool binary, bool deformed);
    void        plotCloud(std::string directory, std::string filename, int counter, bool binary, bool deformed);


    //execute deformation methods
    void         execute();
    darray3E     apply(darray3E & point);
    dvecarr3E     apply(dvecarr3E * point);
    dvecarr3E     apply(livector1D & map);

    virtual bool         build();

    void     apply();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    darray3E    convertDisplToXYZ(darray3E &, int i);
    dvecarr3E    convertDisplToXYZ();
    void         resizeMapDof();
    virtual void    plotOptionalResults();
    void        swap(FFDLattice &) noexcept;
    void         checkFilter();

private:
    //Nurbs Evaluators
    darray3E    nurbsEvaluator(darray3E &);
    dvecarr3E    nurbsEvaluator(livector1D &);
    double        nurbsEvaluatorScalar(darray3E &, int);

    //Nurbs utilities
    dvector1D    basisITS0(int k, int pos, double coord);
    dvector1D    getNodeSpacing(int dir);

    //knots mantenaince utilities
    void         clearKnots();
    void         setKnotsStructure();
    void         setKnotsStructure(int dir, CoordType type);
    int          getKnotInterval(double, int);
    double         getKnotValue(int, int);
    int         getKnotIndex(int,int);
    int         getTheoreticalKnotIndex(int,int);

    //nodal displacement utility
    dvecarr3E    recoverFullGridDispl();
    dvector1D    recoverFullNodeWeights();
    void         setMapNodes(int ind);
    int          accessMapNodes(int,int,int);


    //dimension utilities
    void        orderDimension();
};

/*! Return real global index of a nodal displacement, given its position i,j,k in knots indexing logic*/
inline int FFDLattice::accessMapNodes(int i, int j, int k){
    return(accessPointIndex(m_mapNodes[0][i], m_mapNodes[1][j], m_mapNodes[2][k]));
};

REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT,__FFDLATTICE_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_,__FFDLATTICE_HPP__)
REGISTER_PORT(M_DEG, MC_ARRAY3, MD_INT,__FFDLATTICE_HPP__)
REGISTER_PORT(M_NURBSWEIGHTS, MC_VECTOR, MD_FLOAT,__FFDLATTICE_HPP__)
REGISTER_PORT(M_NURBSCOORDTYPE, MC_ARRAY3, MD_COORDT,__FFDLATTICE_HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_,__FFDLATTICE_HPP__)

REGISTER(BaseManipulation, FFDLattice, "mimmo.FFDLattice")
}

#endif /* __FFDLATTICE_HPP__ */
