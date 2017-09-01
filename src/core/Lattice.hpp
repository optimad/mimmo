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
#ifndef __LATTICE_HPP__
#define __LATTICE_HPP__

#include "BasicMeshes.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class Lattice
 * \ingroup core
 * \brief Structured 3D Cartesian Mesh.
 * 
 *
 * Basically, it builds an elemental 3D shape (box, sphere, cylinder or part of them) around the
 * geometry and set a structured cartesian mesh of control points on it (lattice).
 * No displacements for control points and NO NURBS parameters for FFD are present
 * in this structure, only geometrical information are stored in the object.
 *
 * Ports available in Lattice class: 
 *
 * =========================================================
 *  
 *  
 *  |Port Input  | | | |
    |-|-|-|-|
    |<B>PortID</B>|<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | 99    | M_GEOM      | m_geometry                            | (SCALAR, MIMMO_)      |
    | 24    | M_DIMENSION | setDimension                          | (ARRAY3, INT)         |
    | 25    | M_INFLIMITS | setInfLimits                          | (ARRAY3, FLOAT)       |
    | 22    | M_AXES      | setRefSystem                          | (ARR3ARR3, FLOAT)     |
    | 23    | M_SPAN      | setSpan                               | (ARRAY3, FLOAT)       |
    | 20    | M_POINT     | setOrigin                             | (ARRAY3, FLOAT)       |
    | 26    | M_SHAPE     | setShape(mimmo::ShapeType)            | (SCALAR, SHAPET)      |
    | 27    | M_COPYSHAPE | setShape(const BasicShape*  )         | (SCALAR, SHAPE_)      |
    | 28    | M_SHAPEI    | setShape(int)                         | (SCALAR, INT)         |



    | Port Output | | | |
    |-|-|-|-|
    |<B>PortID</B>|<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | 1     | M_GLOBAL    | getGlobalCoords   | (VECARR3, FLOAT)    |
    | 2     | M_LOCAL     | getLocalCoords    | (VECARR3, FLOAT)    |
    | 20    | M_POINT     | getOrigin         | (ARRAY3, FLOAT)     |
    | 22    | M_AXES      | getRefSystem      | (ARR3ARR3, FLOAT)   |
    | 25    | M_INFLIMITS | getInfLimits      | (ARRAY3, FLOAT)     |
    | 23    | M_SPAN      | getSpan           | (ARRAY3, FLOAT)     |
    | 24    | M_DIMENSION | getDimension      | (ARRAY3, INT)       |
    | 27    | M_COPYSHAPE | getShape          | (SCALAR, SHAPE_)    |
    | 99    | M_GEOM      | getGeometry       | (SCALAR, MIMMO_)    |

* =========================================================
* 
* The xml available parameters, sections and subsections  are the following:
* 
* Inherited from BaseManipulation:
* - <B>ClassName</B>: name of the class as "mimmo.Lattice"
* - <B>Priority</B>: uint marking priority in multi-chain execution;
* - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
* - <B>OutputPlot</B>: target directory for optional results writing. 
*
*  Proper of the class:
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
* Geometry has to be mandatorily passed through port.
*/
class Lattice: public BaseManipulation, public UStructMesh {

protected:
    int         m_np;           /**< Number of control nodes.*/
    ivector1D   m_intMapDOF;    /**< Map of grid nodes -> degrees of freedom of lattice */
    
public:
    Lattice();
    Lattice(const bitpit::Config::Section & rootXML);
    virtual ~Lattice();

    //copy operators/constructors
    Lattice(const Lattice & other);

    void        buildPorts();
    void        clearLattice();
    int         getNNodes();
    dvecarr3E   getGlobalCoords();
    dvecarr3E   getLocalCoords();

    int         accessDOFFromGrid(int index);
    int         accessGridFromDOF(int index);
    ivector1D   accessDOFFromGrid(ivector1D gNindex);
    ivector1D   accessGridFromDOF(ivector1D dofIndex);

    void        plotGrid(std::string directory, std::string filename, int counter, bool binary);
    void        plotCloud(std::string directory, std::string filename, int counter, bool binary);

    virtual     void build();
    void        execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void            resizeMapDof();
    virtual void    plotOptionalResults();

private:
    int             reduceDimToDOF(int,int,int, bvector1D &info);

};

REGISTER(BaseManipulation, Lattice, "mimmo.Lattice")

}

#endif /* __LATTICE_HPP__ */
