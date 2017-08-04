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
 *  |Port Input  | | |
    |-|-|-|
    |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | M_GEOM      | m_geometry                            | (MC_SCALAR, MD_MIMMO_)      |
    | M_DIMENSION | setDimension                          | (MC_ARRAY3, MD_INT)         |
    | M_INFLIMITS | setInfLimits                          | (MC_ARRAY3, MD_FLOAT)       |
    | M_AXES      | setRefSystem                          | (MC_ARR3ARR3, MD_FLOAT)     |
    | M_SPAN      | setSpan                               | (MC_ARRAY3, MD_FLOAT)       |
    | M_POINT     | setOrigin                             | (MC_ARRAY3, MD_FLOAT)       |
    | M_SHAPE     | setShape(mimmo::ShapeType)            | (MC_SCALAR, MD_SHAPET)      |
    | M_COPYSHAPE | setShape(const BasicShape*  )         | (MC_SCALAR, MD_SHAPE_)      |
    | M_SHAPEI    | setShape(int)                         | (MC_SCALAR, MD_INT)         |



    | Port Output | | |
    |-|-|-|
    |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | M_GLOBAL    | getGlobalCoords   | (MC_VECARR3, MD_FLOAT)    |
    | M_LOCAL     | getLocalCoords    | (MC_VECARR3, MD_FLOAT)    |
    | M_POINT     | getOrigin         | (MC_ARRAY3, MD_FLOAT)     |
    | M_AXES      | getRefSystem      | (MC_ARR3ARR3, MD_FLOAT)   |
    | M_INFLIMITS | getInfLimits      | (MC_ARRAY3, MD_FLOAT)     |
    | M_SPAN      | getSpan           | (MC_ARRAY3, MD_FLOAT)     |
    | M_DIMENSION | getDimension      | (MC_ARRAY3, MD_INT)       |
    | M_COPYSHAPE | getShape          | (MC_SCALAR, MD_SHAPE_)    |
    | M_GEOM      | getGeometry       | (MC_SCALAR, MD_MIMMO_)    |

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
    void            swap(Lattice & ) noexcept;
    void            resizeMapDof();
    virtual void    plotOptionalResults();

private:
    int             reduceDimToDOF(int,int,int, bvector1D &info);

};

//Ports
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_)
REGISTER_PORT(M_DIMENSION, MC_ARRAY3, MD_INT)
REGISTER_PORT(M_INFLIMITS, MC_ARRAY3, MD_FLOAT)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT)
REGISTER_PORT(M_SHAPE, MC_SCALAR, MD_SHAPET)
REGISTER_PORT(M_COPYSHAPE, MC_SCALAR, MD_SHAPE_)
REGISTER_PORT(M_SHAPEI, MC_SCALAR, MD_INT)
REGISTER_PORT(M_GLOBAL, MC_VECARR3, MD_FLOAT)
REGISTER_PORT(M_LOCAL, MC_VECARR3, MD_FLOAT)


//ManipBlocks
REGISTER(BaseManipulation, Lattice, "mimmo.Lattice")

}

#endif /* __LATTICE_HPP__ */
