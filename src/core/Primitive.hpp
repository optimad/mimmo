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
#ifndef __PRIMITIVE_HPP__
#define __PRIMITIVE_HPP__

#include "BaseManipulation.hpp"
#include "BasicShapes.hpp"
#include "BasicMeshes.hpp"

namespace mimmo{

/*!
 * \class Primitive
 * \ingroup core
 * \brief Primitive object generation.
 *
 *
 * Basically, it builds an elemental 3D shape.
 * The available shapes are: box, sphere, cylinder, wedge.
 * A non-complete primitive can be defined by using custom span dimensions.
 *
 *
 * Ports available in Primitive class:
 *
 * =========================================================
 *
 *
 *  |Port Input  | | |
    |-|-|-|
    |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | M_DIMENSION | setDimension                          | (MC_ARRAY3, MD_INT)         |
    | M_INFLIMITS | setInfLimits                          | (MC_ARRAY3, MD_FLOAT)       |
    | M_AXES      | setRefSystem                          | (MC_ARR3ARR3, MD_FLOAT)     |
    | M_SPAN      | setSpan                               | (MC_ARRAY3, MD_FLOAT)       |
    | M_POINT     | setOrigin                             | (MC_ARRAY3, MD_FLOAT)       |
    | M_SHAPE     | setShape(mimmo::ShapeType)            | (MC_SCALAR, MD_SHAPET)      |
    | M_SHAPEI    | setShape(int)                         | (MC_SCALAR, MD_INT)         |


    | Port Output | | |
    |-|-|-|
    |<B>PortType</B>|<B>variable/function</B>|<B>DataType</B>|
    | M_COPYSHAPE | getShape          | (MC_SCALAR, MD_SHAPE_)    |

* =========================================================
*
* The xml available parameters, sections and subsections  are the following:
*
*Inherited from BaseManipulation:
* - <B>ClassName</B>: name of the class as <tt>mimmo.Primitive</tt>;
* - <B>Priority</B>: uint marking priority in multi-chain execution;
* - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
* - <B>OutputPlot</B>: target directory for optional results writing.
*
*Proper of the class:
* - <B>Shape</B>: type of basic shape for your lattice. Available choice are CUBE, CYLINDER,SPHERE,WEDGE
* - <B>Origin</B>: 3D point marking the shape barycenter
* - <B>Span</B>: span dimensions of your shape (width-height-depth for CUBE,
                baseRadius-azimuthalspan-height for CYLINDER,
                radius-azimuthalspan-polarspan for SPHERE,
                triangle width-triangle height-depth for WEDGE)
* - <B>RefSystem</B>: axes of current shape reference system. written in XML as:\n\n
            <tt> <B>\<RefSystem\> </B> \n
            &nbsp;&nbsp;&nbsp;<B>\<axis0\></B> 1.0 0.0 0.0 <B>\</axis0\></B> \n
            &nbsp;&nbsp;&nbsp;<B>\<axis1\></B> 0.0 1.0 0.0 <B>\</axis1\></B> \n
            &nbsp;&nbsp;&nbsp;<B>\<axis2\></B> 0.0 0.0 1.0 <B>\</axis2\></B> \n
            <B>\</RefSystem\></B> </tt>\n\n
* - <B>InfLimits</B>: inferior limits for shape coordinates (meaningful only for CYLINDER AND SPHERE curvilinear coordinates)
*
*/
class Primitive: public BaseManipulation, public UStructMesh {

public:
    Primitive();
    Primitive(const bitpit::Config::Section & rootXML);
    virtual ~Primitive();

    //copy operators/constructors
    Primitive(const Primitive & other);

    void        buildPorts();
    void        clearPrimitive();

    void        execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

    darray3E                getSpacing() = delete;
    iarray3E                getDimension() = delete;
    darray3E                getLocalCCell(int) = delete;
    darray3E                getLocalCCell(int, int, int) = delete;
    darray3E                getLocalPoint(int) = delete;
    darray3E                getLocalPoint(int, int, int) = delete;
    darray3E                getGlobalCCell(int) = delete;
    darray3E                getGlobalCCell(int, int, int) = delete;
    darray3E                getGlobalPoint(int) = delete;
    darray3E                getGlobalPoint(int, int, int) = delete;
    ivector1D               getCellNeighs(int) = delete;
    ivector1D               getCellNeighs(int, int, int) = delete;
    dvecarr3E               getLocalCoords() = delete;
    dvecarr3E               getGlobalCoords() = delete;
    dvecarr3E               getGlobalCellCentroids() = delete;
    dvecarr3E               getLocalCellCentroids() = delete;
    void    setDimension(ivector1D dim) = delete;
    void    setDimension(iarray3E dim) = delete;
    void    locateCellByPoint(darray3E & point, int &i, int &j, int &k) = delete;
    void    locateCellByPoint(dvector1D & point, int &i, int &j, int &k) = delete;
    int     accessCellIndex(int i, int j, int k) = delete;
    void    accessCellIndex(int N_, int &i, int &j, int &k) = delete;
    int     accessPointIndex(int i, int j, int k) = delete;
    void    accessPointIndex(int N_, int &i, int &j, int &k) = delete;
    double      interpolateCellData(darray3E & point, dvector1D & celldata) = delete;
    int         interpolateCellData(darray3E & point, ivector1D & celldata) = delete;
    darray3E    interpolateCellData(darray3E & point, dvecarr3E & celldata) = delete;
    double      interpolatePointData(darray3E & point, dvector1D & pointdata) = delete;
    int         interpolatePointData(darray3E & point, ivector1D & pointdata) = delete;
    darray3E    interpolatePointData(darray3E & point, dvecarr3E & pointdata) = delete;
    void        plotCloud( std::string & , std::string, int , bool, const ivector1D & labels, dvecarr3E * extPoints=nullptr) = delete;
    void        plotCloudScalar(std::string, std::string , int, bool, dvector1D & data) = delete;
    void        plotGrid(std::string &, std::string , int, bool, const ivector1D & labels, dvecarr3E * extPoints=nullptr) = delete;
    void        plotGridScalar(std::string, std::string , int, bool, dvector1D & data) = delete;
    bool        isBuilt() = delete;

protected:
    void            swap(Primitive & ) noexcept;
//    virtual void    plotOptionalResults();

    void        resizeMesh() = delete;
    void        destroyNodalStructure() = delete;
    void        reshapeNodalStructure() = delete;


};

REGISTER_PORT(M_DIMENSION, MC_ARRAY3, MD_INT, __PRIMITIVE_HPP__)
REGISTER_PORT(M_INFLIMITS, MC_ARRAY3, MD_FLOAT, __PRIMITIVE_HPP__)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT, __PRIMITIVE_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT, __PRIMITIVE_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT, __PRIMITIVE_HPP__)
REGISTER_PORT(M_SHAPE, MC_SCALAR, MD_SHAPET, __PRIMITIVE_HPP__)
REGISTER_PORT(M_COPYSHAPE, MC_SCALAR, MD_SHAPE_, __PRIMITIVE_HPP__)
REGISTER_PORT(M_SHAPEI, MC_SCALAR, MD_INT, __PRIMITIVE_HPP__)

REGISTER(BaseManipulation, Primitive, "mimmo.Primitive")

}

#endif /* __PRIMITIVE_HPP__ */
