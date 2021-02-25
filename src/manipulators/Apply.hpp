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
#ifndef __APPLYDEFORMATION_HPP__
#define __APPLYDEFORMATION_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{
/*!
 *  \class Apply
 *  \ingroup manipulators
 *  \brief Apply is the class that applies the deformation resulting from a manipulation object to the geometry.
 *
 *    Apply is derived from BaseManipulation class. It uses the base member m_geometry to apply
 *    the result of the parent manipulator to the target mimmo object.
 *    The deformation displacements have to be passed to the input member of base class through
 *    a pin linking or set by the user, i.e. one has to use setInput method to set
 *    the displacements to be applied to the geometry.
 *    The displacements can be isotropically scaled by setting a scaling factor.
 *    The displacements can be passed as scalar field. In this case the scalar field is used as magnitude
 *    of a displacements field with direction along the normal of the surface on each vertex.
 *    The Apply block allows to pass a scalar filter field; the filter field is applied to the
 *    displacements field before the geometry deformation.
 *
 *    After the execution of an object Apply, the original geometry will be modified.
 *    The resulting deformation field is provided as output of the block.
 *
 * \n
 * Ports available in Apply Class :
 *
 *    =========================================================
 *
     | Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GDISPLS | setInput          | (MC_SCALAR,MD_MPVECARR3FLOAT_) |
     | M_SCALARFIELD | setScalarInput    | (MC_SCALAR,MD_MPVECFLOAT_) |
     | M_FILTER | setFilter    | (MC_SCALAR,MD_MPVECFLOAT_) |
     | M_GEOM    | setGeometry       | (MC_SCALAR,MD_MIMMO_) |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>  |
     | M_GEOM          | getGeometry              | (MC_SCALAR,MD_MIMMO_) |
     | M_GDISPLS       | getOutput          | (MC_SCALAR,MD_MPVECARR3FLOAT_) |
     | M_LONGFIELD     | getAnnotatedVertices     | (MC_SCALAR,MD_MPVECLONG_) |
     | M_LONGFIELD2    | getAnnotatedCells        | (MC_SCALAR,MD_MPVECLONG_) |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.Apply</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 *
 * Proper of the class:
   - <B>Annotation</B>             : 0/1 boolean, if 1 provide list of cell/vertices into deformation as annotation list;
 * - <B>CellsAnnotationName</B>    : (string) define label to mark annotated cells (list of cell-ids involved into deformation);
 * - <B>VerticesAnnotationName</B> : (string) define label to mark annotated vertices (list of vertex-ids involved into deformation);
   - <B>AnnotationThreshold</B>    : double val >0, threshold over which all mesh elements with deformation field over it are annotated;

 *

 * Geometry and one of the Inputs have to be mandatorily passed through port.
 *
 */
class Apply: public BaseManipulation{

public:

    Apply();
    Apply(const bitpit::Config::Section & rootXML);

    ~Apply();

    Apply(const Apply & other);

    void buildPorts();

    void setInput(dmpvecarr3E *input);
    void setScalarInput(dmpvector1D* input);
    void setFilter(dmpvector1D* input);

    void setScaling(double alpha);

    void setAnnotation(bool activate);
    void setCellsAnnotationName(const std::string & label);
    void setVerticesAnnotationName(const std::string & label);
    void setAnnotationThreshold(double threshold);

    void execute();

    dmpvecarr3E * getOutput();
    MimmoPiercedVector<long> * getAnnotatedVertices();
    MimmoPiercedVector<long> * getAnnotatedCells();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

    dmpvecarr3E    m_input; /**< storing vector fields of floats */
    dmpvector1D    m_scalarinput; /**< storing scalar fields of floats */
    dmpvector1D    m_filter; /**< storing filter multiplying to the deformation field */
    double		   m_factor; /**< scaling factor of deformation field. */
    dmpvecarr3E    m_output; /**< storing vector fields of floats for output (input modified with filter and factor)*/

    bool m_annotation; /**< boolean to activate annotation providing*/
    std::string m_annCellLabel; /**< label marking cell annotation*/
    std::string m_annVertexLabel; /**< label marking cell annotation*/
    double m_annotationThres; /**< annotation threshold*/
    MimmoPiercedVector<long> m_cellAnnotation; /**< structure to store the list of cells involved into deformation*/
    MimmoPiercedVector<long> m_vertexAnnotation; /**< structure to store the list of vertices involved into deformation*/


    void swap(Apply & x) noexcept;
    void    checkInput();

};

REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_, __APPLYDEFORMATION_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __APPLYDEFORMATION_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_, __APPLYDEFORMATION_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __APPLYDEFORMATION_HPP__)
REGISTER_PORT(M_LONGFIELD, MC_SCALAR, MD_MPVECLONG_, __APPLYDEFORMATION_HPP__)
REGISTER_PORT(M_LONGFIELD2, MC_SCALAR, MD_MPVECLONG_, __APPLYDEFORMATION_HPP__)

REGISTER(BaseManipulation, Apply, "mimmo.Apply")

};

#endif /* __APPLYDEFORMATION_HPP__ */
