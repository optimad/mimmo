/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
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
#ifndef __APPLYFILTER_HPP__
#define __APPLYFILTER_HPP__

#include "Apply.hpp"

namespace mimmo{
/*!
 *  \class ApplyFilterFilter
 *  \ingroup manipulators
 *  \brief ApplyFilterFilter is the class that applies a filter field to a deformation field defined on a geometry.
 *
 *    ApplyFilterFilter is derived from ApplyFilter class. It uses an input deformation field defined
 *    on a linked geometry to apply to it an input filter field defined over the same geometry.
 *    The scalar filter field is applied by a simple modulation of the deformation field.
 *    The displacements can be isotropically scaled by setting a scaling factor.
 *    The displacements can be passed as scalar field. In this case the scalar field is used as magnitude
 *    of a displacements field with direction along the normal of the surface on each vertex.
 *    If a geoemtry is passed through the port the resulting deformation field will be linked to the input geometry,
 *    even if the input fields are not linked to it. If missing data are present, they are set to zero values.
 *    The ApplyFilterFilter block applies the filter to the
 *    displacements field without perform a geometry deformation.
 *    After the execution, an object ApplyFilterFilter provides as result the filtered deformation field.
 *
 * \n
 * Ports available in ApplyFilterFilter Class :
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
     | M_GDISPLS | geOutput          | (MC_SCALAR,MD_MPVECARR3FLOAT_) |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.ApplyFilter</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 *
 * Proper of the class:

 *
 *
 */
class ApplyFilter: public Apply{

public:

    ApplyFilter();
    ApplyFilter(const bitpit::Config::Section & rootXML);

    ~ApplyFilter();

    ApplyFilter(const ApplyFilter & other);

    void buildPorts();

    using Apply::setScaling;

    void setAnnotation(bool activate) = delete;
    void setCellsAnnotationName(const std::string & label) = delete;
    void setVerticesAnnotationName(const std::string & label) = delete;
    void setAnnotationThreshold(double threshold) = delete;

    void execute();

    using Apply::getOutput;

    MimmoPiercedVector<long> * getAnnotatedVertices() = delete;
    MimmoPiercedVector<long> * getAnnotatedCells() = delete;

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

    using   Apply::m_input;         /**< storing vector fields of floats */
    using   Apply::m_scalarinput;   /**< storing scalar fields of floats */
    using   Apply::m_filter;        /**< storing filter multiplying to the deformation field */
    using   Apply::m_factor;        /**< scaling factor of deformation field. */
    using   Apply::m_output;        /**< storing vector fields of floats for output (input modified with filter and factor)*/

    void swap(ApplyFilter & x) noexcept;
    void    checkInput();

};

REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_, __APPLYFILTER_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __APPLYFILTER_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_, __APPLYFILTER_HPP__)

REGISTER(BaseManipulation, ApplyFilter, "mimmo.ApplyFilter")

};

#endif /* __APPLYFILTER_HPP__ */
