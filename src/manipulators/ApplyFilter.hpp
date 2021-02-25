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
#ifndef __APPLYFILTER_HPP__
#define __APPLYFILTER_HPP__

#include "Apply.hpp"

namespace mimmo{
/*!
 *  \class ApplyFilter
 *  \ingroup manipulators
 *  \brief ApplyFilter is a class that applies a filter field to a deformation field defined on a geometry.
 *
 *    ApplyFilter is custom derivation of Apply class, targeted to the manipulation/modulation
 *    of a displacements field by means of a custom filter field, both externally provided,
      both referring to a target geometry.
      After the execution, ApplyFilter provides as result the filtered deformation field.
      BEWARE: Differently from its base class, it does not apply the final displacement field to
      the geometry.To achieve it, use an Apply block.
 *
 * \n
 * Ports available in ApplyFilter Class :
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

    void execute();

    using Apply::getOutput;

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

private:
    void setAnnotation(bool activate) = delete;
    void setCellsAnnotationName(const std::string & label) = delete;
    void setVerticesAnnotationName(const std::string & label) = delete;
    void setAnnotationThreshold(double threshold) = delete;
    MimmoPiercedVector<long> * getAnnotatedVertices() = delete;
    MimmoPiercedVector<long> * getAnnotatedCells() = delete;

};

REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_, __APPLYFILTER_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __APPLYFILTER_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_, __APPLYFILTER_HPP__)

REGISTER(BaseManipulation, ApplyFilter, "mimmo.ApplyFilter")

};

#endif /* __APPLYFILTER_HPP__ */
