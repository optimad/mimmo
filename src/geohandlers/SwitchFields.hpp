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
#ifndef __SWITCHFIELDS_HPP__
#define __SWITCHFIELDS_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class SwitchField
 * \ingroup geohandlers
 * \brief SwitchField is an abstract executable block class capable of
 *         switching a field from a list of fields.
 *
 *
 * SwitchField takes as input the target geometry used to identify the related field to choose from
 * the input list. First it tries to switch the field by investigating the linked geometries
 * in the MimmoPiercedVector input fields and by extracting the first corresponding field.
 * In case of none of the input fields is linked to the target geometry
 * (and the block is set to do it) it performs a mapping
 * by testing if the geometry overlap with the geometries linked to the input fields;
 * for each element of the target geometry is chosen the first related element found in the list of the
 * input geometries.
 *
 *    SwitchField is an abstract class. To use its features take a look to its specializations,
 *  here presented as derived class, SwitchScalarField and SwitchVectorField.
 *
 *
 * Ports available in SwitchField Class :
 *
 *    =========================================================
 *
     |                 Port Input     ||                                                     |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |


     |            Port Output         ||             |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.Switch<...>Fields</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Mapping</B>: boolen 0/1 to force the research by mapping
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.

 * Geometries and fields have to be mandatorily passed through port.
 *
 */
class SwitchField: public BaseManipulation{
protected:
    bool m_mapping; /**< If true force the second research by mapping.*/
    double m_tol;  /**< Tolerance for extraction by patch.*/

public:
    SwitchField();
    virtual ~SwitchField();

    SwitchField(const SwitchField & other);

    void buildPorts();

    void        setGeometry(MimmoObject * geo);
    void        setMapping(bool flag = true);
    void        setTolerance(double tol);

    void         clear();
    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SwitchField &) noexcept;
private:
    /*!
     * Pure virtual method
     */

    virtual bool mswitch() = 0;
};

/*!
 * \class SwitchScalarField
 * \ingroup geohandlers
 * \brief SwitchScalarField is specialized derived class of SwitchField to switch a
 *         scalar field.
 *
 * Ports available in SwitchScalarField Class :
 *
 *    =========================================================

     |                 Port Input   ||                              |
     |--------------|--------------------|----------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VECSFIELDS| setFields           | (MC_VECTOR, MD_MPVECFLOAT_)|
     | M_SCALARFIELD| addField           | (MC_SCALAR, MD_MPVECFLOAT_)|


     |            Port Output   ||                                        |
     |-----------|-------------------|--------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | getSwitchedField     | (MC_SCALAR, MD_MPVECFLOAT_)       |


  Inherited from SwitchField

     |                 Port Input   ||                                                         |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

     |            Port Output  ||                               |
     |-----------|------------------------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SwitchScalarField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SwitchField:
 * - <B>Mapping</B>: boolen 0/1 to force the research by mapping.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 *
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Geometries and fields have to be mandatorily passed through port.
 */

class SwitchScalarField: public SwitchField{
private:
    MPVLocation m_loc;  /**< field data reference location */
    std::vector<dmpvector1D> m_fields;   /**<Input fields to be switch. */
    dmpvector1D m_result;               /**<Result switch fields. */

public:
    SwitchScalarField(MPVLocation loc = MPVLocation::POINT);
    SwitchScalarField(const bitpit::Config::Section & rootXMl);
    SwitchScalarField(const SwitchScalarField & other);
    virtual ~SwitchScalarField();

    void buildPorts();
    dmpvector1D*     getSwitchedField();
    void     setFields(std::vector<dmpvector1D *> fields);
    void     addField(dmpvector1D *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SwitchScalarField &) noexcept;

private:
    bool mswitch();

};

/*!
 *  \class SwitchVectorField
    \ingroup geohandlers
 *    \brief SwitchVectorField is specialized derived class of SwitchField to switch a
 *         vector field.
 *
 * Ports available in SwitchVectorField Class :
 *
 *    =========================================================


 |                 Port Input   ||                              |
 |--------------|--------------------|----------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECVFIELDS| setFields           | (MC_VECTOR, MD_MPVECARR3FLOAT_)|
 | M_VECTORFIELD| addField           | (MC_SCALAR, MD_MPVECARR3FLOAT_)|


 |            Port Output   ||                                        |
 |-----------|-------------------|--------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECTORFIELD  | getSwitchedField     | (MC_SCALAR, MD_MPVECARR3FLOAT_)       |


 Inherited from SwitchField

 |                 Port Input   ||                                                         |
 |------------|------------------------------------|-----------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 |            Port Output  ||                               |
 |-----------|------------------------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SwitchVectorField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SwitchField:
 * - <B>Mapping</B>: boolen 0/1 to force the research by mapping.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
*
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Geometries and fields have to be mandatorily passed through port.
 */
class SwitchVectorField: public SwitchField{
private:
    MPVLocation m_loc;  /**< field data reference location */
    std::vector<dmpvecarr3E>  m_fields;   /**<Input fields to be switch. */
    dmpvecarr3E m_result;                /**<Result switch fields. */

public:

    SwitchVectorField(MPVLocation loc = MPVLocation::POINT);
    SwitchVectorField(const bitpit::Config::Section & rootXMl);
    SwitchVectorField(const SwitchVectorField & other);
    virtual ~SwitchVectorField();

    void buildPorts();
    dmpvecarr3E*     getSwitchedField();
    void     setFields(std::vector<dmpvecarr3E *> fields);
    void     addField(dmpvecarr3E *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SwitchVectorField &) noexcept;

private:
    bool mswitch();

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __SWITCHFIELDS_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __SWITCHFIELDS_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_, __SWITCHFIELDS_HPP__)
REGISTER_PORT(M_VECSFIELDS, MC_VECTOR, MD_MPVECFLOAT_, __SWITCHFIELDS_HPP__)
REGISTER_PORT(M_VECVFIELDS, MC_VECTOR, MD_MPVECARR3FLOAT_, __SWITCHFIELDS_HPP__)

REGISTER(BaseManipulation, SwitchScalarField, "mimmo.SwitchScalarField")
REGISTER(BaseManipulation, SwitchVectorField, "mimmo.SwitchVectorField")
};

#endif /* __SWITCHFIELDS_HPP__ */
