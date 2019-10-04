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
#ifndef __EXTRACTFIELDS_HPP__
#define __EXTRACTFIELDS_HPP__

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{
/*!
 * \enum ExtractMode
 * \ingroup geohandlers
 * \brief Modes available to extract fields.See class ExtractField documentation.
 */
enum class ExtractMode{
    ID = 1 /**< Extract via ID*/,
    PID = 2 /**< Extract via PID*/,
    MAPPING = 3 /**< Extract by proximity mapping.*/
};

/*!
 * \class ExtractField
 * \ingroup geohandlers
 * \brief ExtractField is an abstract executable block class capable of
 *         extracting a portion of an input field related to a geometry
 *         sub-portion of the geometry linked to the input field .
 *
 *
 * ExtractField takes as input the target geometry used to extract a portion of
 * an input field linked to the same or another geometry.
 * Reference data location (POINT, CELLS or INTERFACES) of the input field is used to extract data
 * on vertices, cells or interfaces of the target geometry.
 * Three methods are available:
 * - <B>ID  = 1</B> : extraction of data on reference Location with the same IDs between target and input field linked geometry.
 *                    Return an extracted field referred to the target geometry, with the same data location of the input field.
 * - <B>PID = 2</B> : extraction of data on reference Location within common PIDs between target geometry and field linked geometry;
 *                    Return an extracted field referred to the input field linked geometry, with the same data location of the input field.
 * - <B>MAPPING = 3</B> : extraction of data on reference Location identified by a proximity mapping between target geometry and field linked geometry;
 *                       Return an extracted field referred to the target geometry,with the same data location of the input field.
 *
 * ExtractField is an abstract class. To use its features take a look to its specializations,
 *  here presented as derived class, ExtractScalarField and ExtractVectorField.
 *
 * Ports available in ExtractField Class :
 *
 *    =========================================================
 *
     |                 Port Input     ||                                                     |
     |------------|------------------------------------|-----------------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |


     |            Port Output         ||             |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.Extract< Scalar/Vector >Field"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>ExtractMode</B>: activate extraction mode by ID(1), PID(2) or Mapping(3);
 * - <B>Tolerance</B>: tolerance for extraction by patch, valid only in Mapping mode.
 *
 * Geometries and fields have to be mandatorily passed through port.
 *
 */
class ExtractField: public BaseManipulation{
protected:
    ExtractMode m_mode; /**< Extraction mode.*/
    double      m_tol;  /**< Tolerance for extraction by patch.*/

public:
    ExtractField();
    virtual ~ExtractField();

    ExtractField(const ExtractField & other);
    ExtractField& operator=(const ExtractField & other);

    void buildPorts();

    void        setGeometry(MimmoObject * geo);
    void        setMode(ExtractMode mode);
    void        setMode(int mode);
    void        setTolerance(double tol);

    ExtractMode getMode();
    double      getTolerance();

    void         clear();
    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    /*!
     * Pure virtual method
     */
    virtual bool extract() = 0;

protected:
    void swap(ExtractField & x) noexcept;

};

/*!
 * \class ExtractScalarField
 * \ingroup geohandlers
 * \brief ExtractScalarField is specialized derived class of ExtractField to extract a
 *         scalar field of doubles.
 *
 * Parameters:
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.Extract< Scalar/Vector >Field"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Inherited from ExtractField:
 * - <B>ExtractMode</B>: activate extraction mode by ID(1), PID(2) or Mapping(3);
 * - <B>Tolerance</B>: tolerance for extraction by patch, valid only in Mapping mode.

 * Ports available in ExtractScalarField Class :
 *
 *    =========================================================
 *
     |                 Port Input   ||                              |
     |--------------|--------------------|----------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD| setField           | (MC_SCALAR, MD_MPVECFLOAT_)|


     |            Port Output   ||                                        |
     |-----------|-------------------|--------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | getExtractedField     | (MC_SCALAR, MD_MPVECFLOAT_)       |


  Inherited from ExtractField

     |                 Port Input   ||                                                         |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

     |            Port Output  ||                               |
     |-----------|------------------------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 *
 */

class ExtractScalarField: public ExtractField{
private:
    dmpvector1D m_field;             /**<Input field to be extract. */
    dmpvector1D m_result;               /**<Result extract field. */

public:
    ExtractScalarField();
    ExtractScalarField(const bitpit::Config::Section & rootXMl);
    virtual ~ExtractScalarField();

    ExtractScalarField(const ExtractScalarField & other);
    ExtractScalarField& operator=(const ExtractScalarField & other);

    void buildPorts();
    dmpvector1D* getExtractedField();
    void     setField(dmpvector1D *field);

    dmpvector1D getOriginalField();

    void clear();

    void     plotOptionalResults();
    bool     extract();

protected:
    void swap(ExtractScalarField & x) noexcept;

private:

    void extractID(mimmo::MPVLocation loc);
    void extractPID(mimmo::MPVLocation loc);
    void extractMapping(mimmo::MPVLocation loc);

};

/*!
 *  \class ExtractVectorField
 *    \brief ExtractVectorField is specialized derived class of ExtractField to extract a
 *         scalar field of array<double,3>.
 *
 * Parameters:
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.Extract< Scalar/Vector >Field"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Inherited from ExtractField:
 * - <B>ExtractMode</B>: activate extraction mode by ID(1), PID(2) or Mapping(3);
 * - <B>Tolerance</B>: tolerance for extraction by patch, valid only in Mapping mode.

 * Ports available in ExtractVectorField Class :
 *
 *    =========================================================

 |                 Port Input   ||                              |
 |--------------|--------------------|----------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECTORFIELD| setField           | (MC_SCALAR, MD_MPVECARR3FLOAT_)|


 |            Port Output   ||                                        |
 |-----------|-------------------|--------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECTORFIELD  | getExtractedField     | (MC_SCALAR, MD_MPVECARR3FLOAT_)       |


 Inherited from ExtractField

 |                 Port Input   ||                                                         |
 |------------|------------------------------------|-----------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 |            Port Output  ||                               |
 |-----------|------------------------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


 *    =========================================================
 *
 */
class ExtractVectorField: public ExtractField{
private:
    dmpvecarr3E m_field;         /**<Input field to be extract. */
    dmpvecarr3E m_result;        /**<Result extract field. */

public:

    ExtractVectorField();
    ExtractVectorField(const bitpit::Config::Section & rootXMl);
    virtual ~ExtractVectorField();
    ExtractVectorField(const ExtractVectorField & other);
    ExtractVectorField& operator=(const ExtractVectorField & other);

    void buildPorts();
    dmpvecarr3E  *   getExtractedField();
    void     setField(dmpvecarr3E*field);

    dmpvecarr3E     getOriginalField();
    void clear();

    void     plotOptionalResults();
    bool     extract();

protected:
    void swap(ExtractVectorField & x) noexcept;

private:
    void extractID(mimmo::MPVLocation loc);
    void extractPID(mimmo::MPVLocation loc);
    void extractMapping(mimmo::MPVLocation loc);

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __EXTRACTFIELDS_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__EXTRACTFIELDS_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__EXTRACTFIELDS_HPP__)


REGISTER(BaseManipulation, ExtractScalarField, "mimmo.ExtractScalarField")
REGISTER(BaseManipulation, ExtractVectorField, "mimmo.ExtractVectorField")
};

#endif /* __EXTRACTFIELDS_HPP__ */
