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
 * an input field linked to another geometry.
 * Three methods are available:
 * - <B>ID  = 0</B> : extraction of data on vertices with the same IDs between target and input geometries;
 * - <B>PID = 1</B> : extraction of data on vertices of cells with the PIDs found in the target geometry;
 * - <B>MAPPING = 2<B> : extraction of data on vertices of cells identified by a proximity mapping.
 *
 * ExtractField is an abstract class. To use its features take a look to its specializations,
 *  here presented as derived class, ExtractScalarField and ExtractVectorField.
 *
 * Ports available in ExtractField Class :
 *
 *    =========================================================
 *
     |                 Port Input     |||                                                     |
     |-------|------------|------------------------------------|-----------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 99    | M_GEOM     | setGeometry                        | (SCALAR, MIMMO_)            |


    |            Port Output         |||             |
     |-------|------------|------------------------------------|-----------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as "mimmo.Extract<Scalar/Vector>Fields"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>ExtractMode</B>: activate extraction mode by ID, PID or Mapping;
 *
 * Geometries and fields have to be mandatorily passed through port.
 *
 */
class ExtractField: public BaseManipulation{
protected:
    ExtractMode m_mode; /*< Extraction mode.*/

public:
    ExtractField();
    virtual ~ExtractField();

    ExtractField(const ExtractField & other);
    ExtractField & operator=(const ExtractField & other);

    void buildPorts();

    void        setGeometry(MimmoObject * geo);
    void        setMode(ExtractMode mode);
    void        setMode(int mode);

    void         clear();
    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

private:
    /*!
     * Pure virtual method
     */
    virtual bool extract() = 0;

};

/*!
 * \class ExtractScalarField
 * \ingroup geohandlers
 * \brief ExtractScalarField is specialized derived class of ExtractField to extract a
 *         scalar field of doubles.
 * 
 * Ports available in ExtractScalarField Class :
 *
 *    =========================================================
 *
     |                 Port Input   |||                              |
     |-------|--------------|--------------------|----------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 18    | M_SCALARFIELD| setField           | (MPVECTOR, FLOAT)|


     |            Port Output   |||                                        |
     |-------|-----------|-------------------|--------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 18    | M_SCALARFIELD  | getExtractedField     | (MPVECTOR, FLOAT)       |


  Inherited from ExtractField

     |                 Port Input   |||                                                         |
     |-------|------------|------------------------------------|-----------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 99    | M_GEOM     | setGeometry                        | (SCALAR, MIMMO_)            |

     |            Port Output  |||                               |
     |-------|-----------|------------------------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

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

    void buildPorts();
    dmpvector1D     getExtractedField();
    void     setField(dmpvector1D field);

    void clear();

    void     plotOptionalResults();

private:
    bool extract();

};

/*!
 *  \class ExtractVectorField
 *    \brief ExtractVectorField is specialized derived class of ExtractField to extract a
 *         scalar field of array<double,3>.
 * 
 * Ports available in ExtractVectorField Class :
 *
 *    =========================================================

    |------------------------------------------------------------------|
    |                 Port Input                                       |
     |-------|--------------|-------------------|-----------------------|
     |PortID | PortType     | variable/function | DataType              |
     |-------|--------------|-------------------|-----------------------|
      | 19    | M_VECTORFIELD     | setField     | (MPVECARR3, FLOAT)        |


     |-----------------------------------------------------------------------|
     |            Port Output                                                |
     |-------|-----------|-------------------|-------------------------------|
     |PortID | PortType  | variable/function | DataType                      |
     |-------|-----------|-------------------|-------------------------------|
      | 19    | M_VECTORFIELD   | getExtractedField     | (MPVECARR3, FLOAT)        |
     |-------|-----------|-------------------|-------------------------------|


   Inherited from ExtractField

     |---------------------------------------------------------------------------------------|
     |                 Port Input                                                            |
     |-------|------------|------------------------------------|-----------------------------|
     |PortID | PortType   | variable/function                  | DataType                    |
     |-------|------------|------------------------------------|-----------------------------|
     | 99    | M_GEOM     | setGeometry                        | (SCALAR, MIMMO_)            |


     |--------------------------------------------------------|-----------------------|
     |            Port Output                                 |                       |
     |-------|-----------|------------------------------------|-----------------------|
     |PortID | PortType  | variable/function                  | DataType              |
     |-------|-----------|------------------------------------|-----------------------|
     |-------|-----------|------------------------------------|-----------------------|

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

    void buildPorts();
    dmpvecarr3E     getExtractedField();
    void     setField(dmpvecarr3E field);

    void clear();

    void     plotOptionalResults();

private:
    bool extract();

};

REGISTER(BaseManipulation, ExtractScalarField, "mimmo.ExtractScalarField")
REGISTER(BaseManipulation, ExtractVectorField, "mimmo.ExtractVectorField")
};

#endif /* __EXTRACTFIELDS_HPP__ */
