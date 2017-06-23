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

#include "MimmoObject.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class SwitchField
 * \ingroup geohandlers
 * \brief SwitchField is an abstract executable block class capable of
 *         switchting a field from a list of fields.
 *
 *
 * SwitchField takes as input the target geometry used to indentify the field to choose from
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
 * - <B>ClassName</B>: name of the class as "mimmo.Switch<Scalar/Vector>Fields"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Mapping</B>: boolen 0/1 to force the research by mapping
 *
 * Geometries and fields have to be mandatorily passed through port.
 *
 */
class SwitchField: public BaseManipulation{
protected:
    bool m_mapping; /*< If true force the second research by mapping.*/

public:
    SwitchField();
    virtual ~SwitchField();

    SwitchField(const SwitchField & other);
    SwitchField & operator=(const SwitchField & other);

    void buildPorts();

    void        setGeometry(MimmoObject * geo);
    void        setMapping(bool flag = true);

    void         clear();
    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

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
 *         scalar field of doubles.
 * 
 * Ports available in SwitchScalarField Class :
 *
 *    =========================================================
 *
     |                 Port Input   |||                              |
     |-------|--------------|--------------------|----------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 60    | M_VECFIELDS| setFields           | (VECTOR, MPVECFLOAT)|
     | 19    | M_SCALARFIELD| addField           | (MPVECTOR, FLOAT)|


     |            Port Output   |||                                        |
     |-------|-----------|-------------------|--------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 19    | M_SCALARFIELD  | getSwitchedField     | (MPVECTOR, FLOAT)       |


  Inherited from SwitchField

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

class SwitchScalarField: public SwitchField{
private:
    vector<dmpvector1D> m_fields;   /**<Input fields to be switch. */
    dmpvector1D m_result;               /**<Result switch fields. */

public:
    SwitchScalarField();
    SwitchScalarField(const bitpit::Config::Section & rootXMl);
    virtual ~SwitchScalarField();

    void buildPorts();
    dmpvector1D     getSwitchedField();
    void     setFields(vector<dmpvector1D> fields);
    void     addField(dmpvector1D field);

    void clear();

    void     plotOptionalResults();

private:
    bool mswitch();

};

/*!
 *  \class SwitchVectorField
 *    \brief SwitchVectorField is specialized derived class of SwitchField to switch a
 *         scalar field of array<double,3>.
 * 
 * Ports available in SwitchVectorField Class :
 *
 *    =========================================================
 * ~~~
 *    |------------------------------------------------------------------|
 *    |                 Port Input                                       |
 *    |-------|--------------|-------------------|-----------------------|
 *    |PortID | PortType     | variable/function | DataType              |
 *    |-------|--------------|-------------------|-----------------------|
      | 61    | M_VECGDISPLS | setFields    | (VECTOR, MPVECARR3)       |
      | 11    | M_GDISPL     | addField     | (MPVECARR3, FLOAT)        |
 *
 *
 *    |-----------------------------------------------------------------------|
 *    |            Port Output                                                |
 *    |-------|-----------|-------------------|-------------------------------|
 *    |PortID | PortType  | variable/function | DataType                      |
 *    |-------|-----------|-------------------|-------------------------------|
      | 11    | M_GDISPL   | getSwitchedField     | (MPVECARR3, FLOAT)        |
 *    |-------|-----------|-------------------|-------------------------------|
 * 
 * 
 *  Inherited from SwitchField
 * 
 *    |---------------------------------------------------------------------------------------|
 *    |                 Port Input                                                            |
 *    |-------|------------|------------------------------------|-----------------------------|
 *    |PortID | PortType   | variable/function                  | DataType                    |
 *    |-------|------------|------------------------------------|-----------------------------|
 *    | 99    | M_GEOM     | setGeometry                        | (SCALAR, MIMMO_)            |
 *
 *
 *    |--------------------------------------------------------|-----------------------|
 *    |            Port Output                                 |                       |
 *    |-------|-----------|------------------------------------|-----------------------|
 *    |PortID | PortType  | variable/function                  | DataType              |
 *    |-------|-----------|------------------------------------|-----------------------|
 *    |-------|-----------|------------------------------------|-----------------------|
 * ~~~
 *    =========================================================
 *
 */
class SwitchVectorField: public SwitchField{
private:
    vector<dmpvecarr3E>  m_fields;   /**<Input fields to be switch. */
    dmpvecarr3E m_result;                /**<Result switch fields. */

public:

    SwitchVectorField();
    SwitchVectorField(const bitpit::Config::Section & rootXMl);
    virtual ~SwitchVectorField();

    void buildPorts();
    dmpvecarr3E     getSwitchedField();
    void     setFields(vector<dmpvecarr3E> fields);
    void     addField(dmpvecarr3E field);

    void clear();

    void     plotOptionalResults();

private:
    bool mswitch();

};

REGISTER(BaseManipulation, SwitchScalarField, "mimmo.SwitchScalarField")
REGISTER(BaseManipulation, SwitchVectorField, "mimmo.SwitchVectorField")
};

#endif /* __SWITCHFIELDS_HPP__ */
