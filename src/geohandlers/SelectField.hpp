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
#ifndef __SelectFieldS_HPP__
#define __SelectFieldS_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \ingroup geohandlers
 * \brief Methods available for selecting a field.
 */
enum class SelectType{
			GEOMETRY = 0, /**< Select the field linked to an input geometry.
    					The first found field linked to the geometry is selected. */
     		NAME = 1, /**< Select the field by name. The first found field with the input name is selected. */
		    MAPPING = 2 /**< Select the elements of the result field by geometry mapping. For each element of the target geometry is chosen the first related element found in the list of the
 	 	 	 	 	 	 input geometries. Note: not allowed for point cloud objects.*/
};

/*!
 * \class SelectField
 * \ingroup geohandlers
 * \brief SelectField is an abstract executable block class capable of
 *         Selecting a field from a list of fields.
 *
 * SelectField can work with different selection mode:
 * - SelectType::GEOMETRY, SelectField takes as input the target geometry used to identify the related field to choose from
 * the input list. It tries to Select the field by investigating the linked geometries
 * in the MimmoPiercedVector input fields and by extracting the first corresponding field.
 * - SelectType::NAME, SelectField takes as input a field name used to identify the field to choose from
 * the input list. It tries to Select the field by extracting the first field with name equal to the input name.
 * - SelectType::MAPPING, SelectField performs a mapping by testing if the geometry overlap with the geometries
 * linked to the input fields; for each element of the target geometry is chosen the first related element
 * found in the list of the input geometries.
 *
 *    SelectField is an abstract class. To use its features take a look to its specializations,
 *  here presented as derived class, SelectScalarField and SelectVectorField.
 *
 *
 * Ports available in SelectField Class :
 *
 *    =========================================================
 *
     |                 Port Input     ||                                                     |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |
     | M_NAME     | setName		                        | (MC_SCALAR, MD_STRING)            |


     |            Port Output         ||             |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | getGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.Select<...>Fields</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>SelectType</B>: selection method. Default SelectType::NAME.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 * - <B>FieldName</B>: name of the field to be selected. Default "data".

 * Fields have to be mandatorily passed through port.
 *
 */
class SelectField: public BaseManipulation{
protected:
	SelectType m_mode; /**< Selection method to be used.*/
    double m_tol;  /**< Tolerance for extraction by patch.*/
    std::string m_fieldname; /**<Name of the field to be selected.
    						If SelectType is by geometry, after the execution m_fieldname is the name of the selcted field.
    						If SelectType is by mapping, the user can set the result field name by setFieldName method (or FieldName option in xml input). */

public:
    SelectField();
    virtual ~SelectField();

    SelectField(const SelectField & other);

    void buildPorts();

    void        setGeometry(MimmoObject * geo);
    void        setMode(SelectType mode);
    void        setMode(int mode);
    void        setFieldName(std::string fieldname);
    void        setTolerance(double tol);

    std::string	getFieldName();

    void         clear();
    void         execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SelectField &) noexcept;
private:
    /*!
     * Pure virtual method
     */

    virtual bool mSelect() = 0;
};

/*!
 * \class SelectScalarField
 * \ingroup geohandlers
 * \brief SelectScalarField is specialized derived class of SelectField to Select a
 *         scalar field.
 *
 * Ports available in SelectScalarField Class :
 *
 *    =========================================================

     |                 Port Input   ||                              |
     |--------------|--------------------|----------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VECSCALARFIELDS| setFields           | (MC_VECTOR, MD_MPVECFLOAT_)|
     | M_SCALARFIELD| addField           | (MC_SCALAR, MD_MPVECFLOAT_)|


     |            Port Output   ||                                        |
     |-----------|-------------------|--------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | getSelectedField     | (MC_SCALAR, MD_MPVECFLOAT_)       |


  Inherited from SelectField

     |                 Port Input   ||                                                         |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |
     | M_NAME     | setName		                        | (MC_SCALAR, MD_STRING)            |

     |            Port Output  ||                               |
     |-----------|------------------------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM     | getGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectScalarField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SelectField:
 * - <B>SelectType</B>: selection method. Default SelectType::NAME.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 * - <B>FieldName</B>: name of the field to be selected. Default "data".
 *
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Fields have to be mandatorily passed through port.
 */

class SelectScalarField: public SelectField{
private:
    MPVLocation m_loc;  /**< field data reference location */
    std::vector<dmpvector1D> m_fields;   /**<Input fields to be Select. */
    dmpvector1D m_result;               /**<Result Select fields. */

public:
    SelectScalarField(MPVLocation loc = MPVLocation::POINT);
    SelectScalarField(const bitpit::Config::Section & rootXMl);
    SelectScalarField(const SelectScalarField & other);
    virtual ~SelectScalarField();

    void buildPorts();
    dmpvector1D*     getSelectedField();
    void     setFields(std::vector<dmpvector1D *> fields);
    void     addField(dmpvector1D *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SelectScalarField &) noexcept;

private:
    bool mSelect();

};



/*!
 *  \class SelectVectorField
    \ingroup geohandlers
 *    \brief SelectVectorField is specialized derived class of SelectField to Select a
 *         vector field.
 *
 * Ports available in SelectVectorField Class :
 *
 *    =========================================================


 |                 Port Input   ||                              |
 |--------------|--------------------|----------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECVECTORFIELDS| setFields           | (MC_VECTOR, MD_MPVECARR3FLOAT_)|
 | M_VECTORFIELD| addField           | (MC_SCALAR, MD_MPVECARR3FLOAT_)|


 |            Port Output   ||                                        |
 |-----------|-------------------|--------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECTORFIELD  | getSelectedField     | (MC_SCALAR, MD_MPVECARR3FLOAT_)       |


 Inherited from SelectField

 |                 Port Input   ||                                                         |
 |------------|------------------------------------|-----------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |
 | M_NAME     | setName		                        | (MC_SCALAR, MD_STRING)            |

 |            Port Output  ||                               |
 |-----------|------------------------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | getGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectVectorField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SelectField:
 * - <B>SelectType</B>: selection method. Default SelectType::NAME.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 * - <B>FieldName</B>: name of the field to be selected. Default "data".
 *
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Fields have to be mandatorily passed through port.
 */
class SelectVectorField: public SelectField{
private:
    MPVLocation m_loc;  /**< field data reference location */
    std::vector<dmpvecarr3E>  m_fields;   /**<Input fields to be Select. */
    dmpvecarr3E m_result;                /**<Result Select fields. */

public:

    SelectVectorField(MPVLocation loc = MPVLocation::POINT);
    SelectVectorField(const bitpit::Config::Section & rootXMl);
    SelectVectorField(const SelectVectorField & other);
    virtual ~SelectVectorField();

    void buildPorts();
    dmpvecarr3E*     getSelectedField();
    void     setFields(std::vector<dmpvecarr3E *> fields);
    void     addField(dmpvecarr3E *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SelectVectorField &) noexcept;

private:
    bool mSelect();

};



/*!
 *  \class SelectLongField
    \ingroup geohandlers
 *    \brief SelectLongField is specialized derived class of SelectField to Select a
 *         scalar field of long data.
 *
 * Ports available in SelectLongField Class :
 *
 *    =========================================================


 |                 Port Input   ||                              |
 |--------------|--------------------|----------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECLONGFIELDS| setFields           | (MC_VECTOR, MD_MPVECLONG_)|
 | M_LONGFIELD| addField           | (MC_SCALAR, MD_MPVECLONG_)|


 |            Port Output   ||                                        |
 |-----------|-------------------|--------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_LONGFIELD  | getSelectedField     | (MC_SCALAR, MD_MPVECLONG_)       |


 Inherited from SelectField

 |                 Port Input   ||                                                         |
 |------------|------------------------------------|-----------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |
 | M_NAME     | setName		                        | (MC_SCALAR, MD_STRING)            |

 |            Port Output  ||                               |
 |-----------|------------------------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | getGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectLongField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SelectField:
 * - <B>SelectType</B>: selection method. Default SelectType::NAME.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 * - <B>FieldName</B>: name of the field to be selected. Default "data".
 *
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Fields have to be mandatorily passed through port.
 */
class SelectLongField: public SelectField{
private:
    MPVLocation 							m_loc;  	/**< field data reference location */
    std::vector<MimmoPiercedVector<long> >	m_fields;	/**<Input fields to be Select. */
    MimmoPiercedVector<long>				m_result;   /**<Result Select fields. */

public:

    SelectLongField(MPVLocation loc = MPVLocation::POINT);
    SelectLongField(const bitpit::Config::Section & rootXMl);
    SelectLongField(const SelectLongField & other);
    virtual ~SelectLongField();

    void buildPorts();
    MimmoPiercedVector<long>*     getSelectedField();
    void     setFields(std::vector<MimmoPiercedVector<long> *> fields);
    void     addField(MimmoPiercedVector<long> *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SelectLongField &) noexcept;

private:
    bool mSelect();

};



/*!
 *  \class SelectStringField
    \ingroup geohandlers
 *    \brief SelectStringField is specialized derived class of SelectField to Select a
 *         scalar field of string data.
 *
 * Ports available in SelectStringField Class :
 *
 *    =========================================================


 |                 Port Input   ||                              |
 |--------------|--------------------|----------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_VECSTRINGFIELDS| setFields           | (MC_VECTOR, MD_MPVECSTRING_)|
 | M_STRINGFIELD| addField           | (MC_SCALAR, MD_MPVECSTRING_)|


 |            Port Output   ||                                        |
 |-----------|-------------------|--------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_STRINGFIELD  | getSelectedField     | (MC_SCALAR, MD_MPVECSTRING_)       |


 Inherited from SelectField

 |                 Port Input   ||                                                         |
 |------------|------------------------------------|-----------------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | setGeometry                        | (MC_SCALAR, MD_MIMMO_)            |
 | M_NAME     | setName		                        | (MC_SCALAR, MD_STRING)            |

 |            Port Output  ||                               |
 |-----------|------------------------------------|-----------------------|
 | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 | M_GEOM     | getGeometry                        | (MC_SCALAR, MD_MIMMO_)            |

 *    =========================================================
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.SelectStringField</tt>
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.

 * Inherited from SelectField:
 * - <B>SelectType</B>: selection method. Default SelectType::NAME.
 * - <B>Tolerance</B>: double > 0, set tolerance of mapping. The option is ignored if mapping is not active.
 * - <B>FieldName</B>: name of the field to be selected. Default "data".
 *
 * Proper of the class:
 * - <B>Location</B> set unique data location for all fields 1-POINT, 2-CELL, 3-INTERFACE.
 *
 * Fields have to be mandatorily passed through port.
 */
class SelectStringField: public SelectField{
private:
    MPVLocation 									m_loc;		/**< field data reference location */
    std::vector<MimmoPiercedVector<std::string>> 	m_fields;	/**<Input fields to be Select. */
    MimmoPiercedVector<std::string>					m_result;   /**<Result Select fields. */

public:

    SelectStringField(MPVLocation loc = MPVLocation::POINT);
    SelectStringField(const bitpit::Config::Section & rootXMl);
    SelectStringField(const SelectStringField & other);
    virtual ~SelectStringField();

    void buildPorts();
    MimmoPiercedVector<std::string>*     getSelectedField();
    void     setFields(std::vector<MimmoPiercedVector<std::string> *> fields);
    void     addField(MimmoPiercedVector<std::string> *field);

    void clear();

    void     plotOptionalResults();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(SelectStringField &) noexcept;

private:
    bool mSelect();

};



REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_, __SelectFieldS_HPP__)
REGISTER_PORT(M_NAME, MC_SCALAR, MD_STRING, __SelectFieldS_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_, __SelectFieldS_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_, __SelectFieldS_HPP__)
REGISTER_PORT(M_LONGFIELD, MC_SCALAR, MD_MPVECLONG_, __SelectFieldS_HPP__)
REGISTER_PORT(M_STRINGFIELD, MC_SCALAR, MD_MPVECSTRING_, __SelectFieldS_HPP__)
REGISTER_PORT(M_VECSCALARFIELDS, MC_VECTOR, MD_MPVECFLOAT_, __SelectFieldS_HPP__)
REGISTER_PORT(M_VECVECTORFIELDS, MC_VECTOR, MD_MPVECARR3FLOAT_, __SelectFieldS_HPP__)
REGISTER_PORT(M_VECLONGFIELDS, MC_VECTOR, MD_MPVECLONG_, __SelectFieldS_HPP__)
REGISTER_PORT(M_VECSTRINGFIELDS, MC_VECTOR, MD_MPVECSTRING_, __SelectFieldS_HPP__)

REGISTER(BaseManipulation, SelectScalarField, "mimmo.SelectScalarField")
REGISTER(BaseManipulation, SelectVectorField, "mimmo.SelectVectorField")
REGISTER(BaseManipulation, SelectLongField, "mimmo.SelectLongField")
REGISTER(BaseManipulation, SelectStringField, "mimmo.SelectStringField")
};

#endif /* __SelectFieldS_HPP__ */
