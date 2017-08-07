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
 *    After the execution of an object Apply, the original geometry will be modified.
 *
 * \n
 * Ports available in Apply Class :
 *
 *    =========================================================
 * 
     | Port Input | | | |                                                              
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> | 
     | 11    | M_GDISPLS | setInput          | (MPVECARR3,FLOAT) |
     | 99    | M_GEOM    | setGeometry       | (SCALAR,MIMMO_) |
 
     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |    
     | M_GEOM   | getGeometry       | (SCALAR,MIMMO_) |
 
 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.Apply</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 *
 * Geometry and Input have to be mandatorily passed through port.
 *
 */
class Apply: public BaseManipulation{
public:

    dmpvecarr3E    m_input; /**< storing vector fields of floats */

    Apply();
    Apply(const bitpit::Config::Section & rootXML);

    ~Apply();

    Apply(const Apply & other);

    void buildPorts();

    void setInput(dmpvecarr3E input);

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
protected:
    void swap(Apply & x) noexcept;
    
};

REGISTER(BaseManipulation, Apply, "mimmo.Apply")

};

#endif /* __APPLYDEFORMATION_HPP__ */
