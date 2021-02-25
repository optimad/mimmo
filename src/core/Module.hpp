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
#ifndef __MODULE_HPP__
#define __MODULE_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class Module
 * \ingroup core
 * \brief Module is an executable block class capable of
 *         computing the magnitude field of a vector field.
 *
 * Module takes as input a vector field with any data location (POINT, CELLS or INTERFACES), computes its magnitude
 * and it gives the scalar field as result.
 *
 * Ports available in Module Class :
 *
 *    =========================================================
 *
     |                 Port Input     ||                                                     |
     |------------|------------------------------------|-----------------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_VECTORFIELD| setField           | (MC_SCALAR, MD_MPVECARR3FLOAT_)|


     |            Port Output         ||             |
     |------------|------------------------------------|-----------------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_SCALARFIELD  | getResult       | (MC_SCALAR, MD_MPVECFLOAT_)       |

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.Module</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 */
class Module: public BaseManipulation{
protected:
    dmpvecarr3E m_field;         /**<Input field to be extract. */
    dmpvector1D m_result;               /**<Result extract field. */

public:
    Module();
    Module(const bitpit::Config::Section & rootXML);
    virtual ~Module();

    Module(const Module & other);
    Module& operator=(Module other);

    void buildPorts();

    void     	 setField(dmpvecarr3E*field);
    dmpvector1D* getResult();

    void    plotOptionalResults();
    void	execute();
    void    clear();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(Module & x) noexcept;

};

REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__MODULE_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__MODULE_HPP__)

REGISTER(BaseManipulation, Module, "mimmo.Module")
};

#endif /* __MODULE_HPP__ */
