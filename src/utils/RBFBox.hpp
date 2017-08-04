/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-201/ OPTIMAD engineering Srl
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
#ifndef __RBFBox_HPP__
#define __RBFBox_HPP__

#include "BaseManipulation.hpp"



namespace mimmo{

/*!
 *    \class RBFBox
 *    \ingroup utils
 *    \brief Radial Basis Functions Bounding Box calculator.
 *
 *    Builds the axes aligned bounding box of a set of RBF points.
 *
 * \n
 * Ports available in RBFBox Class :
 *
 *    =========================================================
 * 
     | Port Input |  | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_COORDS    | setNode                               | (MC_VECARR3, MD_FLOAT)      |
     | M_VALUED    | setSupportRadius                      | (MC_SCALAR, MD_FLOAT)       |

     
     
     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>        |
     | M_POINT     | getOrigin         | (MC_ARRAY3, MD_FLOAT)     |
     | M_AXES      | getAxes           | (MC_ARR3ARR3, MD_FLOAT)   |
     | M_SPAN      | getSpan           | (MC_ARRAY3, MD_FLOAT)     |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.RBFBox</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>SupportRadius</B>: Influence Radius value for RBF cloud in input;
 *
 * Nodes list has to be mandatorily passed through port.
 *
 */
class RBFBox: public BaseManipulation {

protected:
    darray3E    m_origin;        /**< Origin of the BB.*/
    darray3E    m_span;         /**< Span of the BB. */
    dmatrix33E    m_axes;            /**<reference system of the bbox (axes aligned AABB) */
    dvecarr3E   m_nodes;        /**<Radial Basis Functions control points.*/
    double      m_suppR;        /**<Support radius value of RBFs.*/


public:
    RBFBox();
    RBFBox(const bitpit::Config::Section & rootXML);
    virtual ~RBFBox();

    //copy operators/constructors
    RBFBox(const RBFBox & other);

    void buildPorts();

    //clean structure;
    void         clearRBFBox();

    //internal methods
    darray3E    getOrigin();
    darray3E    getSpan();
    dmatrix33E    getAxes();
    void        setNode(dvecarr3E);
    void        setSupportRadius(double suppR_);

    //plotting wrappers
    void        plot(std::string directory, std::string filename, int counter, bool binary);

    //building method
    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");


protected:
    virtual void plotOptionalResults();
    void swap(RBFBox & x) noexcept;
};

REGISTER(BaseManipulation,RBFBox, "mimmo.RBFBox")
};

#endif /* __RBFBOX_HPP__ */
