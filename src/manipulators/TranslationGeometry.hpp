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
#ifndef __TRANSLATIONGEOMETRY_HPP__
#define __TRANSLATIONGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *    \class TranslationGeometry
 *    \ingroup manipulators
 *    \brief TranslationGeometry is the class that applies a translation to a given geometry patch.
 *
 *    The used parameters are the translation value and the direction of the translation axis.
 *
 * \n
 * Ports available in TranslationGeometry Class :
 *
 *    =========================================================
 
     |Port Input | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 21    | M_AXIS   | setDirection      | (ARRAY3, FLOAT)       |
     | 30    | M_VALUED | setTranslation    | (SCALAR, FLOAT)       |
     | 12    | M_FILTER | setFilter         | (MPVECTOR, FLOAT)       |
     | 99    | M_GEOM   | setGeometry       | (SCALAR, MIMMO_)      |
 
     |Port Output | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | 11    | M_GDISPLS | getDisplacements  | (MPVECARR3, FLOAT)      |
     | 99    | M_GEOM   | getGeometry       | (SCALAR,MIMMO_) |
 
 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.TranslationGeometry</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 *
 * Proper of the class:
 * - <B>Direction</B>: axis direction coordinates;
 * - <B>Translation</B>: translation value in length unity.

 *
 * Geometry has to be mandatorily passed through port.
 *
 */
class TranslationGeometry: public BaseManipulation{
private:
    //members
    darray3E    m_direction;    /**<Components of the translation axis.*/
    double      m_alpha;        /**<Angle of translation in radiant. */
    dmpvector1D   m_filter;      /**<Filter field for displacements modulation. */
    dmpvecarr3E   m_displ;       /**<Resulting displacements of geometry vertex.*/

public:
    TranslationGeometry(darray3E direction = { {0, 0, 0} });
    TranslationGeometry(const bitpit::Config::Section & rootXML);
    ~TranslationGeometry();

    TranslationGeometry(const TranslationGeometry & other);

    void        buildPorts();

    void        setDirection(darray3E direction);
    void        setTranslation(double alpha);
    void        setFilter(dmpvector1D filter);

    dmpvecarr3E   getDisplacements();

    void         execute();
    void         apply();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
};

REGISTER(BaseManipulation, TranslationGeometry, "mimmo.TranslationGeometry")

};

#endif /* __TRANSLATIONGEOMETRY_HPP__ */
