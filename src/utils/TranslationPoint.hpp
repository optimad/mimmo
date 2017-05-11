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
#ifndef __TRANSLATIONPOINT_HPP__
#define __TRANSLATIONPOINT_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class TranslationPoint
 * \ingroup utils
 * \brief TranslationPoint is the class that applies the a translation to a point.
 *
 * Here the point is supposed to be an origin of a reference system.
 * The point is translated over a direction and for a quantity set by the user or
 * an external input.
 * Result of the translation are saved in result of base class and
 * in the modified member m_origin.
 * 
 * \n
 * Ports available in GenericInput Class :
 *
 *    =========================================================

     |Port Input | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 20    | M_POINT  | m_origin          | (ARRAY3, FLOAT)        |
     | 21    | M_AXIS   | m_direction       | (ARRAY3, FLOAT)        |
     | 30    | M_VALUED | m_alpha           | (SCALAR, FLOAT)        |

     |Port Output | | | |
     |-|-|-|-|
     |<B>PortID</B> | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | 20    | M_POINT   | getOrigin         | (ARRAY3, FLOAT)       |
     | 22    | M_AXIS    | getDirection      | (ARRAY3, FLOAT)       |
     | 30    | M_VALUED  | getTranslation    | (SCALAR, FLOAT)       |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.TranslationPoint</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Origin</B>: point tha need to be translated;
 * - <B>Direction</B>: translation direction;
 * - <B>Translation</B>: entity of translation.
 *
 */
class TranslationPoint: public BaseManipulation{
private:
    //members
    darray3E    m_direction;    /**<Components of the translation axis.*/
    darray3E    m_origin;       /**<Origin of box to be deformed.*/
    double      m_alpha;        /**<Value of translation.*/

public:
    TranslationPoint(darray3E direction = { {0, 0, 0} });
    TranslationPoint(const bitpit::Config::Section & rootXML);
    ~TranslationPoint();

    TranslationPoint(const TranslationPoint & other);
    TranslationPoint & operator=(const TranslationPoint & other);

    void        buildPorts();

    darray3E     getDirection();
    darray3E     getOrigin();
    double         getTranslation();
    void         setDirection(darray3E direction);
    void         setOrigin(darray3E origin);
    void         setTranslation(double alpha);

    void     execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

};

REGISTER(BaseManipulation, TranslationPoint, "mimmo.TranslationPoint")
};

#endif /* __TRANSLATIONPOINT_HPP__ */
