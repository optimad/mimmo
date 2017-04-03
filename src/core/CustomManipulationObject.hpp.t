/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __CUSTOM_HPP__
#define __CUSTOM_HPP__

#include "BaseManipulation.hpp"
#include "<userinclude>.hpp"

namespace mimmo{

/*!
 *  \date           --/---/----
 *  \authors        Developer Name 1
 *  \authors        Developer Name 2
 *
 *  \brief Custom brief description.
 *
 *  Custom description
 *
 */
class CustomClass: public BaseManipulation, public CustomBaseClass{
private:

    type1 m_member1;   /**< Description member 1.*/
    type2 m_member2;   /**< Description member 2.*/


public:
    CustomClass();
    ~CustomClass();

    //copy operators/constructors
    CustomClass(const CustomClass & other);
    CustomClass & operator=(const CustomClass & other);

    //internal methods
    void setMember1(type1 arg1);
    void setMember2(type2 arg2);

    type1 getMember1();
    type2 getMember2();

    //execute deformation methods
    void        execute();

};

}

#endif /* __CUSTOM_HPP__ */
