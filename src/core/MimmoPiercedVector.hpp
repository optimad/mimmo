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
#ifndef __MIMMOPIERCEDVECTOR_HPP__
#define __MIMMOPIERCEDVECTOR_HPP__

#include "piercedVector.hpp"
#include "mimmoTypeDef.hpp"
#include "MimmoObject.hpp"

namespace mimmo{

/*!
 * \ingroup core
 */

/*!
 * \class MimmoPiercedVector
 * \brief MimmoPiercedVector is the basic data container for mimmo library
 *
 * MimmoPiercedVector is the basic container for data attached to a geometry
 * defined as a MimmoObject.
 * It is based on bitpit::PiercedVector container.
 * It supports interface methods to recover the related geometric object.
 * It supports a string name attribute to mark the field.
 */
template<typename value_t, typename id_t = long int>
class MimmoPiercedVector : public bitpit::PiercedVector<value_t, id_t> {

    MimmoObject*                m_geometry;            /**<Pointer to geometry. */
    std::string                 m_name;                /**<Name of the field.*/

public:
    MimmoPiercedVector(MimmoObject* geo = NULL, std::string name = "");
    ~MimmoPiercedVector();

    //copy operators/constructors
    MimmoPiercedVector(const MimmoPiercedVector & other);
    MimmoPiercedVector & operator*(double val);

    void            clear();

    MimmoObject*    getGeometry() const;
    std::string     getName() const;

    void            setGeometry(MimmoObject* geo);
    void            setName(std::string name);

};

};

#include "MimmoPiercedVector.tpp"

#endif /* __MIMMOPIERCEDVECTOR_HPP__ */



