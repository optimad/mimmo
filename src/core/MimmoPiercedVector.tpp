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

using namespace std;

namespace mimmo{

/*!
 * Default constructor of MimmoPiercedVector.
 * \param[in] geo pointer to related geometry
 * \param[in] name field's name
 */
template<typename value_t, typename id_t>
MimmoPiercedVector<value_t, id_t>::MimmoPiercedVector(MimmoObject* geo, std::string name){
    m_geometry = geo;
    m_name = name;
}

/*!
 * Default destructor of MimmoPiercedVector.
 */
template<typename value_t, typename id_t>
MimmoPiercedVector<value_t, id_t>::~MimmoPiercedVector(){
    clear();
}

/*! Copy Constructor
 *\param[in] other MimmoPiercedVector object
 */
template<typename value_t, typename id_t>
MimmoPiercedVector<value_t, id_t>::MimmoPiercedVector(const MimmoPiercedVector & other):bitpit::PiercedVector<value_t, id_t>(){
    *this = other;
};

/*! Copy Operator
 * \param[in] other MimmoPiercedVector object
 */
template<typename value_t, typename id_t>
MimmoPiercedVector<value_t, id_t> & MimmoPiercedVector<value_t, id_t>::operator=(const MimmoPiercedVector & other){

    *(static_cast<bitpit::PiercedVector<value_t, id_t> *>(this))  = *(static_cast<const bitpit::PiercedVector<value_t, id_t> *>(&other));
    m_geometry = other.m_geometry;
    m_name = other.m_name;
    return(*this);
};

/*!
 * Clear MimmoPiercedVector.
 */
template<typename value_t, typename id_t>
void
MimmoPiercedVector<value_t, id_t>::clear(){
    m_geometry = NULL;
    m_name = "";
}

/*!
 * Get the linked MimmoObject.
 * return pointer to linked geometry.
 */
template<typename value_t, typename id_t>
MimmoObject*
MimmoPiercedVector<value_t, id_t>::getGeometry() const{
    return m_geometry;
}

/*!
 * Get the name of the field.
 * return name of the data field.
 */
template<typename value_t, typename id_t>
std::string
MimmoPiercedVector<value_t, id_t>::getName() const{
    return m_name;
}

/*!
 * Set the linked MimmoObject.
 * \param[in] geo pointer to linked geometry.
 */
template<typename value_t, typename id_t>
void
MimmoPiercedVector<value_t, id_t>::setGeometry(MimmoObject* geo){
    m_geometry = geo;
}

/*!
 * Set the name of the field.
 * \param[in] name name of the data field.
 */
template<typename value_t, typename id_t>
void
MimmoPiercedVector<value_t, id_t>::setName(std::string name){
    m_name = name;
}


}
