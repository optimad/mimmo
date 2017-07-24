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

#define __MPV_CONST_REFERENCE__ typename bitpit::PiercedVector<value_t, long int>::const_reference
#define __MPV_REFERENCE__ typename bitpit::PiercedVector<value_t, long int>::reference

namespace mimmo{

/*!
 * \ingroup core
 */


/*!
 * \enum MPVLocation
 * define location of the MimmoPiercedVector data:
 *  - 0-POINT on geometry vertex
 *  - 1-CELL on geometry cells
 *  - 2-INTERFACE on geometry interfaces
 */
enum class MPVLocation{
    POINT=0,
    CELL=1,
    INTERFACE=2
};

/*!
 * \class MimmoPiercedVector
 * \brief MimmoPiercedVector is the basic data container for mimmo library
 *
 * MimmoPiercedVector is the basic container for data attached to a geometry
 * defined as a MimmoObject.
 * It is based on bitpit::PiercedVector container.
 * It supports interface methods to recover the related geometric object.
 * It supports a string name attribute to mark the field as well as a location enum to 
 * understand to which structures of geometry refers the data (POINT-vertices, 
 * CELL-simplicies, INTERFACE-interfaces).
 */
template<typename value_t>
class MimmoPiercedVector{
private:
    MimmoObject*                             m_geometry;            /**<Pointer to geometry. */
    std::string                              m_name;                /**<Name of the field.*/
    MPVLocation                              m_loc;                 /**< MPVLocation enum */
    bitpit::PiercedVector<value_t, long int> m_data;                /**<real data */

    
public:
    MimmoPiercedVector(MimmoObject* geo = NULL, std::string name = "", MPVLocation loc = MPVLocation::POINT);
    virtual ~MimmoPiercedVector();
    //copy constructors and operators
    MimmoPiercedVector(const MimmoPiercedVector & other);
    MimmoPiercedVector & operator=(const MimmoPiercedVector & other);
    
    void            clear();
    void            clearData();
    
    //functional interfaces to inner PiercedVector
    bitpit::PiercedVector<value_t, long int> & data();
    std::size_t            size() const;
    bool                   exists(long int id) const;
    __MPV_CONST_REFERENCE__ operator[](long int id) const;
    __MPV_REFERENCE__       operator[](long int id);
    void                   reserve(std::size_t n);
    void                   resize(std::size_t n);
    
    // get/set methods of the class;
    MimmoObject*            getGeometry() const;
    std::string             getName() const;
    MPVLocation             getDataLocation() const;
    std::vector<long>       getIds(bool ordered=false);
    std::vector<value_t>    getDataAsVector(bool ordered=false);
    
    void                   setGeometry(MimmoObject* geo);
    void                   setName(std::string name);
    void                   setDataLocation(MPVLocation loc);
    void                   setData(bitpit::PiercedVector<value_t, long int >& data);
    void                   setData(std::vector<value_t> &rawdata);

    bool checkDataSizeCoherence();
    bool checkDataIdsCoherence();
    
};

};

#include "MimmoPiercedVector.tpp"

#endif /* __MIMMOPIERCEDVECTOR_HPP__ */



