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
 * \param[in] name field name
 */
template<typename value_t>
MimmoPiercedVector<value_t>::MimmoPiercedVector(MimmoObject* geo, MPVLocation loc):bitpit::PiercedVector<value_t,long int>(){
    m_geometry = geo;
    m_loc = loc;
//     m_name = name;
}

/*!
 * Default destructor of MimmoPiercedVector.
 */
template<typename value_t>
MimmoPiercedVector<value_t>::~MimmoPiercedVector(){
    clear();
}

/*! 
 * Copy Constructor
 *\param[in] other MimmoPiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t>::MimmoPiercedVector(const MimmoPiercedVector<value_t> & other):bitpit::PiercedVector<value_t,long int>(){
    *this = other;
};

/*! 
 * Copy Operator. Values of the Pierced will be copied without holes.
 * \param[in] other MimmoPiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t> & MimmoPiercedVector<value_t>::operator =(const MimmoPiercedVector<value_t> & other){
    this->m_geometry = other.m_geometry;
    this->m_loc = other.m_loc;
//     this->m_name = other.m_name;
    //copy entirely the data in internal pierced vector->
    auto itE = other.end();
    for(auto it  = other.begin(); it!=itE; ++it){
        this->insert(it.getId(),*it);
    }
    return (*this);
};

/*! 
 * Copy Operator for pierced data only. Only Values of the Pierced will be copied without holes.
 * \param[in] other PiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t> & MimmoPiercedVector<value_t>::operator =(const bitpit::PiercedVector<value_t, long int> & other){
    auto itE = other.end();
    for(auto it  = other.begin(); it!=itE; ++it){
        this->insert(it.getId(),*it);
    }
    return (*this);
};


/*!
 * Clear the whole MimmoPiercedVector.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::clear(){
    m_geometry = NULL;
//     m_name = "";
    m_loc = MPVLocation::POINT;
    bitpit::PiercedVector<value_t, long int>::clear();
}

/*!
 * Get the linked MimmoObject.
 * \return pointer to linked geometry.
 */
template<typename value_t>
MimmoObject*
MimmoPiercedVector<value_t>::getGeometry() const{
    return m_geometry;
}

// /*!
//  * Get the name of the field.
//  * \return name of the data field.
//  */
// template<typename value_t>
// std::string
// MimmoPiercedVector<value_t>::getName() const{
//     return m_name;
// }

/*!
 * Get data location w.r.t geometry inner structures.
 * \return MPVLocation enum
 */
template<typename value_t>
MPVLocation
MimmoPiercedVector<value_t>::getDataLocation() const{
    return m_loc;
}

/*!
 * Return data contained in inner pierced vector
 * \param[in] ordered, if true data will be returned in ids ascending order, otherwise they will be returned as 
 * you get iterating the internal m_data PiercedVector from the beginning.
 * \return list of data 
 * 
 */
template<typename value_t>
std::vector<value_t>
MimmoPiercedVector<value_t>::getDataAsVector(bool ordered){
    
     auto ids = this->getIds(ordered);
     std::vector<value_t> result(ids.size());
     int counter= 0;
     for(const auto val: ids){
         result[counter] = (*this)[val];
         ++counter;
     }
     
     return result;
}

/*!
 * Set the linked MimmoObject.
 * \param[in] geo pointer to linked geometry.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setGeometry(MimmoObject* geo){
    m_geometry = geo;
}

// /*!
//  * Set the name of the field.
//  * \param[in] name name of the data field.
//  */
// template<typename value_t>
// void
// MimmoPiercedVector<value_t>::setName(std::string name){
//     m_name = name;
// }

/*!
 * Set the data Location
 * \param[in] loc MPVLocation enum
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setDataLocation(MPVLocation loc){
    m_loc = loc;
}

/*!
 * Set the data of the inner PiercedVector from a row data compound.
 * Ids will be automatically assigned.
 * \param[in] data vector to copy from
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setData(std::vector<value_t>& data){
    bitpit::PiercedVector<value_t, long int>::clear();
    long int  id = 0; 
    for(const auto val: data){
        this->insert( id, val);
        id++;
    }
}

/*!
 * Check data coherence with the geometry linked. Return a coherence boolean flag which is 
 * false if:
 *  - UNDEFINED location is set to the current data
 *  - internal m_data does not match the size of the relative geometry reference structure: vertex, cell or interfaces 
 *  - no geometry is linked 
 * \return boolean coherence flag
 */
template<typename value_t>
bool
MimmoPiercedVector<value_t>::checkDataSizeCoherence(){
    if(getGeometry()==NULL) return false;
    bool check = true;
    switch(m_loc){
        case 2:
            check = (this->size()==m_geometry->getPatch()->getCells().size());
            break;
        case 3:
            check = (this->size()==m_geometry->getPatch()->getInterfaces().size());
            break;
        case 1:
            check = (this->size()==m_geometry->getPatch()->getVertices().size());
            break;
        default:
            check=false;
            break;
    }
    return check;
}

/*!
 * Check data coherence with the geometry linked. Return a coherence boolean flag which is 
 * false if:
 *  - UNDEFINED location is set for the current data.
 *  - all internal m_data ids does not match those available in the relative geometry reference structure: vertex, cell or interfaces 
 *  - empty inner data structure
 *  - no geometry is linked 
 * \return boolean coherence flag
 */
template<typename value_t>
bool
MimmoPiercedVector<value_t>::checkDataIdsCoherence(){
    if(getGeometry()==NULL) return false;
    bool check = this->size() != std::size_t(0);
    switch(m_loc){
        case 2:
            for(auto &el : m_geometry->getPatch()->getCells()){
                check =check && this->exists(el.getId());
            }    
            break;
        case 3:
            for(auto &el : m_geometry->getPatch()->getInterfaces()){
                check =check && this->exists(el.getId());
            }    
            break;
        case 1:
            for(auto &el : m_geometry->getPatch()->getVertices()){
                check =check && this->exists(el.getId());
            }    
            break;
        default:
            check = false;
            break;
    }
    return check;
}

/*!
 * Recover the most probable location of your data with reference to the MimmoObject structures
 * available, i.e. vertices, cells or interfaces. If size does not match any of those structures, 
 * return an MPVLocation::UNDEFINED value.
 * 
 * \return MPVLocation enum value
 */
template<typename value_t>
MPVLocation
MimmoPiercedVector<value_t>::recoverGeometryReferenceLocation(){
    
    if(getGeometry()==NULL) return MPVLocation::UNDEFINED;
    std::size_t datasize = this->size();

    if(datasize == getGeometry()->getPatch()->getVertexCount()){
        return MPVLocation::POINT;
    }
    if(datasize == getGeometry()->getPatch()->getCellCount()){
        return MPVLocation::CELL;
    }
    if(datasize == getGeometry()->getPatch()->getInterfaceCount()){
        return MPVLocation::INTERFACE;
    }
    return MPVLocation::UNDEFINED;
}

}
