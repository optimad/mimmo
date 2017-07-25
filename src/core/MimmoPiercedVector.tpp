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
MimmoPiercedVector<value_t>::MimmoPiercedVector(MimmoObject* geo, MPVLocation loc){
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
MimmoPiercedVector<value_t>::MimmoPiercedVector(const MimmoPiercedVector & other){
    *this = other;
};

/*! 
 * Copy Operator. Values of the Pierced will be copied without holes.
 * \param[in] other MimmoPiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t> & MimmoPiercedVector<value_t>::operator =(const MimmoPiercedVector & other){
    this->m_geometry = other.m_geometry;
    this->m_loc = other.m_loc;
//     this->m_name = other.m_name;
    //copy entirely the data in internal pierced vector->
    auto itE = other.m_data.end();
    for(auto it  = other.m_data.begin(); it!=itE; ++it){
        this->m_data.insert(it.getId(),*it);
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
    clearData();
}

/*!
 * Clear internal Pierced Vector data only
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::clearData(){
    m_data.clear(true);
}

/*!
 * \return reference to internal bitpit::PiercedVector<value_t, long int> data
 */
template<typename value_t>
bitpit::PiercedVector<value_t, long int>&
MimmoPiercedVector<value_t>::data() {
    return m_data;
}

/*!
 * \return const reference to internal bitpit::PiercedVector<value_t, long int> data
 */
template<typename value_t>
const bitpit::PiercedVector<value_t, long int>&
MimmoPiercedVector<value_t>::data() const  {
    return m_data;
}

/*!
 * \return size of the internal PiercedVector data, which are not necessarily equal 
 * to its current capacity.
 */
template<typename value_t>
std::size_t
MimmoPiercedVector<value_t>::size() const{
    return m_data.size();
}

/*!
 * \return true if an element flagged by id is present in the data stock.
 * \param[in] id, long int marker
 */
template<typename value_t>
bool
MimmoPiercedVector<value_t>::exists(long int id) const{
    return m_data.exists(id);
}

/*!
 * Const access operator
 * \return const reference to the element marked as id in internal pv data.
 * \param[in] id, long int marker
 */
template<typename value_t>
__MPV_CONST_REFERENCE__
MimmoPiercedVector<value_t>::operator[](long int id) const{
    return m_data[id];
}

/*!
 * Access operator
 * \return reference to the element marked as id in internal pv data.
 * \param[in] id, long int marker
 */
template<typename value_t>
__MPV_REFERENCE__
MimmoPiercedVector<value_t>::operator[](long int id){
    return m_data[id];
}

/*!
 * Reserve stock of memory for internal data. This method do not size your internal data.
 * Use resize method instead.
 * \param[in] n, number of value_t elements to reserve.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::reserve(std::size_t n){
    return m_data.reserve(n);
}

/*!
 * Resize for internal data pierced vector.
 * \param[in] n, number of value_t elements to resize.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::resize(std::size_t n){
    return m_data.resize(n);
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
 * Get Id-etiquettes attached to each value store in inner data. 
 * \param[in] ordered, if true ids will be returned in ascending order, otherwise they will be returned as 
 * you get iterating the internal m_data PiercedVector from the beginning.
 * \return list of id
 * 
 */
template<typename value_t>
std::vector<long int>
MimmoPiercedVector<value_t>::getIds(bool ordered)const {
    if(ordered) return m_data.getIds(ordered);
    
    std::vector<long int> result(m_data.size());
    auto itE= m_data.end();
    int counter= 0;
    for(auto it = m_data.begin(); it != itE; ++it){
        result[counter] = it.getId();
        ++counter;
    }
    return result;
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
    
     auto ids = getIds(ordered);
     std::vector<value_t> result(ids.size());
     int counter= 0;
     for(const auto val: ids){
         result[counter] = m_data[val];
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
 * Set the data of the inner PiercedVector. Data will be copied and stored internally.
 * All precedent inner data will be cleared.
 * \param[in] data Pierced Vector to copy data form
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setData(bitpit::PiercedVector<value_t, long int>& data){
    m_data.clear();
    auto itE= data.end();
    for(auto it = data.begin(); it !=itE; ++it){
        m_data.insert(it.getId(),*it);
    }
}

/*!
 * Set the data of the inner PiercedVector from a row data compound.
 * Ids will be automatically assigned.
 * \param[in] data vector to copy from
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setData(std::vector<value_t>& data){
    m_data.clear();
    long int  id = 0; 
    for(const auto val: data){
        m_data.insert( id, val);
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
            check = (m_data.size()==m_geometry->getPatch()->getCells().size());
            break;
        case 3:
            check = (m_data.size()==m_geometry->getPatch()->getInterfaces().size());
            break;
        case 1:
            check = (m_data.size()==m_geometry->getPatch()->getVertices().size());
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
    bool check = m_data.size() != std::size_t(0);
    switch(m_loc){
        case 2:
            for(auto &el : m_geometry->getPatch()->getCells()){
                check =check && m_data.exists(el.getId());
            }    
            break;
        case 3:
            for(auto &el : m_geometry->getPatch()->getInterfaces()){
                check =check && m_data.exists(el.getId());
            }    
            break;
        case 1:
            for(auto &el : m_geometry->getPatch()->getVertices()){
                check =check && m_data.exists(el.getId());
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
    auto datasize = m_data.size();

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
