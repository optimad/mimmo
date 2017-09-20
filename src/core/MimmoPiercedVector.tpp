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
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
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
MimmoPiercedVector<value_t>::MimmoPiercedVector(const MimmoPiercedVector<value_t> & other):bitpit::PiercedVector<value_t,long int>(other){
    this->m_geometry = other.m_geometry;
    this->m_loc = other.m_loc;
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
};

/*! 
 * Assignment Operator. 
 * \param[in] other MimmoPiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t> & MimmoPiercedVector<value_t>::operator =(MimmoPiercedVector<value_t> other){
    this->swap(other);
    return *this;
};

/*! 
 * Copy Operator for pierced data only. Values will be stored as is in the inner PiercedVector of the class.
 * \param[in] other PiercedVector object
 */
template<typename value_t>
MimmoPiercedVector<value_t> & MimmoPiercedVector<value_t>::operator =(bitpit::PiercedVector<value_t, long int> other){
    
    this->bitpit::PiercedVector<value_t, long int>::swap(other);
    return *this;
};

/*!
 * Attibutes data of an external object to the current class and vice versa. Containts of both objects will be swapped.
 * Size of the two containers may differ.
 * \param[in] x MimmoPiercedVector to be swapped. Data type of elements of the pierced vector must be of the same type of the current class.
 */
template<typename value_t>
void MimmoPiercedVector<value_t>::swap(MimmoPiercedVector<value_t> & x) noexcept
{
  std::swap(this->m_geometry, x.m_geometry);
  std::swap(this->m_loc, x.m_loc);
  this->bitpit::PiercedVector<value_t, long int>::swap(x);	
}


/*!
 * Clear the whole MimmoPiercedVector.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::clear(){
    m_geometry = NULL;
//     m_name = "";
    m_loc = MPVLocation::UNDEFINED;
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
 * It returns what is stored in the respective member, and does not attempt to
 * recover the reference location from geometry eventually. 
 * \return MPVLocation enum
 */
template<typename value_t>
MPVLocation
MimmoPiercedVector<value_t>::getConstDataLocation() const{
    return m_loc;
}

/*!
 * Get data location w.r.t geometry inner structures.
 * If location is UNDEFINED, attempt to recover reference location from linked
 * geometry employing recoverGeometryReferenceLocation() method.
 * \return MPVLocation enum
 */
template<typename value_t>
MPVLocation
MimmoPiercedVector<value_t>::getDataLocation(){
 /*   if( m_loc == MPVLocation::UNDEFINED){
        return recoverGeometryReferenceLocation();
    }else{
 */       return m_loc;
 //   }
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
 * Set the data Location through integer. If out of MPVLocation enum, 
 * set MPVLocation::UNDEFINED by default.
 * \param[in] loc int 0 to 3 to identify location.
 */
template<typename value_t>
void
MimmoPiercedVector<value_t>::setDataLocation(int loc){
    loc = std::max(0, loc);
    if (loc > 3) loc =0;
    setDataLocation(static_cast<MPVLocation>(loc));
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
    bool check = false;
    switch(m_loc){
        case MPVLocation::CELL:
            check = (this->size()==m_geometry->getPatch()->getCells().size());
            break;
        case MPVLocation::INTERFACE:
            {
                size_t sizeInterfaces = m_geometry->getPatch()->getInterfaces().size();
                check = (this->size()==sizeInterfaces);
                if(sizeInterfaces == 0){
                    (*m_log)<<"Warning: Asked Data Size Coherence in MimmoPiercedVector for INTERFACES, but linked geometry may not have them built."<<std::endl;
                }
            }
            break;
        case MPVLocation::POINT:
            check = (this->size()==m_geometry->getPatch()->getVertices().size());
            break;
        default:
            check=false;
            (*m_log)<<"NO suitable location data found to perform data size coherence check"<<std::endl;
            break;
    }
    return check;
}

/*!
 * Check data coherence with the geometry linked. Return a coherence boolean flag which is 
 * false if:
 *  - UNDEFINED location is set for the current data.
 *  - all internal m_data ids does not match those available in the relative geometry reference structure: vertex, cell or interfaces 
 *  - totally empty vector
 *  - no geometry is linked 
 * \return boolean coherence flag
 */
template<typename value_t>
bool
MimmoPiercedVector<value_t>::checkDataIdsCoherence(){
    if(getGeometry()==NULL) return false;
    auto ids = this->getIds();
    bool check = this->size() > 0;
    switch(m_loc){
        case MPVLocation::CELL:
            {
                auto vcell = m_geometry->getPatch()->getCells();
                for(auto el : ids){
                    check =check && vcell.exists(el);
                }
            }
            break;
        case MPVLocation::INTERFACE:
            {
                size_t sizeInterfaces = m_geometry->getPatch()->getInterfaces().size();
                if(sizeInterfaces == 0){
                    (*m_log)<<"Warning: Asked Data Ids Coherence in MimmoPiercedVector for INTERFACES, but linked geometry may not have them built."<<std::endl;
                }
                auto vint = m_geometry->getPatch()->getInterfaces();
                for(auto el : ids){
                    check =check && vint.exists(el);
                }
            }
        break;
        case MPVLocation::POINT:
            {
                auto vvert = m_geometry->getPatch()->getVertices();
                for(auto el : ids){
                    check =check && vvert.exists(el);
                }
            }
        break;
        default:
            check = false;
            (*m_log)<<"NO suitable location data found to perform ids coherence check"<<std::endl;
            break;
    }
    return check;
}

///*!
// * Recover the most probable location of your data with reference to the MimmoObject structures
// * available, i.e. vertices, cells or interfaces. If size does not match any of those structures, 
// * return an MPVLocation::UNDEFINED value.
// * 
// * \return MPVLocation enum value
// */
//template<typename value_t>
// MPVLocation
// MimmoPiercedVector<value_t>::recoverGeometryReferenceLocation(){
//     
//     if(getGeometry()==NULL) {
//         m_loc = MPVLocation::UNDEFINED;
//         return MPVLocation::UNDEFINED;
//     }    
//     std::size_t datasize = this->size();
// 
//     if(datasize == getGeometry()->getPatch()->getVertexCount()){
//         m_loc = MPVLocation::POINT;
//         return MPVLocation::POINT;
//     }
//     if(datasize == getGeometry()->getPatch()->getCellCount()){
//         m_loc = MPVLocation::CELL;
//         return MPVLocation::CELL;
//     }
//     if(datasize == getGeometry()->getPatch()->getInterfaceCount()){
//         m_loc = MPVLocation::INTERFACE;
//         return MPVLocation::INTERFACE;
//     }
//     m_loc = MPVLocation::UNDEFINED;
//     return MPVLocation::UNDEFINED;
// }

/*!
 * Check if a random integer number is a valid MPVLocation for the current class.
 * \return true if valid.
 */
template<typename value_t>
bool
MimmoPiercedVector<value_t>::intIsValidLocation(int &value){
    return !(value<0 && value>3) ;
}

}
