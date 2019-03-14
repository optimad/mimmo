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


namespace mimmo{

/*!
 * Default constructor of MimmoPiercedVector.
 * \param[in] geo pointer to related geometry
 * \param[in] loc reference location of field. see MPVLocation.
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t>::MimmoPiercedVector(MimmoObject* geo, MPVLocation loc):bitpit::PiercedVector<mpv_t,long int>(){
    m_geometry = geo;
    m_loc = loc;
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
}

/*!
 * Default destructor of MimmoPiercedVector.
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t>::~MimmoPiercedVector(){
    clear();
}

/*! 
 * Copy Constructor
 *\param[in] other MimmoPiercedVector object
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t>::MimmoPiercedVector(const MimmoPiercedVector<mpv_t> & other):bitpit::PiercedVector<mpv_t,long int>(other){
    this->m_geometry = other.m_geometry;
    this->m_loc = other.m_loc;
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
};

/*! 
 * Assignment Operator. 
 * \param[in] other MimmoPiercedVector object
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> & MimmoPiercedVector<mpv_t>::operator =(MimmoPiercedVector<mpv_t> other){
    this->swap(other);
    return *this;
};

/*! 
 * Copy Operator for pierced data only. Values will be stored as is in the inner PiercedVector of the class.
 * \param[in] other PiercedVector object
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> & MimmoPiercedVector<mpv_t>::operator =(bitpit::PiercedVector<mpv_t, long int> other){
    
    this->bitpit::PiercedVector<mpv_t, long int>::swap(other);
    return *this;
};

/*!
 * Attibutes data of an external object to the current class and vice versa. Containts of both objects will be swapped.
 * Size of the two containers may differ.
 * \param[in] x MimmoPiercedVector to be swapped. Data type of elements of the pierced vector must be of the same type of the current class.
 */
template<typename mpv_t>
void MimmoPiercedVector<mpv_t>::swap(MimmoPiercedVector<mpv_t> & x) noexcept
{
  std::swap(this->m_geometry, x.m_geometry);
  std::swap(this->m_loc, x.m_loc);
  this->bitpit::PiercedVector<mpv_t, long int>::swap(x);	
}


/*!
 * Clear the whole MimmoPiercedVector.
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::clear(){
    m_geometry = NULL;
//     m_name = "";
    m_loc = MPVLocation::UNDEFINED;
    bitpit::PiercedVector<mpv_t, long int>::clear();
}

/*!
 * Get the linked MimmoObject.
 * \return pointer to linked geometry.
 */
template<typename mpv_t>
MimmoObject*
MimmoPiercedVector<mpv_t>::getGeometry() const{
    return m_geometry;
}

/*!
 * Get data location w.r.t geometry inner structures.
 * It returns what is stored in the respective member, and does not attempt to
 * recover the reference location from geometry eventually. 
 * \return MPVLocation enum
 */
template<typename mpv_t>
MPVLocation
MimmoPiercedVector<mpv_t>::getConstDataLocation() const{
    return m_loc;
}

/*!
 * Get data location w.r.t geometry inner structures.
 * If location is UNDEFINED, attempt to recover reference location from linked
 * geometry employing recoverGeometryReferenceLocation() method.
 * \return MPVLocation enum
 */
template<typename mpv_t>
MPVLocation
MimmoPiercedVector<mpv_t>::getDataLocation(){
        return m_loc;
}

/*!
 * Return data contained in inner pierced vector. Sequence follows that of reference location in
 * geometry(vertices, cells or interfaces). If no geometry is provided, return empty result.
 * \param[in] ordered if true data will be returned in ids ascending order, otherwise they will be returned as 
 * you get iterating the internal location reference geometry PiercedVector from the beginning.
 * \return list of data 
 * 
 */
template<typename mpv_t>
std::vector<mpv_t>
MimmoPiercedVector<mpv_t>::getDataAsVector(bool ordered){
     if(getGeometry() == NULL) return std::vector<mpv_t>(0);
     livector1D ids = getGeometryIds(ordered);
     std::vector<mpv_t> result(this->size());
     int counter= 0;
     for(const auto val: ids){
         if(this->exists(val)){
            result[counter] = (*this)[val];
            ++counter;
         }   
     }
     
     return result;
}

/*!
 * Return only raw data contained in inner pierced vector. Sequence follows the internal pierced vector id-ing, 
 * without any reference to geometry structure ordering.
 * \param[in] ordered if true data will be returned in ids ascending order, otherwise they will be returned as 
 * you get iterating the class object itself from the beginning.
 * \return list of data 
 * 
 */
template<typename mpv_t>
std::vector<mpv_t>
MimmoPiercedVector<mpv_t>::getRawDataAsVector(bool ordered){
    livector1D ids = this->getIds(ordered);
    std::vector<mpv_t> result(this->size());
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
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::setGeometry(MimmoObject* geo){
    m_geometry = geo;
}

/*!
 * Set the data Location
 * \param[in] loc MPVLocation enum
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::setDataLocation(MPVLocation loc){
    m_loc = loc;
}

/*!
 * Set the data Location through integer. If out of MPVLocation enum, 
 * set MPVLocation::UNDEFINED by default.
 * \param[in] loc int 0 to 3 to identify location.
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::setDataLocation(int loc){
    loc = std::max(0, loc);
    if (loc > 3) loc =0;
    setDataLocation(static_cast<MPVLocation>(loc));
}

/*!
 * Set the data of the inner PiercedVector from a row data compound.
 * Ids will be automatically assigned.
 * \param[in] data vector to copy from
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::setData(std::vector<mpv_t>& data){
    bitpit::PiercedVector<mpv_t, long int>::clear();
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
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::checkDataSizeCoherence(){
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
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::checkDataIdsCoherence(){
    if(getGeometry()==NULL) return false;
    auto ids = this->getIds();
    bool check = !this->isEmpty();
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

/*!
 * Check if a random integer number is a valid MPVLocation for the current class.
 * \return true if valid.
 */
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::intIsValidLocation(int &value){
    return !(value<0 && value>3) ;
}

/*!
 * Return linked Geoemtry Ids of MPVLocation assigned: POINT-geometry vertices, CELL-geometry cells
 * and INTERFACE-geometry interfaces if any.
 * \param[in] ordered to force ascending ordering of ids.
 * \return list of ids
 */
template<typename mpv_t>
livector1D
MimmoPiercedVector<mpv_t>::getGeometryIds(bool ordered){
    if(getGeometry()==NULL) return livector1D(0);
    switch(m_loc){
        case MPVLocation::POINT:
            return getGeometry()->getVertices().getIds(ordered);
            break;
        case MPVLocation::CELL:
            return getGeometry()->getCells().getIds(ordered);
            break;
        case MPVLocation::INTERFACE:
            {
                size_t sizeInterfaces = m_geometry->getPatch()->getInterfaces().size();
                if(sizeInterfaces == 0){
                    (*m_log)<<"Warning: Asked list of geometry Ids in MimmoPiercedVector for INTERFACES, but linked geometry may not have them built."<<std::endl;
                }
                return getGeometry()->getInterfaces().getIds(ordered);
            }    
            break;
        default:
            return livector1D(0);
        break;
    }
}

/*!
 * \return true if current pierced container has no element in it.
 */
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::isEmpty(){
    return this->size() == size_t(0);
}

/*!
 * Check if container current data are coherent with the geometry linked. If it is and
 * current data size does not match the size of the reference geometry container,
 * attempt to complete all values in the missing ids of reference location
 * geometry structure with a User-assigned reference value. 
 * \param[in] defValue User-assigned reference value
 * \return true if the vector is coherent and full values aligned with geoemtry reference structure. 
 */
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::completeMissingData(const mpv_t & defValue){
    
    if(!this->checkDataIdsCoherence()) return false;
    if(!this->checkDataSizeCoherence()){
        
        livector1D ids = this->getGeometryIds();
        for(auto id: ids){
            if(!this->exists(id)) this->insert(id, defValue);
        }
    }    
    return true;
}

/*!
 * Initialize your container with reference data.
 * \param[in] geo target geometry
 * \param[in] loc  data location
 * \param[in] data reference data
 * 
 * if a valid geometry and coherent location are specified, create 
 * a container on all elements of specified location with constant reference data attached.
 * Any pre-existent data will be destroyed.
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::initialize(MimmoObject * geo, MPVLocation loc, const mpv_t & data){
    
    switch(loc){
        case MPVLocation::POINT :
            if(geo->getNVertex() == 0){
                (*m_log)<<"MimmoPiercedVector warning: initialization failed"<<std::endl;
                return;
            }else{
                this->clear();
                this->reserve(geo->getNVertex());
                m_geometry = geo;
                m_loc = loc;
                for(const auto & vertex: geo->getVertices()){
                    this->insert(vertex.getId(), data);
                }
            }
            break;
        case MPVLocation::CELL :
            if(geo->getNCells() == 0){
                (*m_log)<<"MimmoPiercedVector warning: initialization failed"<<std::endl;
                return;
            }else{
                this->clear();
                this->reserve(geo->getNCells());
                m_geometry = geo;
                m_loc = loc;
                for(const auto & cell: geo->getCells()){
                    this->insert(cell.getId(), data);
                }
            }
            break;
        case MPVLocation::INTERFACE :
            if(!geo->areInterfacesBuilt()){
                (*m_log)<<"MimmoPiercedVector warning: initialization failed"<<std::endl;
                return;
            }else{
                this->clear();
                this->reserve(geo->getPatch()->getInterfaceCount());
                m_geometry = geo;
                m_loc = loc;
                for(const auto & interf: geo->getInterfaces()){
                    this->insert(interf.getId(), data);
                }
            }
            break;
        default:
            //do nothing
            break;
    }
}

/*!
 * Point data to Cell data interpolation. Average of point data is set on cell center.
 * \param[in] pointData MimmoPiercedVector object located on MPVLocation::POINT
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::pointDataToCellData(){
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> cellData(geo, MPVLocation::CELL);
	for (bitpit::Cell & cell : geo->getCells()){
		long idcell = cell.getId();
		mpv_t data;
		bool init = false;
		for (long idvertex : cell.getVertexIds()){
			if (!init){
				data = this->at(idvertex);
			}
			else{
				data = data + this->at(idvertex);
			}
		}
		data = data / cell.getVertexCount();
		cellData.insert(idcell, data);
	}
	return cellData;
};

/*!
 * Cell data to Point data interpolation. Average of cell center data is set on point.
 * \param[in] cellData MimmoPiercedVector object located on MPVLocation::CELL
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::cellDataToPointData(){
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> pointData(geo, MPVLocation::POINT);
	MimmoPiercedVector<int> countCells(geo, MPVLocation::POINT);
	for (bitpit::Cell & cell : geo->getCells()){
		long idcell = cell.getId();
		for (long idvertex : cell.getVertexIds()){
			if (!pointData.exists(idvertex)){
				pointData.insert(idvertex, this->at(idvertex));
				countCells.insert(idvertex, 1);
			}
			else{
				pointData[idvertex] = pointData[idvertex] + this->at(idcell);
				countCells[idvertex] = countCells[idvertex] + 1;
			}
		}
	}
	for (bitpit::Vertex & vertex : geo->getVertices()){
		long idvertex = vertex.getId();
		pointData[idvertex] = pointData[idvertex] / countCells[idvertex];
	}
	return pointData;
};

}
