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
MimmoPiercedVector<mpv_t>::~MimmoPiercedVector(){}

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
 * \return list of data. Data on ghost cells are retained.
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
 * Return data contained in inner pierced vector located on internal cells. Sequence follows that of reference location in
 * geometry(vertices, cells, interfaces). If no geometry is provided, return empty result.
 * \param[in] ordered if true data will be returned in ids ascending order, otherwise they will be returned as
 * you get iterating the internal location reference geometry PiercedVector from the beginning.
 * \return list of data. Only data on structures of internal cells are retained.
 *\param[in]	squeeze		if true the result container is squeezed once full
 */
template<typename mpv_t>
std::vector<mpv_t>
MimmoPiercedVector<mpv_t>::getInternalDataAsVector(bool ordered, bool squeeze){
	if(getGeometry() == NULL) return std::vector<mpv_t>(0);
	std::vector<mpv_t> result;
	livector1D ids;
	std::unordered_set<long> ids_;
	switch (getDataLocation())
	{
	case MPVLocation::CELL:
		result.reserve(getGeometry()->getPatch()->getInternalCount());
		ids = getGeometryIds(ordered);
		for(const auto val: ids){
			if(this->exists(val) && getGeometry()->getPatch()->getCell(val).isInterior()){
				result.push_back((*this)[val]);
			}
		}
		break;
	case MPVLocation::POINT:
		result.reserve(getGeometry()->getPatch()->getVertexCount());
		for (auto & cell : getGeometry()->getCells()){
			if (cell.isInterior()){
				for (long id : cell.getVertexIds())
					ids_.insert(id);
			}
		}
		ids = getGeometryIds(ordered);
		for(const auto val: ids){
			if(this->exists(val) && ids_.count(val)){
				result.push_back((*this)[val]);
			}
		}
		break;
	case MPVLocation::INTERFACE:
		result.reserve(getGeometry()->getPatch()->getInterfaceCount());
		//Interfaces ghost/ghost are not stored in bitpit::PathcKernel, so use all the interfaces in the geometry
		ids = getGeometryIds(ordered);
		for(const auto val: ids){
			if(this->exists(val)){
				result.push_back((*this)[val]);
			}
		}
		break;
	default:
		break;
	}
	if (squeeze)
		result.shrink_to_fit();
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
 *  - no geometry is linked
 * \return boolean coherence flag
 */
template<typename mpv_t>
bool
MimmoPiercedVector<mpv_t>::checkDataIdsCoherence(){
	if(getGeometry()==NULL) return false;
	if (this->isEmpty()) return true;
	auto ids = this->getIds();
	bool check = true;
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
 * Return a copy of the current MPV retaining data coherent with the geometry linked.
 *  If reference geometry is nullptr or none of the items are coherent with geometry return
 * an empty structure;
 * \return coherent version of current MPV.
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t>
MimmoPiercedVector<mpv_t>::resizeToCoherentDataIds(){

    MimmoPiercedVector<mpv_t> result(this->getGeometry(), m_loc);
	if(!this->getGeometry()) return result;

    switch(m_loc){

    case MPVLocation::CELL:
	{
		bitpit::PiercedVector<bitpit::Cell> & cells = m_geometry->getCells();
        result.reserve(cells.size());
        long id;
        for(auto it=this->begin(); it!=this->end(); ++it){
            long id = it.getId();
            if(cells.exists(id)){
                result.insert(id, *it);
            }
        }
	}
	break;
	case MPVLocation::INTERFACE:
	{
        bitpit::PiercedVector<bitpit::Interface> & interfaces = m_geometry->getInterfaces();
		size_t sizeInterfaces = interfaces.size();
		if(sizeInterfaces == 0){
			(*m_log)<<"Warning: Asked Data Ids Coherence in MimmoPiercedVector for INTERFACES, but linked geometry may not have them built."<<std::endl;
		}
        result.reserve(interfaces.size());
        long id;
        for(auto it=this->begin(); it!=this->end(); ++it){
            long id = it.getId();
            if(interfaces.exists(id)){
                result.insert(id, *it);
            }
        }
	}
	break;
	case MPVLocation::POINT:
	{
        bitpit::PiercedVector<bitpit::Vertex> & verts = m_geometry->getVertices();
        result.reserve(verts.size());
        long id;
        for(auto it=this->begin(); it!=this->end(); ++it){
            long id = it.getId();
            if(verts.exists(id)){
                result.insert(id, *it);
            }
        }
	}
	break;
	default:
		(*m_log)<<"NO suitable location data found to perform ids coherent resizing "<<std::endl;
		break;
	}
	return result;
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
		if(geo->getNVertices() == 0){
			(*m_log)<<"MimmoPiercedVector warning: geometry is empty"<<std::endl;
		}
		{
			this->clear();
			this->reserve(geo->getNVertices());
			m_geometry = geo;
			m_loc = loc;
			for(const auto & vertex: geo->getVertices()){
				this->insert(vertex.getId(), data);
			}
		}
		break;
	case MPVLocation::CELL :
		if(geo->getNCells() == 0){
			(*m_log)<<"MimmoPiercedVector warning: geometry is empty"<<std::endl;
		}
		{
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
			(*m_log)<<"MimmoPiercedVector warning: geometry interfaces are not built"<<std::endl;
		}
		{
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
 * The current object is a MimmoPiercedVector object with data located on MPVLocation::POINT
 * \param[in] p exponent value of inverse distance exponential weight w = (1/d^p)
 * \return cell data MimmoPiercedVector object located on MPVLocation::CELL
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::pointDataToCellData(double p){
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> cellData(geo, MPVLocation::CELL);
	MimmoPiercedVector<double> sumWeights(geo, MPVLocation::POINT);
	for (bitpit::Cell & cell : geo->getCells()){
		long idcell = cell.getId();
		std::array<double,3> center = geo->getPatch()->evalCellCentroid(idcell);
		mpv_t data{};
		bool init = false;
		for (long idvertex : cell.getVertexIds()){
			std::array<double,3> point = geo->getPatch()->getVertexCoords(idvertex);
			double weight = 1. / std::pow(norm2(center-point),p);
			if (!init){
				data = this->at(idvertex)*weight;
				sumWeights.insert(idcell, weight);
				init = true;
			}
			else{
				data = data + this->at(idvertex)*weight;
				sumWeights[idcell] = sumWeights[idcell] + weight;
			}
		}
		data = data / sumWeights[idcell];
		cellData.insert(idcell, data);
	}
	return cellData;
};

/*!
 * Cell data to Point data interpolation. Average of cell center data is set on point.
 * The current object is a MimmoPiercedVector object with data located on MPVLocation::CELL
 * \param[in] p exponent value of inverse distance exponential weight w = (1/d^p)
 * \return point data MimmoPiercedVector object located on MPVLocation::POINT
 *
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::cellDataToPointData(double p){
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> pointData(geo, MPVLocation::POINT);
	MimmoPiercedVector<int> countCells(geo, MPVLocation::POINT);
	MimmoPiercedVector<double> sumWeights(geo, MPVLocation::POINT);
	for (bitpit::Cell & cell : geo->getCells()){
		long idcell = cell.getId();
		if (this->exists(idcell)){
			std::array<double,3> center = geo->getPatch()->evalCellCentroid(idcell);
			for (long idvertex : cell.getVertexIds()){
				std::array<double,3> point = geo->getPatch()->getVertexCoords(idvertex);
				double weight = 1. / std::pow(norm2(center-point),p);
				if (!pointData.exists(idvertex)){
					pointData.insert(idvertex, this->at(idcell)*weight);
					sumWeights.insert(idvertex, weight);
				}
				else{
					pointData[idvertex] = pointData[idvertex] + this->at(idcell)*weight;
					sumWeights[idvertex] = sumWeights[idvertex] + weight;
				}
			}
		}
	}
	for (bitpit::Vertex & vertex : geo->getVertices()){
		long idvertex = vertex.getId();
		if (pointData.exists(idvertex))
			pointData[idvertex] = pointData[idvertex] / sumWeights[idvertex];
	}
	return pointData;
};

/*!
 * Cell data to Point data interpolation by using gradients on cell centers. A inverse cell volume
 * weighted extrapolation of reconstructed values is placed on points.
 * The current object is a MimmoPiercedVector object with data located on MPVLocation::CELL
 * \param[in] cellGradients gradients on cell centers
 * \return point data MimmoPiercedVector object located on MPVLocation::POINT
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::cellDataToPointData(const MimmoPiercedVector<mpv_t> & cellGradientsX, const MimmoPiercedVector<mpv_t> & cellGradientsY, const MimmoPiercedVector<mpv_t> & cellGradientsZ, bool maximum){
	double p = 3.;
	if (maximum){
		p = 0.;
	}
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> pointData(geo, MPVLocation::POINT);
	MimmoPiercedVector<double> sumWeights(geo, MPVLocation::POINT);
	for (bitpit::Cell & cell : geo->getCells()){
		long idcell = cell.getId();
		if (this->exists(idcell)){
			double volume = geo->evalCellVolume(idcell);
			double weight = 1. / std::pow(volume,p);
			std::array<double,3> center = geo->getPatch()->evalCellCentroid(idcell);
			for (long idvertex : cell.getVertexIds()){
				std::array<double,3> point = geo->getPatch()->getVertexCoords(idvertex);
				std::array<mpv_t,3> grad({cellGradientsX.at(idcell), cellGradientsY.at(idcell), cellGradientsZ.at(idcell)});
				mpv_t value;
				value = this->at(idcell) + this->at(idcell)*(grad[0]*(point-center)[0]);
				value = value + this->at(idcell)*(grad[1]*(point-center)[1]);
				value = value + this->at(idcell)*(grad[2]*(point-center)[2]);
				value = value * weight;
				if (!pointData.exists(idvertex)){
					pointData.insert(idvertex, value);
					sumWeights.insert(idvertex, weight);
				}
				else{
					if (maximum){
						if (norm2(value) > norm2(pointData[idvertex]))
							pointData[idvertex] = value;
					}
					else{
						pointData[idvertex] = pointData[idvertex] + value;
						sumWeights[idvertex] = sumWeights[idvertex] + weight;
					}
				}
			}
		}
	}
	if (!maximum){
		for (bitpit::Vertex & vertex : geo->getVertices()){
			long idvertex = vertex.getId();
			if (pointData.exists(idvertex)){
				pointData[idvertex] = pointData[idvertex] / sumWeights[idvertex];
			}
		}
	}
	return pointData;
};

/*!
 * Point data to boundary Interface data interpolation. Average of point data is set on interface center only for border interfaces.
 * The current object is a MimmoPiercedVector object with data located on MPVLocation::POINT
 * \param[in] p exponent value of inverse distance exponential weight w = (1/d^p)
 * \return boundary interface data MimmoPiercedVector object located on MPVLocation::INTERFACE
 */
template<typename mpv_t>
MimmoPiercedVector<mpv_t> MimmoPiercedVector<mpv_t>::pointDataToBoundaryInterfaceData(double p){
	MimmoObject* geo = this->getGeometry();
	MimmoPiercedVector<mpv_t> interfaceData(geo, MPVLocation::INTERFACE);
	MimmoPiercedVector<double> sumWeights(geo, MPVLocation::POINT);
	for (bitpit::Interface & interface : geo->getInterfaces()){
		//Check if interface is border
		if (interface.isBorder()){
			long idinterface = interface.getId();
			//Check if interface has at least one node presents in mimmo pierced vector point data
			std::size_t found = 0;
			for (long idvertex : interface.getVertexIds()){
				if (this->exists(idvertex))
					found++;
			}
			if (found == interface.getVertexCount()){
				//Interpolate value
				std::array<double,3> center = geo->getPatch()->evalInterfaceCentroid(idinterface);
				mpv_t data{};
				bool init = false;
				for (long idvertex : interface.getVertexIds()){
					std::array<double,3> point = geo->getPatch()->getVertexCoords(idvertex);
					double weight = 1. / std::pow(norm2(center-point),p);
					if (!init){
						if (this->exists(idvertex)){
							data = this->at(idvertex)*weight;
						}
						sumWeights.insert(idinterface, weight);
						init = true;
					}
					else{
						if (this->exists(idvertex)){
							data = data + this->at(idvertex)*weight;
						}
						sumWeights[idinterface] = sumWeights[idinterface] + weight;
					}
				}
				data = data / sumWeights[idinterface];
				interfaceData.insert(idinterface, data);
			}
		}
	}
	return interfaceData;
};

/*!
 * Copy data from another target MimmoPiercedVector of the same type.
 * - strict true: data transfer is made only on the shared ids of the structure.
 * - strict false: ids of target not present in the current MPV are copied also.
 *
 * \param[in] other target structure to copy from.
 * \param[in] strict boolean controlling the data transfer type
 * \return number of updates or brand new insertions done.
 */
template<typename mpv_t>
std::size_t
MimmoPiercedVector<mpv_t>::getDataFrom(const MimmoPiercedVector<mpv_t> & other, bool strict){
    long id;
    std::size_t counter(0);
    for(auto it= other.begin(); it != other.end(); ++it){
        id = it.getId();
        if(this->exists(id)){
            this->at(id) = *it;
            counter++;
        }else if(!strict){
            this->insert(id, *it);
            counter++;
        }
    }
    return counter;
}

/*!
 * Erase all data not contained in the target list and squeeze the structure.
 *
 * \param[in] list of ids to save from cleaning
 * \param[in] keepOrder if true keep the survivors in the same order as they appear in the original list.
 */
template<typename mpv_t>
void
MimmoPiercedVector<mpv_t>::squeezeOutExcept(const std::vector<long int> & list, bool keepOrder){

    if(keepOrder){
        // erase elements in the list then squeeze the piercedvector.
        std::unordered_set<long> maplist(list.begin(), list.end());
        MimmoPiercedVector<mpv_t> result (m_geometry, m_loc);
        result.reserve(maplist.size());
        for(auto it = this->begin(); it != this->end(); ++it){
            if(maplist.count(it.getId()) > 0){
                result.insert(it.getId(), *it);
            }
        }
        this->swap(result);
    }else{
        // if the order does not matter, create a new mpv and swap with the current
        MimmoPiercedVector<mpv_t> result (m_geometry, m_loc);
        for(long id: list){
            if(this->exists(id)){
                result.insert(id, this->at(id));
            }
        }
        this->swap(result);
    }
}


}
