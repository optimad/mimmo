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

#include "MimmoObject.hpp"
#include "Operators.hpp"
#include "bitpit.hpp"
#include <set>

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!Default constructor of MimmoObject.
 * It sets to zero/null each member/pointer.
 * \param[in] type Type of linked Patch 1 = surface (default value), 2 = volume).
 */
MimmoObject::MimmoObject(int type){
	m_type = max(type,1);
	const int id = 0;
	if (m_type == 2){
		m_patch = new VolUnstructured(id, 3);
		dynamic_cast<VolUnstructured*>(m_patch)->setExpert(true);
	}else{
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*>(m_patch)->setExpert(true);
	}
	m_internalPatch = true;
}

/*!Custom constructor of MimmoObject. This constructor builds a generic patch from given vertex and connectivity.
 * \param[in] type Type of linked Patch (1 = surface (default value), 2 = volume).
 * \param[in] vertex Coordinates of vertices of the geometry.
 * \param[in] connectivity Pointer to connectivity strucutre of the surface/volume mesh.
 */
MimmoObject::MimmoObject(int type, dvecarr3E & vertex, livector2D * connectivity){
	m_type = max(1,type);
	m_internalPatch = true;
	const int id = 0;
	if (m_type == 2){
		m_patch = new VolUnstructured(id, 3);
		dynamic_cast<VolUnstructured*> (m_patch)->setExpert(true);
	}else{
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*> (m_patch)->setExpert(true);
	}
	setVertex(vertex);
	if (connectivity != NULL)		setConnectivity(connectivity);	
};

/*!Custom constructor of MimmoObject.
 * This constructor links a given patch of given type.
 * \param[in] type Type of linked Patch (0 = generic (default value), 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class PatchKernel to be linked.
 */
MimmoObject::MimmoObject(int type, bitpit::PatchKernel* geometry){
	m_type 			= type;
	m_patch 		= geometry;
	m_internalPatch = false;

	setMapData();
	if(m_patch->getCellCount() != 0)	setMapCell();
}

/*!Default destructor of MimmoObject.
 * It calls clear() method of MimmoObject.
 */
MimmoObject::~MimmoObject(){
	clear();
};

/*!Copy constructor of MimmoObject. See MimmoObject::operator=.
 */
MimmoObject::MimmoObject(const MimmoObject & other){
	*this = other;
};

/*!
 * Assignement operator of MimmoObject. Please be careful when using it, because
 * internal patch is only copied by their pointer (soft copy). If you destroy the original MimmoObject
 * the copied class will have an internal patch pointing to nullptr.
 */
MimmoObject & MimmoObject::operator=(const MimmoObject & other){
	m_type 			= other.m_type;
	m_patch 		= other.m_patch;
	m_internalPatch = false;
	m_mapData		= other.m_mapData;
	m_mapCell		= other.m_mapCell;
	m_mapDataInv	= other.m_mapDataInv;
	m_mapCellInv	= other.m_mapCellInv;
	m_pids			= other.m_pids;
	m_pidsType		= other.m_pidsType;
	return *this;
};

/*!Clears the object. The pointer to the geometry is set to NULL and the
 *  mesh is deleted only if internally created.
 */
void
MimmoObject::clear(){
	if (m_internalPatch){
		delete m_patch;
	}
	m_patch = NULL;
	m_mapData.clear();
	m_mapCell.clear();
	m_mapDataInv.clear();
	m_mapCellInv.clear();
	m_pids.clear();
	m_pidsType.clear();
};
/*!Is the object empty?
 * \return True/false if the pointer to the geometry is NULL.
 */
bool
MimmoObject::isEmpty(){
	return (m_patch == NULL);
};

/*!It gets the type of the geometry Patch.
 * \return Type of geometry mesh (1 = surface, 2 = volume).
 */
int
MimmoObject::getType(){
	return m_type;
};

/*!It gets the number of vertices of the geometry Patch.
 * \return Number of vertices of geometry mesh.
 */
long
MimmoObject::getNVertex() const {
	return m_patch->getVertexCount();
};

/*!It gets the number of cells of the geometry Patch.
 * \return Number of cells of geometry mesh.
 */
long
MimmoObject::getNCells() const {
	return m_patch->getCellCount();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \return Coordinates of vertices of geometry mesh.
 */
dvecarr3E
MimmoObject::getVertex(){
	dvecarr3E result(getNVertex());
	int  i = 0;
	for (auto & vertex : m_patch->getVertices()){
		result[i] = vertex.getCoords();
		++i;
	}
	return result;
};

/*!It gets the coordinates of the vertices of the geometry Patch. Const overloading
 * \return Coordinates of vertices of geometry mesh.
 */
dvecarr3E
MimmoObject::getVertex() const {
	dvecarr3E result(getNVertex());
	int  i = 0;
	for (auto & vertex : m_patch->getVertices()){
		result[i] = vertex.getCoords();
		++i;
	}
	return result;
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \param[in] i bitpit::PatchKernel ID of the vertex in geometry mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 */
darray3E
MimmoObject::getVertex(long i){
	return 	m_patch->getVertexCoords(i);
};

/*!It gets the connectivity of a cell, in local class sequential indexing, of the linked geometry .
 * \param[in] i bitpit::PatchKernel ID of the cell of geometry mesh.
 * \return vertex connectivity of the i-th cell of geometry mesh, in local sequential indexing.
 */
ivector1D
MimmoObject::getConnectivity(long i){
	
	liimap & vmap = getMapDataInv();
	
	if (isEmpty()) return ivector1D(0);
	
	bitpit::Cell & cell = m_patch->getCell(i); 
	
	int np = cell.getVertexCount();
	ivector1D connecti(np);
	
	for (int i=0; i<np; i++){
		connecti[i] = vmap[cell.getVertex(i)];
	}
	return connecti;
};

/*!It gets the connectivity of the cells of the linked geometry, in local class sequential indexing
 * \return vertex connectivity of the cells of geometry mesh, in local sequential indexing.
 */
ivector2D
MimmoObject::getConnectivity(){
	
	liimap & vmap = getMapDataInv();
	
	if (isEmpty()) return ivector2D(0);
	
	ivector2D connecti(m_patch->getCellCount());
	int np, counter =0;
	
	for(auto & cell : m_patch->getCells()){
		np = cell.getVertexCount();
		connecti[counter].resize(np);
		for (int i=0; i<np; ++i){
			connecti[counter][i] = vmap[cell.getVertex(i)];
		}
		++counter;
	}
	
	return connecti;
};

/*!It gets the connectivity of the cells of the linked geometry, in local class sequential indexing.
 * Const overloading
 * \return vertex connectivity of the cells of geometry mesh, in local sequential indexing.
 */
ivector2D
MimmoObject::getConnectivity() const {
	
	const liimap & vmap = getMapDataInv();
	
	if (m_patch == NULL) return ivector2D(0);
	
	ivector2D connecti(m_patch->getCellCount());
	int np, counter =0;
	
	for(auto & cell : m_patch->getCells()){
		np = cell.getVertexCount();
		connecti[counter].resize(np);
		for (int i=0; i<np; ++i){
			liimap::const_iterator it = vmap.find(cell.getVertex(i));
			connecti[counter][i] = it->second;
		}
		++counter;
	}
	
	return connecti;
};


/*!It gets the geometry Patch linked by Mimmo Object.
 * \return Pointer to geometry mesh.
 */
PatchKernel*
MimmoObject::getPatch(){
	return m_patch;
};

/*!It gets the vertex ids.
 * \return Reference to Map data with vertex ids.
 */
livector1D&
MimmoObject::getMapData(){
	return m_mapData;
};

/*!It gets the i-th vertex id.
 * \param[in] i Index in a sequential vector of target vertex.
 * \return ID of target vertex.
 */
long
MimmoObject::getMapData(int i){
	return m_mapData[i];
};


/*!It gets the vertex ids.
 * \return Reference to inverse of Map data with vertex ids.
 */
liimap&
MimmoObject::getMapDataInv(){
	return m_mapDataInv;
};

/*!It gets the vertex ids.
 * \return Reference to inverse of Map data with vertex ids. Const overloading
 */
const liimap&
MimmoObject::getMapDataInv()const{
	return m_mapDataInv;
};

/*!It gets the id-th vertex index in a sequential vector.
 * \param[in] id ID of target vertex.
 * \return Index in a sequential vector of target vertex.
 */
int
MimmoObject::getMapDataInv(long id){
	return m_mapDataInv[id];
};



/*!It gets the cell ids.
 * \return Reference to Map Cell with cell ids.
 */
livector1D&
MimmoObject::getMapCell(){
	return m_mapCell;
};

/*!It gets the i-th cell id.
 * \param[in] i Index in a sequential vector of target cells.
 * \return ID of target cell.
 */
long
MimmoObject::getMapCell(int i){
	return m_mapCell[i];
};

/*!It gets the cell ids.
 * \return Reference to inverse of Map data with cell ids.
 */
liimap&
MimmoObject::getMapCellInv(){
	return m_mapCellInv;
};

/*!It gets the id-th cell index in a sequential vector.
 * \param[in] id ID of target cell.
 * \return Index in a sequential vector of target cell.
 */
int
MimmoObject::getMapCellInv(long id){
	return m_mapCellInv[id];
};

/*!Return the list of pid types actually present in your geometry
 * If empty list is provided, pidding is actually not supported
 */
const shivector1D 	&
MimmoObject::getPidTypeList() const{
	return m_pidsType;
};

/*!Return the list of pid associated to each cell of tessellation in compact 
 * sequential ordering
 * If empty list is provided, pidding is not supported for this geometry
 */
const shivector1D 	&
MimmoObject::getPid() const{
	return m_pids;
};

/*!
 * Return the pointer to the actual class, as constant one. 
 */
const MimmoObject * 
MimmoObject::getCopy(){
	return this;
}

//TODO is this an ADDVERTEX instead of SETVERTEX?
/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(dvecarr3E & vertex){
	
	if (isEmpty()) return false;
	int mapsize = m_mapData.size();
	PatchKernel::VertexIterator index;
	
	for (auto && val : vertex){
		index = m_patch->addVertex(val);
		m_mapData.push_back(index->getId());
		m_mapDataInv[index->getId()] = mapsize;
		++mapsize;
	}
	return true;
};

//TODO is this an ADDVERTEX instead of SETVERTEX?
/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(dvecarr3E * vertex){
	return(setVertex(*vertex));
};

/*!It resets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::resetVertex(dvecarr3E & vertex){
	if (isEmpty()) return false;
	m_mapData.clear();
	m_mapDataInv.clear();
	m_patch->resetVertices();
	return (setVertex(vertex));
};

/*!It resets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::resetVertex(dvecarr3E * vertex){
	return (resetVertex(*vertex));
};

/*!It adds and it sets the coordinates of one vertex of the geometry Patch.
 * If unique tag is specified for the vertex, assign it, otherwise provide itself
 * to get a unique tag for the added vertex. The latter option is default.
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \param[in] idTag  unique ID associated to the vertex	
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(const darray3E & vertex, long idtag){
	if (isEmpty()) return false;
	if(idtag != bitpit::Vertex::NULL_ID && m_patch->getVertices().exists(idtag))	return false;
	long checkedID;
	
	PatchKernel::VertexIterator it;
	if(idtag == bitpit::Vertex::NULL_ID){
		it =m_patch->addVertex(vertex);
		checkedID = it->getId();
	}else{
		it =m_patch->addVertex(vertex, idtag);
		checkedID = idtag;
	}
	
	m_mapData.push_back(checkedID);
	m_mapDataInv[checkedID] = m_mapData.size()-1;
	return true;
};


/*!It modifies the coordinates of the vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex of geometry mesh.
 * \param[in] id bitpit::PatchKernel ID of vertex of geometry mesh to modify.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::modifyVertex(const darray3E & vertex, long id){
	if (isEmpty()) return false;
	bitpit::Vertex &vert = m_patch->getVertex(id);
	vert.setCoords(vertex);
	return true;
};

/*!It sets/adds stock of cells with their connectivity by vertex IDS on the geometry Patch.
 * \param[in] connectivity Connectivity of cells of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setConnectivity(livector2D * connectivity){
	if (isEmpty() || connectivity == NULL) return false;
	
	int nv = (int)(*connectivity)[0].size();
	long index;
	ElementInfo::Type eltype;

	switch (m_type){
		case 1 :
			if (nv == 3) eltype = ElementInfo::TRIANGLE;
			if (nv == 4) eltype = ElementInfo::QUAD;
			break;
		case 2 :
			if (nv == 4) eltype = ElementInfo::TETRA;
			if (nv == 8) eltype = ElementInfo::QUAD;
			break;
		default: 
			return false;
			break;		
	}
	
	PatchKernel::CellIterator it;
	for (auto & conn : *connectivity){
		
		std::unique_ptr<long[]> connecti = std::unique_ptr<long[]>(new long[nv]);
		for (int j=0; j<nv; j++)	connecti[j] = conn[j];
		
		it = m_patch->addCell(eltype, true);
		it->setConnect(move(connecti));
		it->setType(eltype);
	}
	//create inverse map of cells
	setMapCell();
	return true;
};

/*!It sets/adds one cell with its connectivity by vertex IDs and its own ID 
 * on the geometry Patch. If no ID is specified, assign it automatically
 * \param[in] connectivity Connectivity of target cell of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setConnectivity(const livector1D & conn, long idtag){
	if (isEmpty()) return false;
	if(idtag != bitpit::Cell::NULL_ID && m_patch->getCells().exists(idtag)) return false;
	
	int nv = conn.size();
	long index;
	ElementInfo::Type eltype;
	
	switch (m_type){
		case 1 :
			if (nv == 3) eltype = ElementInfo::TRIANGLE;
			if (nv == 4) eltype = ElementInfo::QUAD;
			break;
		case 2 :
			if (nv == 4) eltype = ElementInfo::TETRA;
			if (nv == 8) eltype = ElementInfo::QUAD;
			break;
		default: 
			return false;
			break;		
	}
	
	PatchKernel::CellIterator it;
	
	std::unique_ptr<long[]> connecti = std::unique_ptr<long[]>(new long[nv]);
	for (int j=0; j<nv; ++j)	connecti[j] = conn[j];
	long checkedID;
	if(idtag == bitpit::Cell::NULL_ID){
		it = m_patch->addCell(eltype, true);
		checkedID = it->getId();
	}else{
		it = m_patch->addCell(eltype, true,idtag);
		checkedID = idtag;
	}
	
	it->setConnect(move(connecti));
	it->setType(eltype);
	
	//create inverse map of cells
	m_mapCell.push_back(checkedID);
	m_mapCellInv[checkedID] = m_mapCell.size()-1;
	return true;
};



/*!It sets the geometry Patch.
 * \param[in] type Type of linked Patch (1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class Patch to be linked.
 * \return False if the argument pointer is NULL or not correct type.
 */
bool
MimmoObject::setPatch(int type, PatchKernel* geometry){
	if (geometry == NULL ) return false;
	if (type<1 || type >2 ) return false;
	m_patch = geometry;
	m_type = type;
	setMapData();
	setMapCell();
	return true;
};


//TODO enrich cleaning of geometry with other useful utilities as double cells removal,
//		zero area/volume cells removal, isolated cells/vertices.

/*!It cleans the geometry Patch.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::cleanGeometry(){
	if (isEmpty()) return false;
	m_patch->deleteCoincidentVertices();
	setMapData();
	setMapCell();
	return true;
};


/*!It sets the mapper external data/vertices.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::setMapData(){
	if (isEmpty()) return false;
	long nv = getNVertex();
	m_mapData.clear();
	m_mapData.resize(nv);
	m_mapDataInv.clear();
	PatchKernel::VertexIterator it;
	PatchKernel::VertexIterator itend = m_patch->vertexEnd();
	int i = 0;
	for (it = m_patch->vertexBegin(); it != itend; ++it){
		m_mapData[i] = it->getId();
		m_mapDataInv[ it->getId()] = i;
		i++;
	}
	return true;
};

/*!It sets the mapper external data/cells.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::setMapCell(){
	if (isEmpty()) return false;
	m_mapCell.clear();
	m_mapCell.resize(getNCells());
	m_mapCellInv.clear();
	PatchKernel::CellIterator it;
	PatchKernel::CellIterator itend = m_patch->cellEnd();
	int i = 0;
	for (it = m_patch->cellBegin(); it != itend; ++it){
		m_mapCell[i] = it->getId();
		m_mapCellInv[ it->getId()] = i;
		i++;
	}
	return true;
};

/*!
 * Set your current class as a "soft" copy of the argument other.
 * All data are replaced by those provided by the argument.
 * Soft means that the m_patch member of the class is copied only by its pointer
 * and not allocated internally. 
 * \param[in] other pointer to another MimmoObject where copy from.
 */
void MimmoObject::setSOFTCopy(const MimmoObject * other){
	clear();
	*this = *other; 
};

/*!
 * Set your current class as a "hard" copy of the argument other.
 * All data are replaced by those provided by the argument.
 * Hard means that the m_patch member of the class is allocated internally
 * as an exact and independent copy of the m_patch member of the argument. 
 * \param[in] other pointer to another MimmoObject where copy from.
 */
void MimmoObject::setHARDCopy(const MimmoObject * other){
	clear();
	m_type 			= other->m_type;
	m_pids			= other->m_pids;
	m_pidsType		= other->m_pidsType;
	
	m_internalPatch = true;
	const int id = 0;
	if (m_type == 2){
		m_patch = new VolUnstructured(id, 3);
		dynamic_cast<VolUnstructured*> (m_patch)->setExpert(true);
	}else{
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*> (m_patch)->setExpert(true);
	}
	//copy data without ids
	
	const dvecarr3E points = other->getVertex();
	int counter = 0;
	long idtag;
	for(auto & val : points){
		idtag = other->m_mapData[counter];
		setVertex(val,idtag);
		++counter;
	}
	
	const ivector2D conn  = other->getConnectivity();
	
	if(!conn.empty()){
		livector1D connID;
		counter = 0;
		for(auto & val : conn){
			connID	= convertLocaltoVertexID(val);
			idtag = other->m_mapCell[counter];
			setConnectivity(connID, idtag);
			connID.clear();
			++counter;
		}
	}
	//it's all copied(maps are update in the loops)
};


/*!It sets the PIDs of all the cells of the geometry Patch.
 * \param[in] pids PIDs of the cells of geometry mesh.
 */
void
MimmoObject::setPID(shivector1D pids){
	if((int)pids.size() != getNCells())	return;
	m_pids.clear();
	m_pids = pids;
	
	//find different type of pids.
	std::unordered_set<short> map;
	map.insert(pids.begin(),pids.end());
	m_pidsType.clear();
	m_pidsType.insert(m_pidsType.end(), map.begin(), map.end());
};


/*! Extract Vertex List from an ensamble of geometry Simplicies.
 *\param[in] cellList List of bitpit::PatchKernel IDs identifying cells.
 *\return List of bitpit::PatchKernel IDs of involved vertices.
 */  

livector1D MimmoObject::getVertexFromCellList(livector1D cellList){
	
	livector1D result;
	set<long int> ordV;
	livector1D::iterator itCList;
	livector1D::iterator itCListEnd=  cellList.end();
	//get conn from each cell of the list
	for(itCList=cellList.begin(); itCList != itCListEnd; ++itCList){
		Cell cell = m_patch->getCell(*itCList);
		long * conn = cell.getConnect();
		int nVloc = cell.getVertexCount();
		for(int i=0; i<nVloc; ++i){
			ordV.insert(conn[i]);
		}
	}
	
	result.resize(ordV.size());
	set<long int>::iterator itS;
	set<long int>::iterator itSEnd=ordV.end();
	
	int counter =0;
	for(itS=ordV.begin(); itS != itSEnd; ++itS){
		result[counter] = *itS;
		++counter;
	}
	
	return(result);
}

/*! Convert Vertex List of bitpit::PatchKernel IDs to local ordered MimmoObject Vertex List.
 * \param[in] vertexList List of bitpit::PatchKernel IDs identifying vertices.
 * \return List of local ids of vertices according m_mapData ordering.
 */  
ivector1D MimmoObject::convertVertexIDtoLocal(livector1D vertexList){
	
	livector1D::iterator it;
	livector1D::iterator itEnd=vertexList.end();
	ivector1D result(vertexList.size());
	
	int counter=0;
	for(it=vertexList.begin(); it != itEnd; ++it){
		result[counter]=m_mapDataInv[*it];
		++counter;
	}
	return result;
}

/*! Convert local ordered MimmoObject Vertex List to Vertex List of bitpit::PatchKernel IDs.
 * \param[in] vList List of local ids of vertices according m_mapData ordering
 * \return list of bitpit::PatchKernel IDs identifying vertices.
 */  
livector1D MimmoObject::convertLocaltoVertexID(ivector1D vList){
	
	ivector1D::iterator it;
	ivector1D::iterator itEnd=vList.end();
	livector1D result(vList.size());
	
	int counter=0;
	for(it=vList.begin(); it != itEnd; ++it){
		result[counter]=m_mapData[*it];
		++counter;
	}
	return result;
}


/*! Convert Cell List of bitpit::PatchKernel IDs to local ordered MimmoObject Cell List.
 * \param[in] cellList list of bitpit::PatchKernel IDs identifying cells.
 * \return list of local ids of cells according m_mapCell ordering.
 */  
ivector1D MimmoObject::convertCellIDtoLocal(livector1D cellList){
	
	livector1D::iterator it;
	livector1D::iterator itEnd=cellList.end();
	ivector1D result(cellList.size());
	
	int counter=0;
	for(it=cellList.begin(); it != itEnd; ++it){
		result[counter]=m_mapCellInv[*it];
		++counter;
	}
	return result;
}

/*! Convert local ordered MimmoObject Cell List to Cell List of bitpit::PatchKernel IDs.
 * \param[in] cList List of local ids of cells according m_mapCell ordering
 * \return list of bitpit::PatchKernel IDs identifying cells.
 */  
livector1D MimmoObject::convertLocaltoCellID(ivector1D cList){
	
	ivector1D::iterator it;
	ivector1D::iterator itEnd=cList.end();
	livector1D result(cList.size());
	
	int counter=0;
	for(it=cList.begin(); it != itEnd; ++it){
		result[counter]=m_mapCell[*it];
		++counter;
	}
	return result;
}


/*!
 * Extract ids of all vertices at mesh boundaries.
 * \return list of vertex IDs.
 */
livector1D 	MimmoObject::extractBoundaryVertexID(){
	
	std::set<long> container;
	std::set<long>::iterator it;
	std::set<long>::iterator itEnd = container.end();
	
	for (const auto & cell : m_patch->getCells()){
		const long * conn = cell.getConnect();
		int size = cell.getFaceCount();
		for(int face=0; face<size; ++face){
			
			if(cell.isFaceBorder(face)){
				ivector1D list = cell.getFaceLocalConnect(face);
				for(auto && index : list ){
					container.insert(conn[index]);
				}
			}//endif
		}// end loop on face
		conn=NULL;
	}
	
	livector1D result(container.size());
	int counter = 0;
	for(it = container.begin(); it !=itEnd; ++it){
		result[counter] = *it;
		++counter;
	}
	return result;
};

/*!
 * Extract all cells by their bitpit::PatchKernel unique ID, associated to 
 * PID flag.
 * \param[in]	flag	PID for extraction
 * \return		list of cell ID marked as PID flag
 */
livector1D	MimmoObject::extractPIDCells(short flag){
	
	int size = m_pids.size();
	livector1D result(size);
	int counter = 0;
	for(int i=0; i<size; ++i){
		if (m_pids[i] == flag)	{
			result[counter] = getMapCell(i);	 
			++counter;
		}	
	}
	result.resize(counter);
	return	result;
};

/*!
 * Extract all cells by their bitpit::PatchKernel unique ID, associated to 
 * PID flag.
 * \param[in]	flag	list of PIDs for extraction
 * \return		list of cell ID marked as PID flags
 */
livector1D	MimmoObject::extractPIDCells(shivector1D flag){
	livector1D result;
	std::unordered_set<long> map;
	for(auto && id : flag){
		livector1D partial = extractPIDCells(id);
		map.insert(partial.begin(), partial.end());
	}
	result.insert(result.end(), map.begin(), map.end());
	return(result);
};
