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
 * \param[in] type Type of linked Patch 1 = surface (default value), 2 = volume, 3= pointCloud).
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
	m_bvTree.setPatch(m_patch);
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	m_bvTreeSupported = (m_type != 3);
	m_retrackBvTree = false;
	m_retrackKdTree = false;
}

/*!Custom constructor of MimmoObject. This constructor builds a generic patch from given vertex list and its related
 * local connectivity. Connectivities supported must be homogeneus w/ a single element type. Element available area
 * triangles or quads for surface geometries of tetrahedrons or hexahedrons for volume ones. 
 * Cloud points are always supported (providing null connectivity)
 * \param[in] type Type of linked Patch (0 = pointCloud, 1 = surface (default value), 2 = volume).
 * \param[in] vertex Coordinates of vertices of the geometry.
 * \param[in] connectivity Pointer to connectivity structure of the surface/volume mesh.
 */
MimmoObject::MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity){
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
	
	bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::UNDEFINED;
	int sizeVert, sizeCell;
	
	sizeVert = vertex.size();
	m_patch->reserveVertices(sizeVert);
	
	for(auto & vv : vertex)	addVertex(vv);
	
	if(connectivity == NULL)	m_type = 3;
	if (m_type != 3){
		int sizeConn = (*connectivity)[0].size();
		sizeCell = connectivity->size();
		
		switch(m_type){
			case 1: 
				if(sizeConn == 3)	eltype = bitpit::ElementInfo::TRIANGLE;
				if(sizeConn == 4)	eltype = bitpit::ElementInfo::QUAD;	
				break;
			case 2: 
				if(sizeConn == 4)	eltype = bitpit::ElementInfo::TETRA;
				if(sizeConn == 8)	eltype = bitpit::ElementInfo::HEXAHEDRON;	
				break;
			default: 
				 // never been reached
				break;
		}
		
		if(eltype != bitpit::ElementInfo::UNDEFINED){
			
			m_patch->reserveCells(sizeCell);
			
			for(auto & cc : *connectivity){
				
				livector1D temp(cc.size());
				int counter=0;
				
				for(auto && val : cc)	{
					temp[counter] = val; 
					++counter;
				}	
				
				addConnectedCell(temp, eltype);
			}
			
			m_pidsType.insert(0);
			
		}else{
			std::cout<<"Not supported connectivity found for MimmoObject"<<std::endl;
			std::cout<<"Proceeding as Point Cloud geometry"<<std::endl;
		}	
	}
	m_bvTree.setPatch(m_patch);
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	m_bvTreeSupported = (m_type != 3);
	m_retrackBvTree = false;
	m_retrackKdTree = false;
};

/*!Custom constructor of MimmoObject.
 * This constructor links a given patch of given type.
 * \param[in] type Type of linked Patch (0 = point cloud, 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class PatchKernel to be linked.
 */
MimmoObject::MimmoObject(int type, bitpit::PatchKernel* geometry){
	setPatch(type,geometry);
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
	clear();
	m_type 			= other.m_type;
	m_patch 		= other.m_patch;
	m_internalPatch = false;
	m_mapData		= other.m_mapData;
	m_mapCell		= other.m_mapCell;
	m_mapDataInv	= other.m_mapDataInv;
	m_mapCellInv	= other.m_mapCellInv;
	m_pidsType		= other.m_pidsType;
	
	m_bvTreeSupported = other.m_bvTreeSupported;
	m_bvTreeBuilt	= other.m_bvTreeBuilt;
	m_kdTreeBuilt   = other.m_kdTreeBuilt;

	if(m_bvTreeSupported && m_bvTreeBuilt)	buildBvTree();
	if(m_kdTreeBuilt)	buildKdTree();
	
	return *this;
};

/*!Clears the object. The pointer to the geometry is set to NULL and the
 *  mesh is deleted only if internally created.
 */
void
MimmoObject::clear(){
	m_type=1;
	if (m_internalPatch){
		delete m_patch;
	}
	m_patch = NULL;
	m_mapData.clear();
	m_mapCell.clear();
	m_mapDataInv.clear();
	m_mapCellInv.clear();
	m_pidsType.clear();
	m_bvTree.clean();
	cleanKdTree();
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	m_bvTreeSupported = true;
	m_retrackKdTree = false;
	m_retrackBvTree = false;
};

/*!Is the object empty?
 * \return True/false if the pointer to the geometry is NULL.
 */
bool
MimmoObject::isEmpty(){
	return (m_patch == NULL);
};

/*! Is the bvTree ordering supported w/ the current geometry?
 * \return True for connectivity-based meshes, false por point clouds
 */
bool
MimmoObject::isBvTreeSupported(){
	return m_bvTreeSupported;
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
MimmoObject::getVertexCoords(){
	dvecarr3E result(getNVertex());
	int  i = 0;
	for (auto & vertex : m_patch->getVertices()){
		result[i] = vertex.getCoords();
		++i;
	}
	return result;
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * If Id is not found, return a [1.e18,1.e18,1.e18] default value.
 * \param[in] i bitpit::PatchKernel ID of the vertex in geometry mesh.
 * \return Coordinates of the id-th vertex of geometry mesh.
 */
darray3E
MimmoObject::getVertexCoords(long i){
	if(!(getVertices().exists(i)))	return darray3E({{1.e18,1.e18,1.e18}});
	return 	m_patch->getVertexCoords(i);
};

/*!
 * Return reference to vertex structure of PatchKernel member inside the class
 * \return 	vertex structure of the Patch
 */
bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices(){
	return m_patch->getVertices();
}

/*!
 * Return const reference to vertex structure of PatchKernel member inside the class.
 * const overloading
 * \return 	vertex structure of the Patch
 */
const bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices() const {
	return m_patch->getVertices();
}

/*!
 * Get the compact connectivity of geometry. "Compact" means that in connectivity matrix the
 * vertex indexing is referred to the local, compact and sequential numbering of vertices, as you 
 * get them from the internal method getVertexCoords.
 * \return 	local connectivity matrix
 */
ivector2D
MimmoObject::getCompactConnectivity(){

	livector2D conn  =	getConnectivity();
	ivector2D result(conn.size());
	int counter = 0;
	for(auto & vv : conn){
		result[counter] = convertVertexIDToLocal(vv);
		++counter;
	}
	return result;
}

/*!
 * It gets the connectivity of the cells of the linked geometry. Index of vertices in
 * connectivity matrix are returned acoording to PatchKernel unique indexing.
 * \return vertex connectivity of the cells of geometry mesh
 */
livector2D
MimmoObject::getConnectivity(){
	
	if (isEmpty()) return livector2D(0);
	
	livector2D connecti(getNCells());
	int np, counter =0;
	
	for(auto & cell : getCells()){
		np = cell.getVertexCount();
		connecti[counter].resize(np);
		for (int i=0; i<np; ++i){
			connecti[counter][i] = cell.getVertex(i);
		}
		++counter;
	}
	
	return connecti;
};


/*!It gets the connectivity of a cell, with vertex Ids in PatchKernel unique indexing, of the linked geometry .
 * \param[in] i bitpit::PatchKernel ID of the cell of geometry mesh.
 * \return vertex connectivity of the id-th cell of geometry mesh, in PatchKernel indexing.
 */
livector1D
MimmoObject::getCellConnectivity(long i){
	
	if (isEmpty()) 					return livector1D(0);
	if (!(getCells().exists(i)))	return livector1D(0);	

	bitpit::Cell & cell = m_patch->getCell(i); 
	int np = cell.getVertexCount();
	livector1D connecti(np);
	
	for (int j=0; j<np; j++){
		connecti[j] = cell.getVertex(j);
	}
	return connecti;
};

/*!
 * Return reference to cell structure of PatchKernel member inside the class
 * \return 	cell structure of the Patch
 */
bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells(){
	return m_patch->getCells();
}

/*!
 * Return const reference to cell structure of PatchKernel member inside the class
 * const overloading
 * \return 	cell structure of the Patch, const
 */
const bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells()  const{
	return m_patch->getCells();
}

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
 * \return ID of target vertex. Return -1 if index is not found.
 */
long
MimmoObject::getMapData(int i){
	if(i<0 || i>=getNVertex())	return -1;
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
 * \return Index in a sequential vector of target vertex. Return -1 if index is not found.
 */
int
MimmoObject::getMapDataInv(long id){
	if(!(getVertices().exists(id)))	return -1;
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
	if(i<0 || i>=getNCells())	return -1;
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
	if(!(getCells().exists(id)))	return -1;
	return m_mapCellInv[id];
};

/*!Return the list of pid types actually present in your geometry
 * If empty list is provided, pidding is actually not supported
 */
std::unordered_set<short> 	&
MimmoObject::getPIDTypeList(){
	return m_pidsType;
};

/*!Return the list of PID associated to each cell of tessellation in compact 
 * sequential ordering
 * If empty list is provided, pidding is not supported for this geometry
 */
shivector1D 	
MimmoObject::getCompactPID() {
	if(!m_bvTreeSupported || m_pidsType.empty())	return shivector1D(0);
	shivector1D result(getNCells());
	int counter=0;
	for(auto & cell : getCells()){
		result[counter] = (short)cell.getPID();
		++counter;
	}
	return(result);
};

/*!Return the list of PID associated to each cell of tessellation as an unordered map
 * reporting as first index(key) the cell ID and as argument(value) the PID of the cell.
 * If empty list is provided, pidding is not supported for this geometry
 */
std::unordered_map<long,short> 	
MimmoObject::getPID() {
	if(!m_bvTreeSupported || m_pidsType.empty())	return std::unordered_map<long,short>();
	std::unordered_map<long,short> 	result;
	for(auto & cell : getCells()){
		result[cell.getId()] = (short) cell.getPID();
	}
	return(result);
};

/*!
 * Return if the bVTree ordering structure for simplicies is built/sync'd w/ your current geometry
 */
bool
MimmoObject::isBvTreeBuilt(){
	if (!m_bvTreeBuilt || !m_bvTreeSupported) return(false);
	return (m_bvTree.m_nelements == getNCells());
}

/*!
 * Return pointer to your BvTree internal structure
 */
BvTree*
MimmoObject::getBvTree(){
	if (!m_bvTreeSupported) return NULL;
	return &m_bvTree;
}

/*!
 * Return if the kdTree vertices ordering structure is built/sync'd w/ your current geometry
 */
bool
MimmoObject::isKdTreeBuilt(){
	if (!m_kdTreeBuilt ) return(false);
	return (m_kdTree.n_nodes == getNVertex());
}

/*!
 * Return pointer to your kdTree internal structure
 */
bitpit::KdTree<3, bitpit::Vertex, long >*
MimmoObject::getKdTree(){
	return &m_kdTree;
}

/*!
 * Return true if class trees are still usable but not sync'ed w/ the current geometry modifications
 * (as in the case of deformations, preserving number of nodes & cells if any)
 */
bool
MimmoObject::retrackTrees(){
	return (m_retrackBvTree || m_retrackKdTree);
}

/*!
 * Return the pointer to the actual class, as constant one. 
 */
const MimmoObject * 
MimmoObject::getCopy(){
	return this;
}

/*!Sets the vertices structure of the geometry Patch, clearing any previous vertex list stored.
 * \param[in] vertex vertex structure of geometry mesh.
 * \return False if no geometry is linked, not all vertices inserted or empty argument .
 */
bool
MimmoObject::setVertices(const bitpit::PiercedVector<bitpit::Vertex> & vertices){
	
	if (isEmpty() || vertices.size()==0 ) return false;
	
	m_mapData.clear();
	m_mapDataInv.clear();
	m_patch->resetVertices();
	
	int sizeVert = vertices.size();
	m_patch->reserveVertices(sizeVert);
	
	long id;
	darray3E coords;
	bool checkTot = true, check;
	for (auto && val : vertices){
		id = val.getId();
		coords = val.getCoords();
		check =  addVertex(coords, id);
		checkTot = checkTot && check;
	}	

	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	
	return checkTot;
};

/*!It adds the coordinates of one vertex of the geometry Patch.
 * If unique tag is specified for the vertex, assign it, otherwise provide itself
 * to get a unique tag for the added vertex. The latter option is default.
 * If tag is already assigned, return with unsuccessful insertion. 
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \param[in] idTag  unique ID associated to the vertex	
 * \return boolean true if the vertex is successful inserted.
 */
bool
MimmoObject::addVertex(const darray3E & vertex, const long idtag){
	if (isEmpty()) return false;
	if(idtag != bitpit::Vertex::NULL_ID && m_patch->getVertices().exists(idtag))	return false;
	long checkedID;
	
	bitpit::PatchKernel::VertexIterator it;
	if(idtag == bitpit::Vertex::NULL_ID){
		it = m_patch->addVertex(vertex);
		checkedID = it->getId();
	}else{
		it =m_patch->addVertex(vertex, idtag);
		checkedID = idtag;
	}
	
	m_mapData.push_back(checkedID);
	m_mapDataInv[checkedID] = m_mapData.size()-1;
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	return true;
};


/*!It modifies the coordinates of the vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex of geometry mesh.
 * \param[in] id bitpit::PatchKernel ID of vertex of geometry mesh to modify.
 * \return False if no geometry is linked or vertex id does not exist.
 */
bool
MimmoObject::modifyVertex(const darray3E & vertex, long id){
	if (isEmpty()) return false;
	if(!(getVertices().exists(id)))	return false;
	bitpit::Vertex &vert = m_patch->getVertex(id);
	vert.setCoords(vertex);
	if(m_bvTreeBuilt)  m_retrackBvTree = true;
	if(m_kdTreeBuilt)  m_retrackKdTree = true;
	return true;
};

/*!Sets the cell structure of the geometry Patch, clearing any previous cell list stored.
 * Does not do anything if class type is a pointcloud one (type 3);
 * \param[in] cells cell structure of geometry mesh.
 * \return False if no geometry is linked, not all cells inserted, empty argument.
 */
bool
MimmoObject::setCells(const bitpit::PiercedVector<Cell> & cells){
	if (isEmpty() || cells.size()==0 || !m_bvTreeSupported) return false;
	
	m_mapCell.clear();
	m_mapCellInv.clear();
	m_pidsType.clear();
	
	getPatch()->resetCells();
	
	int sizeCell = cells.size();
	getPatch()->reserveCells(sizeCell);
	
	long idc;
	int nVert;
	ElementInfo::Type eltype;
	
	for (const auto & cell : cells){
		
		// get ID
		idc = cell.getId();
		//check on element type
		nVert = checkCellType(cell.getType());
		//pass info to my local cell
		if(nVert > 0)	m_patch->addCell(cell, idc);
		
		short pid = (short)	cell.getPID();
		m_pidsType.insert(pid);
		m_patch->getCell(idc).setPID((int)pid);	
	}
	
	//create inverse map of cells
	setMapCell();
	m_bvTreeBuilt = false;
	return true;
};

/*!It adds one cell with its connectivity by vertex PatchKernel unique IDs, the type of cell to
 * be added and its own unique ID on the geometry Patch. If no ID is specified, assign it automatically.
 * The cell type is directly recovered from the size of local connectivity array. If connectivity dimension
 * mismatches with the type of PatchKernel internally defined into the class, return false. 
 * Does nothing if class type is a pointcloud one (type 3);
 * \param[in] connectivity	Connectivity of target cell of geometry mesh.
 * \param[in] type          type of element to be added, according to bitpit ElementInfo enum.
 * \param[in] idtag			id of the cell
 * \return False if no geometry is linked, idtag already assigned, mismatched connectivity 
 */
bool
MimmoObject::addConnectedCell(const livector1D & conn, bitpit::ElementInfo::Type type, long idtag){

	if (isEmpty() || conn.empty() || !m_bvTreeSupported) return false;
	if(idtag != bitpit::Cell::NULL_ID && m_patch->getCells().exists(idtag)) return false;
	
	int sizeElement = checkCellType(type); 
	if(sizeElement < 0)  return false;
		
	bitpit::PatchKernel::CellIterator it;
	
	livector1D conn_dum = conn;
	conn_dum.resize(sizeElement, 0);
	
	long checkedID;
	if(idtag == bitpit::Cell::NULL_ID){
		it = m_patch->addCell(type, true, conn_dum);
		checkedID = it->getId();
	}else{
		it = m_patch->addCell(type, true,conn_dum, idtag);
		checkedID = idtag;
	}
	
	m_pidsType.insert(0);		
	
	//create inverse map of cells
	m_mapCell.push_back(checkedID);
	m_mapCellInv[checkedID] = m_mapCell.size()-1;
	m_bvTreeBuilt = false;
	return true;
};

/*!It adds one cell with its connectivity by vertex PatchKernel unique IDs, the type of cell to
 * be added and its own unique ID on the geometry Patch. If no ID is specified, assign it automatically.
 * Overloading give chance to assign directly a geometry Part IDentifier (PID) to the cell.
 * The cell type is directly recovered from the size of local connectivity array. If connectivity dimension
 * mismatches with the type of PatchKernel internally defined into the class, return false. 
 * Does nothing if class type is a pointcloud one (type 3);
 * \param[in] connectivity	Connectivity of target cell of geometry mesh.
 * \param[in] type          type of element to be added, according to bitpit ElementInfo enum.
 * \param[in] PID           Part Identifier
 * \param[in] idtag			id of the cell
 * \return False if no geometry is linked, idtag already assigned, mismatched connectivity 
 */
bool
MimmoObject::addConnectedCell(const livector1D & conn, bitpit::ElementInfo::Type type, short PID, long idtag){
	
	if (isEmpty() || conn.empty() || !m_bvTreeSupported) return false;
	if(idtag != bitpit::Cell::NULL_ID && m_patch->getCells().exists(idtag)) return false;
	
	int sizeElement = checkCellType(type); 
	if(sizeElement < 0)  return false;
	
	bitpit::PatchKernel::CellIterator it;
	
	livector1D conn_dum = conn;
	conn_dum.resize(sizeElement, 0);
	
	long checkedID;
	if(idtag == bitpit::Cell::NULL_ID){
		it = m_patch->addCell(type, true, conn_dum);
		checkedID = it->getId();
	}else{
		it = m_patch->addCell(type, true,conn_dum, idtag);
		checkedID = idtag;
	}
	
	m_pidsType.insert(PID);
	
	//create inverse map of cells
	m_mapCell.push_back(checkedID);
	m_mapCellInv[checkedID] = m_mapCell.size()-1;
	m_bvTreeBuilt = false;
	return true;
};


/*!It sets the geometry Patch.
 * \param[in] type Type of linked Patch (1 = surface, 2 = volume, 3 = pointCloud).
 * \param[in] geometry Pointer to a geometry of class Patch to be linked.
 * \return False if the argument pointer is NULL or not correct type.
 */
bool
MimmoObject::setPatch(int type, PatchKernel* geometry){
	if (geometry == NULL ) return false;
	if (type<1 || type >3 ) return false;
	m_type 			= type;
	m_patch 		= geometry;
	m_internalPatch = false;
	
	setMapData();

	m_pidsType.clear();

	if(m_patch->getCellCount() != 0){
		m_bvTreeSupported = true;
		m_bvTree.setPatch(m_patch);
		
		for(auto & cell : geometry->getCells()){
			m_pidsType.insert(cell.getPID());
		}
	}	
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	m_retrackKdTree = false;
	m_retrackBvTree = false;
	return true;
};


/*!It sets the mapper external data/vertices.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::setMapData(){
	
	if (isEmpty() ) return false;
	
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
	
	if (isEmpty() || !m_bvTreeSupported) return false;
	
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


/*!It sets the PIDs of all the cells of the geometry Patch. 
 * \param[in] pids PIDs of the cells of geometry mesh, in compact sequential form. Vector lenght must be coherent with cells number in the object. 
 */
void
MimmoObject::setPID(shivector1D pids){
	if((int)pids.size() != getNCells() || !m_bvTreeSupported)	return;

	m_pidsType.clear();
	int counter = 0;
	for(auto & cell: getCells()){
		m_pidsType.insert(pids[counter]);
		cell.setPID(pids[counter]);
		++counter;
	}	
};

/*!It sets the PIDs of part of/all the cells in the geometry Patch. 
 * \param[in] pidsMap PIDs of the cells of geometry mesh, in map form. First argument is the cell  PatchKernel::ID, second argument is the PID. 
 */
void
MimmoObject::setPID(std::unordered_map<long, short> pidsMap){
	if(getNCells() == 0 || !m_bvTreeSupported)	return;
	
	m_pidsType.clear();
	auto & cells = getCells();
	for(auto & val: pidsMap){
		if(cells.exists(val.first)){
			cells[val.first].setPID(val.second);
			m_pidsType.insert(val.second);
		}	
	}	
};

/*!Sets the PID of single cell in the geometry Patch. 
 * \param[in] id id of the Cell
 * \param[in] pid pid to assign on cell
 */
void
MimmoObject::setPIDCell(long id, short pid){
	auto & cells = getCells();
	if(cells.exists(id)){
		cells[id].setPID((int)pid);
		m_pidsType.insert(pid);
	}	
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
	
	m_internalPatch = true;
	const int id = 0;
	if (m_type == 2){
		m_patch = new VolUnstructured(id, 3);
		dynamic_cast<VolUnstructured*> (m_patch)->setExpert(true);
	}else{
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*> (m_patch)->setExpert(true);
	}
	
	//copy data 
	const bitpit::PiercedVector<bitpit::Vertex> & pvert = other->getVertices();
	setVertices(pvert);
	
	m_bvTreeSupported = other->m_bvTreeSupported;
	m_bvTreeBuilt	= other->m_bvTreeBuilt;
	m_kdTreeBuilt   = other->m_kdTreeBuilt;
	
	if(m_bvTreeSupported){
		const bitpit::PiercedVector<bitpit::Cell> & pcell = other->getCells();
		setCells(pcell);
		if(m_bvTreeBuilt)	buildBvTree();
	}	
	if(m_kdTreeBuilt)	buildKdTree();
	m_retrackKdTree = false;
	m_retrackBvTree = false;
	//it's all copied(maps are update in the loops, pids if exists), trees are rebuilt and sync'ed is they must be
};

/*!It cleans the geometry Patch.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::cleanGeometry(){
	if (isEmpty()) return false;
	m_patch->deleteCoincidentVertices();
	if(m_bvTreeSupported)	m_patch->deleteOrphanVertices();
	
	setMapData();
	setMapCell();
	m_bvTreeBuilt = false;
	m_kdTreeBuilt = false;
	m_retrackKdTree = false;
	m_retrackBvTree = false;
	return true;
};

/*! Extract Vertex List from an ensamble of geometry Cells.
 *\param[in] cellList List of bitpit::PatchKernel IDs identifying cells.
 *\return List of bitpit::PatchKernel IDs of involved vertices.
 */  

livector1D MimmoObject::getVertexFromCellList(livector1D cellList){
	if(isEmpty())	return livector1D(0);
	
	livector1D result;
	set<long int> ordV;
	//get conn from each cell of the list
	for(auto && id : cellList){
		if(getCells().exists(id)){
			Cell  & cell = m_patch->getCell(id);
			long * conn = cell.getConnect();
			int nVloc = cell.getVertexCount();
			for(int i=0; i<nVloc; ++i)				ordV.insert(conn[i]);
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
 *  Unexistent vertex ID are marked in return vector as -1.
 * \param[in] vertexList List of bitpit::PatchKernel IDs identifying vertices.
 * \return List of local ids of vertices according m_mapData ordering. 
 */  
ivector1D MimmoObject::convertVertexIDToLocal(livector1D vertexList){
	
	ivector1D result(vertexList.size());
	int counter=0;
	for(auto && id : vertexList){
		result[counter]=getMapDataInv(id);
		++counter;
	}
	return result;
}

/*! Convert local ordered MimmoObject Vertex List to Vertex List of bitpit::PatchKernel IDs.
 *  Unexistent local vertex index are marked in return vector as -1.
 * \param[in] vList List of local ids of vertices according m_mapData ordering
 * \return list of bitpit::PatchKernel IDs identifying vertices.
 */  
livector1D MimmoObject::convertLocalToVertexID(ivector1D vList){
	
	livector1D result(vList.size());
	int counter=0;
	for(auto && i : vList){
		result[counter]=getMapData(i);
		++counter;
	}
	return result;
}


/*! Convert Cell List of bitpit::PatchKernel IDs to local ordered MimmoObject Cell List.
 *  Unexistent cell ID are marked in return vector as -1.
 * \param[in] cellList list of bitpit::PatchKernel IDs identifying cells.
 * \return list of local ids of cells according m_mapCell ordering.
 */  
ivector1D MimmoObject::convertCellIDToLocal(livector1D cellList){
	
	ivector1D result(cellList.size());
	
	int counter=0;
	for(auto && id : cellList){
		result[counter]=getMapCellInv(id);
		++counter;
	}
	return result;
}

/*! Convert local ordered MimmoObject Cell List to Cell List of bitpit::PatchKernel IDs.
 * Unexistent local cell index are marked in return vector as -1.
 * \param[in] cList List of local ids of cells according m_mapCell ordering
 * \return list of bitpit::PatchKernel IDs identifying cells.
 */  
livector1D MimmoObject::convertLocalToCellID(ivector1D cList){
	
	livector1D result(cList.size());
	
	int counter=0;
	for(auto && i : cList){
		result[counter]=getMapCell(i);
		++counter;
	}
	return result;
}


/*!
 * Extract ids of all vertices at mesh boundaries.
 * \return list of vertex IDs.
 */
livector1D 	MimmoObject::extractBoundaryVertexID(){
	
	if(isEmpty())	return livector1D(0);
	getPatch()->buildAdjacencies();
	
	std::unordered_set<long> container;
	std::unordered_set<long>::iterator it;
	std::unordered_set<long>::iterator itEnd;
	
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
	
	itEnd = container.end();
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
	
	if(m_pidsType.count(flag) < 1)	return livector1D(0);
	
	livector1D result(getNCells());
	int counter = 0;
	for(auto & cell : getCells()){
		if ( cell.getPID() == flag)	{
			result[counter] = cell.getId();	 
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
	result.reserve(getNCells());
	for(auto && id : flag){
		livector1D partial = extractPIDCells(id);
		result.insert(result.end(),partial.begin(), partial.end());
	}
	return(result);
};


/*!
 * Check if specified ElementType for a given cell is supported by the current PatchKernel mesh
 * internally allocated in the current class. Supported types are for superficial mesh 
 * bitpit::ElementInfo::TRIANGLE and bitpit::ElementInfo::QUAD, for volumetric mesh bitpit::ElementInfo::TETRA
 * and bitpit::ElementInfo::HEXAHEDRON
 * \param[in] type	element tyoe ti check, from bitpit::ElementInfo enum
 * \return integer with the dimension of the element supported. -1 flag the unsupported element;
 */
int MimmoObject::checkCellType(bitpit::ElementInfo::Type type){
	int check = -1;
	int patchType =getType();
	
	switch(patchType){
		case 1:
			if	(type == bitpit::ElementInfo::TRIANGLE) check = 3;
			if  (type == bitpit::ElementInfo::QUAD)		check = 4;
			break;
		case 2:
			if	(type == bitpit::ElementInfo::TETRA) 		check = 4;
			if  (type == bitpit::ElementInfo::HEXAHEDRON)	check = 8;
			break;
		default:	//do nothing
			break;
	}
	return check;
};

/*!
 * Reset and build again simplex bvTree of your geometry (if supports connectivity elements).
 *\param[in] value build the minimum leaf of the tree as a bounding box containing value elements at most.
 */
void MimmoObject::buildBvTree(int value){
	if(!m_bvTreeSupported || m_patch == NULL)	return;
	
	if (!m_bvTreeBuilt || m_retrackBvTree){
		m_bvTree.clean();
		m_bvTree.setup();
		m_bvTree.setMaxLeafSize(value);
		m_bvTree.buildTree();
		m_bvTreeBuilt = true;
	}
	
	m_retrackBvTree = false;
	return;
}

/*!
 * Reset and build again vertex kdTree of your geometry .
 */
void MimmoObject::buildKdTree(){
	if( m_patch == NULL)	return;
	long label;
	
	if (!m_kdTreeBuilt || m_retrackKdTree){
		cleanKdTree();
		m_kdTree.nodes.resize(getNVertex() + m_kdTree.MAXSTK);
		
		for(auto & val : getVertices()){
			label = val.getId();
			m_kdTree.insert(&val, label);
		}
		m_kdTreeBuilt = true;
	}
	m_retrackKdTree = false;
	return;
}

/*!
 * Clean the KdTree of the class
 */
void	MimmoObject::cleanKdTree(){
	m_kdTree.n_nodes = 0;
	m_kdTree.nodes.clear();
}

