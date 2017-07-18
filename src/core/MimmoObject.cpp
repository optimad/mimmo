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

#include "MimmoObject.hpp"
#include "Operators.hpp"
#include <set>

using namespace std;
using namespace bitpit;

namespace mimmo{

/*!
* Default constructor of MimmoObject.
* It requires a int flag identifying the type of mesh meant to be created:
*  - surface unstructured mesh = 1
*  - volume unstructured mesh  = 2
*  - 3D Cloud Point            = 3
*  - 3D tessellated Curve      = 4
* 
* \param[in] type type of mesh
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
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
}

/*!
* Custom constructor of MimmoObject. 
* This constructor builds a generic mesh from given vertex list and its related
* local connectivity. The type of mesh is desumed by local cell-vertices connectivity size. 
* Connectivities supported must be homogeneus w/ a single element type, no mixed cell types are allowed. 
* Cell element types supported are:
*  - triangles or quads for surface meshes
*  - tetrahedrons or hexahedrons for volume meshes
*  - vertex for 3D point clouds
*  - lines for 3D curves
* 
* Cloud points are always supported (no connectivity field must be provided)
* Type of meshes supported are described in the default constructor MimmoObject(int type) documentation. 
* \param[in] type type of meshes.
* \param[in] vertex Coordinates of geometry vertices.
* \param[in] connectivity pointer to mesh connectivity list (optional).
*/
MimmoObject::MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity){
    m_type = max(1,type);
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
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
            case 4: 
                if(sizeConn == 2)	eltype = bitpit::ElementInfo::LINE;	
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
            (*m_log)<<"Not supported connectivity found for MimmoObject"<<std::endl;
            (*m_log)<<"Proceeding as Point Cloud geometry"<<std::endl;
        }	
    }
    m_bvTree.setPatch(m_patch);
    m_bvTreeBuilt = false;
    m_kdTreeBuilt = false;
    m_bvTreeSupported = (m_type != 3);
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
};

/*!
* Custom constructor of MimmoObject.
* This constructor builds a geometry data structure soft linking to an external bitpit::PatchKernel object.
* The mesh type needs to be specified (see default constructor MimmoObject(int type) doc).
* \param[in] type type of mesh
* \param[in] geometry pointer to a geometry of class PatchKernel to be linked.
*/
MimmoObject::MimmoObject(int type, bitpit::PatchKernel* geometry){
    m_patch = NULL;
    setPatch(type,geometry);
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
}

/*!
* Default destructor
*/
MimmoObject::~MimmoObject(){
    clear();
};

/*!
 * Copy constructor of MimmoObject. 
 */
MimmoObject::MimmoObject(const MimmoObject & other){
    m_patch = NULL;
    *this = other;
};

/*!
* Assignement operator of MimmoObject. 
* Please be careful when using it, because internal bitpit::PatchKernel structure 
* is only copied by their pointer (soft copy). If you destroy the original MimmoObject
* the copied class will have an internal patch pointing to nullptr.
* \param[in] other reference to another MimmoObject
*/
MimmoObject & MimmoObject::operator=(const MimmoObject & other){
    m_type 			= other.m_type;
    if(m_patch != NULL){
        if (m_internalPatch)    delete m_patch;
        m_patch = NULL;
    }
    m_patch 		= other.m_patch;
    m_internalPatch = false;
    
    m_mapData		= other.m_mapData;
    m_mapCell		= other.m_mapCell;
    m_mapDataInv	= other.m_mapDataInv;
    m_mapCellInv	= other.m_mapCellInv;
    m_pidsType		= other.m_pidsType;
    
    m_bvTreeSupported = other.m_bvTreeSupported;
    m_bvTreeBuilt   = other.m_bvTreeBuilt;
    m_kdTreeBuilt   = other.m_kdTreeBuilt;
    m_bvTreeSync   = other.m_bvTreeSync;
    m_kdTreeSync   = other.m_kdTreeSync;

    if(m_bvTreeSupported && m_bvTreeBuilt)	buildBvTree();
    if(m_kdTreeBuilt)	buildKdTree();

    m_AdjBuilt = other.m_AdjBuilt;
    
    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
    return *this;
};

/*!
 * It clears all contents of the class. 
 * The internal pointer to bitpit::PatchKernel geometry is set to NULL and deleted eventually,
 * if it was internally created.
 */
void
MimmoObject::clear(){
    m_type=1;
    if(m_patch != NULL){
        if (m_internalPatch)    delete m_patch;
        m_patch = NULL;
    }    
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
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;

};

/*!
 * Is the object empty?
 * \return True/false if the geometry data structure pointer is NULL.
 */
bool
MimmoObject::isEmpty(){
    if(this == NULL)	return true;
    return (m_patch == NULL);
};

/*! 
 * Is the bvTree ordering supported w/ the current geometry?
 * \return True for connectivity-based meshes, false por point clouds
 */
bool
MimmoObject::isBvTreeSupported(){
    return m_bvTreeSupported;
};

/*!
 * Return the type of mesh currently hold by the class 
 * (for type of mesh allowed see MimmoObject(int type) documentation).
 * \return integer flag for mesh type
 */
int
MimmoObject::getType(){
    return m_type;
};

/*!
 * Return the total number of vertices within the data structure.
 * \return number of mesh vertices 
 */
long
MimmoObject::getNVertex() const {
    return m_patch->getVertexCount();
};

/*!
 * Return the total number of cells within the data structure.
 * \return number of mesh cells.
 */
long
MimmoObject::getNCells() const {
    return m_patch->getCellCount();
};

/*!
 * Return the compact list of vertices hold by the class
 * \return coordinates of mesh vertices 
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

/*!
 * Gets the vertex coordinates marked with a unique label i (bitpit::PatchKernel id);
 * If i is not found, return a [1.e18,1.e18,1.e18] default value.
 * \param[in] i bitpit::PatchKernel id of the vertex in the mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 */
darray3E
MimmoObject::getVertexCoords(long i){
    if(!(getVertices().exists(i)))	return darray3E({{1.e18,1.e18,1.e18}});
    return 	m_patch->getVertexCoords(i);
};

/*!
 * Return reference to the PiercedVector structure of vertex hold 
 * by bitpit::PatchKernel class member
 * \return  Pierced Vector structure of vertices
 */
bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices(){
    return m_patch->getVertices();
}

/*!
 * Return const reference to the PiercedVector structure of vertex hold 
 * by bitpit::PatchKernel class member
 * \return  const Pierced Vector structure of vertices
 */
const bitpit::PiercedVector<bitpit::Vertex> &
MimmoObject::getVertices() const {
    return m_patch->getVertices();
}

/*!
 * Get the geometry compact connectivity.
 * "Compact" means that in the connectivity matrix the vertex indexing is referred to 
 * the local, compact and sequential numbering of vertices, as it 
 * gets them from the internal method getVertexCoords().
 * \return local connectivity list 
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
 * connectivity matrix are returned acoording to bitpit::PatchKernel unique indexing.
 * \return cell-vertices connectivity 
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

/*!
 * It gets the connectivity of a cell, with vertex id's in bitpit::PatchKernel unique indexing.
 * \param[in] i bitpit::PatchKernel ID of the cell.
 * \return i-th cell-vertices connectivity, in bitpit::PatchKernel vertex indexing.
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
 * Return reference to the PiercedVector structure of cells hold 
 * by bitpit::PatchKernel class member
 * \return  Pierced Vector structure of cells
 */
bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells(){
    return m_patch->getCells();
}

/*!
 * Return const reference to the PiercedVector structure of cells hold 
 * by bitpit::PatchKernel class member
 * \return  const Pierced Vector structure of cells
 */
const bitpit::PiercedVector<bitpit::Cell> &
MimmoObject::getCells()  const{
    return m_patch->getCells();
}

/*!
 * \return pointer to bitpit::PatchKernel structure hold by the class.
 */
PatchKernel*
MimmoObject::getPatch(){
    return m_patch;
};

/*!
 * Return the indexing vertex map, to pass from local, compact indexing 
 * to bitpit::PatchKernel unique-labeled indexing
 * \return local/unique-id map
 */
livector1D&
MimmoObject::getMapData(){
    return m_mapData;
};

/*!
 * Return the unique id label of a vertex, given its local, compact index i.
 * \param[in] i index in a sequential vector of target vertex.
 * \return unique-id of target vertex. Return -1 if index is not found.
 */
long
MimmoObject::getMapData(int i){
    if(i<0 || i>=getNVertex())	return -1;
    return m_mapData[i];
};


/*!
 * Return the indexing vertex map, to pass from bitpit::PatchKernel unique-labeled indexing
 * to local, compact indexing;
 * \return unique-id/local map
 */
liimap&
MimmoObject::getMapDataInv(){
    return m_mapDataInv;
};

/*!
 * Return the const reference to the indexing vertex map, to pass from 
 * bitpit::PatchKernel unique-labeled indexing to local, compact indexing.
 * \return unique-id/local map as const reference
 */
const liimap&
MimmoObject::getMapDataInv()const{
    return m_mapDataInv;
};

/*!
 * Return the local compact index i of a vertex, given its unique id label. 
 * \param[in] id unique-id of the vertex.
 * \return local index of the vertex. Return -1 if index is not found.
 */
int
MimmoObject::getMapDataInv(long id){
    if(!(getVertices().exists(id)))	return -1;
    return m_mapDataInv[id];
};


/*!
 * Return the indexing cell map, to pass from local, compact indexing 
 * to bitpit::PatchKernel unique-labeled indexing
 * \return local/unique-id map
 */
livector1D&
MimmoObject::getMapCell(){
    return m_mapCell;
};

/*!
 * Return the unique id label of a cell, given its local, compact index i.
 * \param[in] i index in a sequential vector of target cells.
 * \return unique-id of target cell. Return -1 if index is not found.
 */
long
MimmoObject::getMapCell(int i){
    if(i<0 || i>=getNCells())	return -1;
    return m_mapCell[i];
};

/*!
 * Return the indexing cell map, to pass from bitpit::PatchKernel unique-labeled indexing
 * to local, compact indexing;
 * \return unique-id/local map
 */
liimap&
MimmoObject::getMapCellInv(){
    return m_mapCellInv;
};

/*!
 * Return the local compact index i of a cell, given its unique id label. 
 * \param[in] id unique-id of the cell.
 * \return local index of the cell. Return -1 if index is not found.
 */
int
MimmoObject::getMapCellInv(long id){
    if(!(getCells().exists(id)))	return -1;
    return m_mapCellInv[id];
};

/*!
 * \return the list of PID types actually present in your geometry.
 * If empty list is returned, pidding is actually not supported for this geometry
 */
std::unordered_set<short> 	&
MimmoObject::getPIDTypeList(){
    return m_pidsType;
};

/*!
 * \return the list of PID associated to each cell of tessellation in compact 
 * sequential ordering
 * If empty list is returned, pidding is not supported for this geometry
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

/*!
 * \return the map of the PID (argument) associated to each cell unique-id (key) in tessellation.
 * If empty map is returned, pidding is not supported for this geometry.
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
 * \return true if the bVTree ordering structure for cells is built/synchronized 
 * with your current geometry
 */
bool
MimmoObject::isBvTreeBuilt(){
    if (!m_bvTreeBuilt || !m_bvTreeSupported) return(false);
    return (isBvTreeSync());
}

/*!
 * \return pointer to geometry BvTree internal structure
 */
BvTree*
MimmoObject::getBvTree(){
    if (!m_bvTreeSupported) return NULL;
    return &m_bvTree;
}

/*!
 * \return true if the kdTree vertices ordering structure is 
 * built/synchronized with your current geometry
 */
bool
MimmoObject::isKdTreeBuilt(){
    if (!m_kdTreeBuilt ) return(false);
    return (isKdTreeSync());
}

/*!
 * \return pointer to geometry KdTree internal structure
 */
bitpit::KdTree<3, bitpit::Vertex, long >*
MimmoObject::getKdTree(){
    return &m_kdTree;
}

/*!
 * \return true if the bVTree is synchronized
 * with your current geometry
 */
bool
MimmoObject::isBvTreeSync(){
    if (!m_bvTreeBuilt || !m_bvTreeSupported) return(false);
    return (m_bvTreeSync);
}

/*!
 * \return true if the kdTree is synchronized
 * with your current geometry
 */
bool
MimmoObject::isKdTreeSync(){
    if (!m_kdTreeBuilt) return(false);
    return (m_kdTreeSync);
}

/*!
 * \return the pointer to the actual class, as constant one. 
 */
const MimmoObject * 
MimmoObject::getCopy(){
    return this;
}

/*!
 * Set the vertices structure of the class, clearing any previous vertex list stored.
 * Be careful: any connectivity information stored in an existent cell list will be erased too.
 * The cell list will survive, but carrying no connectivity information.
 * Local/unique-id vertex maps are updated automatically. 
 * \param[in] vertices geometry vertex structure .
 * \return false if no geometry is linked, not all vertices inserted or empty argument.
 */
bool
MimmoObject::setVertices(const bitpit::PiercedVector<bitpit::Vertex> & vertices){
    
    if (vertices.empty()) return false;
    
    m_mapData.clear();
    m_mapDataInv.clear();
    m_patch->resetVertices();
    
    int sizeVert = vertices.size();
    m_patch->reserveVertices(sizeVert);
    
    long id;
    darray3E coords;
    bool checkTot = true;
    for (auto && val : vertices){
        id = val.getId();
        coords = val.getCoords();
        checkTot = checkTot && addVertex(coords, id);
    }	

    return checkTot;
};

/*!
 *It adds one vertex to the mesh.
 * If unique-id is specified for the vertex, assign it, otherwise provide itself
 * to get a unique-id for the added vertex. The latter option is the default.
 * If the unique-id is already assigned, return with unsuccessful insertion.
 * 
 * Local/unique-id vertex maps are updated automatically. 
 * 
 * 
 * \param[in] vertex vertex coordinates to be added 
 * \param[in] idtag  unique id associated to the vertex	
 * \return true if the vertex is successful inserted.
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
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    m_AdjBuilt = false;
    return true;
};


/*!
 * It modifies the coordinates of a pre-existent vertex.
 * \param[in] vertex new vertex coordinates
 * \param[in] id unique-id of the vertex meant to be modified
 * \return false if no geometry is present or vertex id does not exist.
 */
bool
MimmoObject::modifyVertex(const darray3E & vertex, long id){
    if (isEmpty()) return false;
    if(!(getVertices().exists(id)))	return false;
    bitpit::Vertex &vert = m_patch->getVertex(id);
    vert.setCoords(vertex);
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    return true;
};

/*!
 * Sets the cell structure of the geometry, clearing any previous cell list stored.
 * Does not do anything if class type is a point cloud (mesh type 3).
 *
 * Local/unique-id cell maps are updated automatically. 
 * 
 * \param[in] cells cell structure of geometry mesh.
 * \return false if no geometry is linked, not all cells are inserted or empty argument.
 */
bool
MimmoObject::setCells(const bitpit::PiercedVector<Cell> & cells){
    
    if (cells.empty() || !m_bvTreeSupported) return false;

    m_mapCell.clear();
    m_mapCellInv.clear();
    m_pidsType.clear();

    getPatch()->resetCells();

    int sizeCell = cells.size();
    getPatch()->reserveCells(sizeCell);
    
    long idc;
    int nVert;
    short pid;
    bitpit::ElementInfo::Type eltype;
    bool checkTot = true;
    livector1D connectivity;
    
    for (const auto & cell : cells){
        // get ID
        idc = cell.getId();
        //check on element type
        eltype = cell.getType();
        //check info PID 
        pid = (short)cell.getPID();
        nVert = cell.getVertexCount();
        connectivity.resize(nVert);
        auto conn = cell.getConnect();
        for(int i=0; i<nVert; ++i){
            connectivity[i] = conn[i];
        }
        checkTot = checkTot && addConnectedCell(connectivity, eltype, pid, idc);
    }

    return checkTot;
};

/*!
 * It adds one cell with its vertex-connectivity (vertex in bitpit::PatchKernel unique id's), the type of cell to
 * be added and its own unique id. If no id is specified, teh method assigns it automatically.
 * Cell type and connectivity dimension of the cell are cross-checked with the mesh type of the class: if mismatch, the method
 * does not add the cell and return false. 
 * The method does nothing, if class type is a pointcloud one (type 3).
 * 
 * Local/unique-id cell maps are updated automatically. 
 * 
 * 
 * \param[in] conn  connectivity of target cell of geometry mesh.
 * \param[in] type  type of element to be added, according to bitpit::ElementInfo enum.
 * \param[in] idtag id of the cell
 * \return false if no geometry is linked, idtag already assigned or mismatched connectivity/element type 
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
    m_bvTreeSync = false;
    m_AdjBuilt = false;
    return true;
};

/*!
 * It adds one cell with its vertex-connectivity (vertex in bitpit::PatchKernel unique id's), the type of cell to
 * be added and its own unique id. If no id is specified, teh method assigns it automatically.
 * Cell type and connectivity dimension of the cell are cross-checked with the mesh type of the class: if mismatch, the method
 * does not add the cell and return false. 
 * The method does nothing, if class type is a pointcloud one (type 3).
 * A part identifier PID mark can be associated to the cell.
 * 
 * Local/unique-id cell maps are updated automatically. 
 * 
 * \param[in] conn  connectivity of target cell of geometry mesh.
 * \param[in] type  type of element to be added, according to bitpit ElementInfo enum.
 * \param[in] PID   part identifier
 * \param[in] idtag id of the cell
 * \return false if no geometry is linked, idtag already assigned or mismatched connectivity/element type 
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
    
    setPIDCell(checkedID, PID);
    //create inverse map of cells
    m_mapCell.push_back(checkedID);
    m_mapCellInv[checkedID] = m_mapCell.size()-1;
    m_bvTreeBuilt = false;
    m_bvTreeSync = false;
    m_AdjBuilt = false;
    return true;
};

/*!
 *  Link as geometry an external bitpit::PatchKernel  data structure;
 *  If the external data structure is not compatible with the mesh type internally set, does not link anything.
 *  Through this method, the class only copy the pointer to the external geometry.
 * 
 * \param[in] type mesh type see MimmoObject(int type) constructor documentation
 * \param[in] geometry pointer to an external geoemtry data structure
 * \return false if the argument pointer is NULL or not compatible data structure
 */
bool
MimmoObject::setPatch(int type, PatchKernel* geometry){
    if (geometry == NULL ) return false;
    if (type<1 || type >4 ) return false;
    m_type 			= type;
    
    if(m_patch != NULL){
        if (m_internalPatch)    delete m_patch;
        m_patch = NULL;
    }    
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
        m_AdjBuilt = false;
    }	
    m_bvTreeBuilt = false;
    m_kdTreeBuilt = false;
    m_bvTreeSync = false;
    m_kdTreeSync = false;
    return true;
};

/*!
 * It builds the vertex map of local/unique-id indexing and its inverse.
 * \return false if no geometry is present in the class.
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

/*!
 * It builds the cell map of local/unique-id indexing and its inverse.
 * \return false if no geometry is present in the class.
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

/*!
 * Set PIDs for all geometry cells available. 
 * The PID list must be referred to the compact local/indexing of the cells in the class.
 * \param[in] pids PID list.
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

/*!
 * Set PIDs for all geometry cells available. 
 * The PID must be provided as a map with cell unique-id as key and pid associated to it as argument.
 * \param[in] pidsMap PID amp.
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

/*!
 * Set the PID of a target cell
* \param[in] id unique-id of the cell
* \param[in] pid PID to assign on cell
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
 * Set your current class as a "soft" copy of another MimmoObject.
 * All data are replaced by those provided by the argument.
 * Soft means that the m_patch member of the class is copied only by its pointer
 * and not allocated internally. 
 * \param[in] other MimmoObject class .
 */
void MimmoObject::setSOFTCopy(const MimmoObject * other){
    clear();
    *this = *other; 
};

/*!
 * Set your current class as a "hard" copy of another MimmoObject.
 * All data are replaced by those provided by the argument.
 * Hard means that the m_patch member of the class is allocated internally
 * as an exact and independent copy of the m_patch member of the argument. 
 * \param[in] other MimmoObject class.
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
    m_bvTreeSync = true;
    m_kdTreeSync = true;

    m_AdjBuilt = other->m_AdjBuilt;;

    //it's all copied(maps are update in the loops, pids if exists), trees are rebuilt and sync'ed is they must be
};

/*!
 * It cleans geometry duplicated and, in case of connected tessellations, all orphan/isolated vertices.
 * \return false if the geometry member pointer is NULL.
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
    m_bvTreeSync = false;
    m_bvTreeSync = false;
    m_AdjBuilt = false;
    return true;
};

/*! 
 * Extract vertex list from an ensamble of geometry cells.
 *\param[in] cellList list of bitpit::PatchKernel ids identifying cells.
 *\return the list of bitpit::PatchKernel ids of involved vertices.
 */  
livector1D MimmoObject::getVertexFromCellList(livector1D cellList){
    if(isEmpty())	return livector1D(0);
    if(getType() == 3)  return  livector1D(0);
    
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

/*! 
 * Extract cell list from an ensamble of geometry vertices.
 * \param[in] vertexList list of bitpit::PatchKernel IDs identifying vertices.
 * \return list of bitpit::PatchKernel IDs of involved 1-Ring cells.
 */

livector1D MimmoObject::getCellFromVertexList(livector1D vertexList){
    if(isEmpty())   return livector1D(0);
    if(getType() == 3)  return  vertexList;
    
    livector1D result;
    std::unordered_set<long int> ordV, ordC;
    ordV.insert(vertexList.begin(), vertexList.end());
    //get conn from each cell of the list
    for(auto & cell : m_patch->getCells()){
        int nVloc = cell.getVertexCount();
        int i=0;
        bool check = false;
        while(i< nVloc && !check){
            check = ordV.count(cell.getVertex(i));
            ++i;
        }
        if(check) ordC.insert(cell.getId());
    }

    result.resize(ordC.size());
    std::unordered_set<long int>::iterator itS;
    std::unordered_set<long int>::iterator itSEnd=ordC.end();

    int counter =0;
    for(itS=ordC.begin(); itS != itSEnd; ++itS){
        result[counter] = *itS;
        ++counter;
    }

    return(result);
}

/*! 
 * Convert a list of bitpit::PatchKernel vertex ids to their local/compact index value.
 * Unexistent vertex ids are marked in return vector as -1.
 * \param[in] vertexList list of bitpit::PatchKernel IDs identifying vertices.
 * \return list of local index of vertices according to compact ordering
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

/*! 
 * Convert a list of local/compact vertex index in their bitpit::PatchKernel unique id value.
 * Unexistent vertex index are marked in return vector as -1.
 * \param[in] vList list of vertex index in local/compact ordering 
 * \return list bitpit::PatchKernel ids identifying vertices.
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

/*! 
 * Convert a list of bitpit::PatchKernel cell ids to their local/compact index value.
 * Unexistent cell ids are marked in return vector as -1.
 * \param[in] cellList list of bitpit::PatchKernel IDs identifying cells.
 * \return list of local index of cells according to compact ordering
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

/*! 
 * Convert a list of local/compact cell index in their bitpit::PatchKernel unique id value.
 * Unexistent cell index are marked in return vector as -1.
 * \param[in] cList list of cell index in local/compact ordering 
 * \return list bitpit::PatchKernel ids identifying cells.
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
 * Extract vertices at the mesh boundaries, if any.
 * \return list of vertex unique-ids.
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
 * Extract all cells marked with a target PID flag.
 * \param[in]	flag	PID for extraction
 * \return		list of cells as unique-ids
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
 * Extract all cells marked with a series of target PIDs.
 * \param[in]   flag    list of PID for extraction
 * \return      list of cells as unique-ids
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
 * Check if a specified cell topology is currently supported by the class, and return the typical size of its connectivity.
 * 
 * \param[in] type  cell type to check, as bitpit::ElementInfo enum
 * \return integer with the dimension of the element supported. -1 flag the unsupported element;
 */
int MimmoObject::checkCellType(bitpit::ElementInfo::Type type){
    int check = -1;
    int patchType =getType();
    
    switch(patchType){
        case 1:
            if  (type == bitpit::ElementInfo::TRIANGLE)     check = 3;
            if  (type == bitpit::ElementInfo::QUAD)         check = 4;
            break;
        case 2:
            if  (type == bitpit::ElementInfo::TETRA)        check = 4;
            if  (type == bitpit::ElementInfo::PYRAMID)      check = 5;
            if  (type == bitpit::ElementInfo::WEDGE)        check = 6;
            if  (type == bitpit::ElementInfo::HEXAHEDRON)   check = 8;
            break;
        case 3:
            if (type == bitpit::ElementInfo::VERTEX)        check = 1;
            break;
        case 4:
            if	(type == bitpit::ElementInfo::LINE)         check = 2;
            break;
        default:	//do nothing
            break;
    }
    return check;
};

/*!
 * Evaluate axis aligned bounding box of the current MimmoObject 
 * \param[out] pmin lowest bounding box point
 * \param[out] pmax highest bounding box point
 */
void MimmoObject::getBoundingBox(std::array<double,3> & pmin, std::array<double,3> & pmax){
    if(m_patch == NULL)   return;
    m_patch->getBoundingBox(pmin,pmax);
    return;
}

/*!
 * Reset and build again simplex bvTree of your geometry (if supports connectivity elements).
 *\param[in] value build the minimum leaf of the tree as a bounding box containing value elements at most.
 */
void MimmoObject::buildBvTree(int value){
    if(!m_bvTreeSupported || m_patch == NULL)	return;
    
    if (!m_bvTreeBuilt || !m_bvTreeSync){
        m_bvTree.clean();
        m_bvTree.setup();
        m_bvTree.setMaxLeafSize(value);
        m_bvTree.buildTree();
        m_bvTreeBuilt = true;
        m_bvTreeSync = true;
    }
    return;
}

/*!
 * Reset and build again vertex kdTree of your geometry.
 */
void MimmoObject::buildKdTree(){
    if( m_patch == NULL)	return;
    long label;
    
    if (!m_kdTreeBuilt || !m_bvTreeSync){
        cleanKdTree();
        m_kdTree.nodes.resize(getNVertex() + m_kdTree.MAXSTK);
        
        for(auto & val : getVertices()){
            label = val.getId();
            m_kdTree.insert(&val, label);
        }
        m_kdTreeBuilt = true;
        m_bvTreeSync = true;
    }
    return;
}

/*!
 * Clean the KdTree of the class
 */
void	MimmoObject::cleanKdTree(){
    m_kdTree.n_nodes = 0;
    m_kdTree.nodes.clear();
}

/*!
 * \return true if cell-cell adjacency is built for your current mesh.
 */
bool MimmoObject::areAdjacenciesBuilt(){
    
    return(m_AdjBuilt);

    /*
    bool check = true;
    
    auto itp = getCells().cbegin();
    
    int cAdj = itp->getAdjacencyCount();
    const long * adj = itp->getAdjacencies();
    
    for(int i=0; i<cAdj; ++i){
        check = check && (adj[i] == -1);
    }
    
    return !check;
    */
};

/*!
 * \return false if your mesh has open edges/faces. True otherwise
 */
bool MimmoObject::isClosedLoop(){
    
    if(!areAdjacenciesBuilt())	buildAdjacencies();
    bool check = true;
    
    auto itp = getCells().cbegin();
    auto itend = getCells().cend();
    
    while(itp != itend && check){
        
        int cAdj = itp->getAdjacencyCount();
        const long * adj = itp->getAdjacencies();
        
        for(int i=0; i<cAdj; ++i){
            check = check && (adj[i] != -1);
        }
        
        itp++;
    }
    
    return check;
};

/*!
 * Force the class to build cell-cell adjacency connectivity.
 */
void MimmoObject::buildAdjacencies(){
    getPatch()->buildAdjacencies();
    m_AdjBuilt = true;
};

/*!
 * Desume Element type of your current mesh. 
 * Please note MimmoObject is handling meshes with homogeneous elements.
 * Return undefined type for unexistent or unsupported element, or mixed element type connectivity.
 * \return cell type hold by the mesh
 */
bitpit::VTKElementType	MimmoObject::desumeElement(){
    bitpit::VTKElementType result = bitpit::VTKElementType::UNDEFINED;
    
    if(getPatch() == NULL)	return result;	
    livector1D conn;
    switch(m_type){
        case	1:
            if(getNCells() == 0) 		return result;
            conn = getCellConnectivity((*(getCells().begin())).getId());
            if(conn.size() == 3)		result = bitpit::VTKElementType::TRIANGLE;
            if(conn.size() == 4)		result = bitpit::VTKElementType::QUAD;
            break;
        case	2:
            if(getNCells() == 0) 		return result;
            conn = getCellConnectivity((*(getCells().begin())).getId());
            if(conn.size() == 4)		result = bitpit::VTKElementType::TETRA;
            if(conn.size() == 8)		result = bitpit::VTKElementType::HEXAHEDRON;
            break;
        case	3:
            result = bitpit::VTKElementType::VERTEX;
            break;
        case	4:
            result = bitpit::VTKElementType::LINE;
            break;
        default : 
            break;
    }
    
    return result;
};


}



