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

#include <set>
#include "MimmoObject.hpp"
#include "Operators.hpp"
#include "bitpit.hpp"

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
		m_geometry = new VolUnstructured(id, 3);
		dynamic_cast<VolUnstructured*>(m_geometry)->setExpert(true);
	}else{
		m_geometry = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*>(m_geometry)->setExpert(true);
	}
	m_internalPatch = true;
}

/*!Custom constructor of MimmoObject. This constructor builds a generic patch from given vertex and connectivity.
 * \param[in] type Type of linked Patch (1 = surface (default value), 2 = volume).
 * \param[in] vertex Coordinates of vertices of the geometry.
 * \param[in] connectivity Pointer to connectivity strucutre of the surface/volume mesh.
 */
MimmoObject::MimmoObject(int type, dvecarr3E & vertex, ivector2D * connectivity){
	m_type = max(1,type);
	m_internalPatch = true;
	const int id = 0;
	if (m_type == 2){
//		m_geometry = new VolUnstructured(id, 3);
//		dynamic_cast<VolUnstructured*> (m_geometry)->setExpert(true);
	}else{
		m_geometry = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured*> (m_geometry)->setExpert(true);
	}
	setVertex(vertex);
	if (connectivity != NULL)
		setConnectivity(connectivity);
	setMapData();
	setMapCell();
};

/*!Custom constructor of MimmoObject.
 * This constructor links a given patch of given type.
 * \param[in] type Type of linked Patch (0 = generic (default value), 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class PatchKernel to be linked.
 */
MimmoObject::MimmoObject(int type, bitpit::PatchKernel* geometry){
	m_type 			= type;
	m_geometry 		= geometry;
	m_internalPatch = false;
	setMapData();
	setMapCell();
}

/*!Default destructor of MimmoObject.
 * It calls clear() method of MimmoObject.
 */
MimmoObject::~MimmoObject(){
	clear();
};

/*!Copy constructor of MimmoObject.
 */
MimmoObject::MimmoObject(const MimmoObject & other){
	m_type 			= other.m_type;
	m_geometry 		= other.m_geometry;
	m_internalPatch = false;
	m_mapData		= other.m_mapData;
	m_mapDataInv	= other.m_mapDataInv;
	m_mapCellInv	= other.m_mapCellInv;
};

/*!Assignement operator of MimmoObject.
 */
MimmoObject & MimmoObject::operator=(const MimmoObject & other){
	m_type 			= other.m_type;
	m_geometry 		= other.m_geometry;
	m_internalPatch = false;
	m_mapData		= other.m_mapData;
	m_mapDataInv	= other.m_mapDataInv;
	m_mapCellInv	= other.m_mapCellInv;
	return *this;
};

/*!It Clears the object. The pointer to the geometry is set to NULL and the
 *  mesh is deleted only if internally created.
 */
void
MimmoObject::clear(){
	if (m_internalPatch){
		delete m_geometry;
	}
	m_geometry = NULL;
};
/*!Is the object empty?
 * \return True/false if the pointer to the geometry is NULL.
 */
bool
MimmoObject::isEmpty(){
	return (m_geometry == NULL);
};

/*!It gets the type of the geometry Patch.
 * \return Type of geometry mesh (0 = generic (deprecated), 1 = surface, 2 = volume).
 */
int
MimmoObject::getType(){
	return m_type;
};

/*!It gets the number of vertices of the geometry Patch.
 * \return Number of vertices of geometry mesh.
 */
long
MimmoObject::getNVertex(){
	return m_geometry->getVertexCount();
};

/*!It gets the number of cells of the geometry Patch.
 * \return Number of cells of geometry mesh.
 */
long
MimmoObject::getNCells(){
	return m_geometry->getCellCount();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \return Coordinates of vertices of geometry mesh.
 */
dvecarr3E
MimmoObject::getVertex(){
	long nv = getNVertex();
	dvecarr3E result(nv);
	long i = 0;
	for (const Vertex &vertex : m_geometry->getVertices()){
		result[i++] = vertex.getCoords();
	}
	return result;
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \param[in] i Index of the vertex of geometry mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 */
darray3E
MimmoObject::getVertex(long i){
	return 	m_geometry->getVertexCoords(i);
};

/*!It gets the connectivity of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \return Connectivity of the i-th cell of geometry mesh.
 */
ivector1D
MimmoObject::getConnectivity(long i){
	if (m_geometry == NULL) return ivector1D();
	int np = m_geometry->getCell(i).getVertexCount();
	ivector1D connect(np);
	const long * connectivity = m_geometry->getCell(i).getConnect();
	for (int i=0; i<np; i++){
		connect[i] = connectivity[i];
	}
	return connect;
};

/*!It gets the connectivity of the cells of the geometry Patch.
 * \return Connectivity of the cells of geometry mesh.
 */
ivector2D
MimmoObject::getConnectivity(){
	if (m_geometry == NULL) return ivector2D();
	int np = m_geometry->getCell(0).getVertexCount();
	int nc = m_geometry->getCellCount();
	ivector2D connect(nc, ivector1D(np));
	for (int i=0; i<nc; i++){
		const long * connectivity = m_geometry->getCell(i).getConnect();
		for (int j=0; j<np; j++){
			connect[i][j] = connectivity[j];
		}
	}
	return connect;
};

/*!It gets the geometry Patch linked by Mimmo Object.
 * \return Pointer to geometry mesh.
 */
PatchKernel*
MimmoObject::getGeometry(){
	return m_geometry;
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

/*!It gets the id-th vertex index in a sequential vector.
 * \param[in] id ID of target vertex.
 * \return Index in a sequential vector of target vertex.
 */
int
MimmoObject::getMapDataInv(long id){
	return m_mapDataInv[id];
};


/*!It gets the cell ids.
 * \return Reference to inverse of Map data with cell ids.
 */
liimap&
MimmoObject::getMapCellInv(){
	return m_mapCellInv;
};


/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(dvecarr3E & vertex){
	
	if (m_geometry == NULL) return false;
	long nv = vertex.size();
	int mapsize = m_mapData.size();
	PatchKernel::VertexIterator index;
	for (long i=0; i<nv; i++){
		index = m_geometry->addVertex(vertex[i]);
		m_mapData.push_back(index->getId());
		m_mapDataInv[index->getId()] = mapsize + i;
	}
	return true;
};

/*!It adds and it sets the coordinates of one vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(darray3E & vertex){
	if (m_geometry == NULL) return false;
	PatchKernel::VertexIterator it = m_geometry->addVertex(vertex);
	m_mapData.push_back(it->getId());
	m_mapDataInv[it->getId()] = m_mapData.size();
	return true;
};

/*!It modifies the coordinates of the vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex of geometry mesh.
 * \param[in] id ID of vertex of geometry mesh to modify.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::modifyVertex(darray3E & vertex, long id){
	if (m_geometry == NULL) return false;
	Vertex &vert = m_geometry->getVertex(id);
	vert.setCoords(vertex);
	return true;
};

/*!It sets the connectivity of the cells of the geometry Patch.
 * \param[in] connectivity Connectivity of cells of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setConnectivity(ivector2D * connectivity){
	if (m_geometry == NULL) return false;
	int nv;
	if (m_type == 1) nv = 3;
	if (m_type == 2) nv = 4;  //only tetrahedra
	long nc = connectivity->size();
	long index;
	for (long i=0; i<nc; i++){
		unique_ptr<long[]> connecti = std::unique_ptr<long[]>(new long[nv]);
		for (int j=0; j<nv; j++){
			connecti[j] = (*connectivity)[i][j];
		}
		index = i;
		PatchKernel::CellIterator it;
		if (m_type == 1)  it = m_geometry->addCell(ElementInfo::TRIANGLE, true, index);
		if (m_type == 2)  it = m_geometry->addCell(ElementInfo::TETRA, true, index);
		it->setConnect(move(connecti));
		//DEBUG FORCE SET_TYPE
		if (m_type == 1)  it->setType(ElementInfo::TRIANGLE);
		if (m_type == 2)  it->setType(ElementInfo::TETRA);
	}
	
	//create inverse map of cells
	setMapCell();
	return true;
};

/*!It sets the geometry Patch.
 * \param[in] type Type of linked Patch (0 = generic (default value), 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class Patch to be linked.
 * \return False if the argument pointer is NULL.
 */
bool
MimmoObject::setGeometry(int type, PatchKernel* geometry){
	if (geometry == NULL) return false;
	m_geometry = geometry;
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
	if (m_geometry == NULL) return false;
	m_geometry->deleteCoincidentVertices();
	setMapData();
	setMapCell();
	return true;
};


/*!It sets the mapper external data/vertices.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::setMapData(){
	if (m_geometry == NULL) return false;
	long nv = getNVertex();
	m_mapData.clear();
	m_mapData.resize(nv);
	m_mapDataInv.clear();
	PatchKernel::VertexIterator it;
	PatchKernel::VertexIterator itend = m_geometry->vertexEnd();
	long i = 0;
	for (it = m_geometry->vertexBegin(); it != itend; ++it){
		m_mapData[i] = it->getId();
		m_mapDataInv[m_mapData[i]] = i;
		i++;
	}
	return true;
};

/*!It sets the mapper external data/cells.
 * \return False if the geometry member pointer is NULL.
 */
bool
MimmoObject::setMapData(){
	if (m_geometry == NULL) return false;
	m_mapCellInv.clear();
	PatchKernel::cellIterator it;
	PatchKernel::CellIterator itend = m_geometry->cellEnd();
	long i = 0;
	for (it = m_geometry->cellBegin(); it != itend; ++it){
		m_mapDataInv[it->getId()] = i;
		i++;
	}
	return true;
};

/*!It writes the mesh geometry on an output file.
 * \param[in] filename Name of the output file.
 */
void
MimmoObject::write(string filename){
	m_geometry->write(filename);
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
		Cell cell = m_geometry->getCell(*itCList);
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


