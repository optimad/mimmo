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

using namespace std;
using namespace bitpit;

/*!Default constructor of MimmoObject.
 * It sets to zero/null each member/pointer.
 * \param[in] type Type of linked Patch 1 = surface (default value), 2 = volume).
 */
MimmoObject::MimmoObject(int type){
	m_type = max(type,1);
	const int id = 0;
	if (m_type == 2){
		m_geometry = new VolTriPatch(id);
		dynamic_cast<VolTriPatch*> (m_geometry)->setExpert(true);
	}else{
		m_geometry = new SurfTriPatch(id);
		dynamic_cast<SurfTriPatch*> (m_geometry)->setExpert(true);
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
		m_geometry = new VolTriPatch(id);
		dynamic_cast<VolTriPatch*> (m_geometry)->setExpert(true);
	}else{
		m_geometry = new SurfTriPatch(id);
		dynamic_cast<SurfTriPatch*> (m_geometry)->setExpert(true);
	}
	setVertex(vertex);
	if (connectivity != NULL)
		setConnectivity(connectivity);
};

/*!Custom constructor of MimmoObject.
 * This constructor links a given patch of given type.
 * \param[in] type Type of linked Patch (0 = generic (default value), 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class Patch to be linked.
 */
MimmoObject::MimmoObject(int type, Patch* geometry){
	m_type 			= type;
	m_geometry 		= geometry;
	m_internalPatch = false;
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
//	m_internalPatch = other.m_internalPatch;
	m_internalPatch = false;
};

/*!Assignement operator of MimmoObject.
 */
MimmoObject & MimmoObject::operator=(const MimmoObject & other){
	m_type 			= other.m_type;
	m_geometry 		= other.m_geometry;
//	m_internalPatch = other.m_internalPatch;
		m_internalPatch = false;
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
	for (const Vertex &vertex : m_geometry->vertices()){
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

/*!It gets the geometry Patch linked by Mimmo Object.
 * \return Pointer to geometry mesh.
 */
Patch*
MimmoObject::getGeometry(){
	return m_geometry;
};

/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(dvecarr3E & vertex){
	
	if (m_geometry == NULL) return false;
	long nv = vertex.size();
	Patch::VertexIterator index;
	for (long i=0; i<nv; i++){
		index = m_geometry->addVertex();
		index->setCoords(vertex[i]);
		
	}
	return true;
};

/*!It adds and it sets the coordinates of one vertex of the geometry Patch.
 * \param[in] index Index of vertex to be added to the geometry mesh.
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoObject::setVertex(int index, darray3E & vertex){
	if (m_geometry == NULL) return false;
	Patch::VertexIterator it = m_geometry->addVertex();
	it->setCoords(vertex);
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
		unique_ptr<long[]> connect = std::unique_ptr<long[]>(new long[nv]);
		for (int j=0; j<nv; j++){
			connect[j] = (*connectivity)[i][j];
		}
		index = i;
		Patch::CellIterator it;
		if (m_type == 1)  it = m_geometry->addCell(ElementInfo::TRIANGLE, true, index);
		if (m_type == 2)  it = m_geometry->addCell(ElementInfo::TETRA, true, index);
		it->setConnect(move(connect));
		//DEBUG FORCE SET_TYPE
		if (m_type == 1)  it->setType(ElementInfo::TRIANGLE);
		if (m_type == 2)  it->setType(ElementInfo::TRIANGLE);
		//
	}
};

/*!It sets the geometry Patch.
 * \param[in] type Type of linked Patch (0 = generic (default value), 1 = surface, 2 = volume).
 * \param[in] geometry Pointer to a geometry of class Patch to be linked.
 * \return False if the argument pointer is NULL.
 */
bool
MimmoObject::setGeometry(int type, Patch* geometry){
	if (geometry == NULL) return false;
	m_geometry = geometry;
	m_type = type;
};

/*!It writes the mesh geometry on an output file.
 * \param[in] filename Name of the output file.
 */
void
MimmoObject::write(string filename){
	m_geometry->write(filename);
};




