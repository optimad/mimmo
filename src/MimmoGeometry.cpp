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

#include "MimmoGeometry.hpp"

using namespace std;
using namespace mimmo;

/*!Default constructor of MimmoGeometry.
 */
MimmoGeometry::MimmoGeometry(){
	m_write		= false;
	m_filename	= "mimmoGeometry";
}

/*!Default destructor of MimmoGeometry.
 */
MimmoGeometry::~MimmoGeometry(){};

/*!Copy constructor of MimmoGeometry.
 */
MimmoGeometry::MimmoGeometry(const MimmoGeometry & other){
	*this = other;
};

/*!Assignement operator of MimmoGeometry.
 */
MimmoGeometry & MimmoGeometry::operator=(const MimmoGeometry & other){
	m_write = other.m_write;
	m_filename = other.m_filename;
	return *this;
};

/*!It gets the type of the geometry Patch.
 * \return Type of geometry mesh (0 = generic (deprecated), 1 = surface, 2 = volume).
 */
int
MimmoGeometry::getType(){
	return getGeometry()->getType();
};

/*!It gets the number of vertices of the geometry Patch.
 * \return Number of vertices of geometry mesh.
 */
long
MimmoGeometry::getNVertex(){
	return getGeometry()->getNVertex();
};

/*!It gets the number of cells of the geometry Patch.
 * \return Number of cells of geometry mesh.
 */
long
MimmoGeometry::getNCells(){
	return getGeometry()->getNCells();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \return Coordinates of vertices of geometry mesh.
 */
dvecarr3E
MimmoGeometry::getVertex(){
	return getGeometry()->getVertex();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \param[in] i Index of the vertex of geometry mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 */
darray3E
MimmoGeometry::getVertex(long i){
	return 	getGeometry()->getVertex(i);
};

/*!It gets the connectivity of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \return Connectivity of the i-th cell of geometry mesh.
 */
ivector1D
MimmoGeometry::getConnectivity(long i){
	return getGeometry()->getConnectivity(i);
};


/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setVertex(dvecarr3E & vertex){
	return getGeometry()->setVertex(vertex);
};

/*!It adds and it sets the coordinates of one vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setVertex(darray3E & vertex){
	return getGeometry()->setVertex(vertex);
};

/*!It modifies the coordinates of the vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex of geometry mesh.
 * \param[in] id ID of vertex of geometry mesh to modify.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::modifyVertex(darray3E & vertex, long id){
	return getGeometry()->modifyVertex(vertex,id);
};

/*!It sets the connectivity of the cells of the geometry Patch.
 * \param[in] connectivity Connectivity of cells of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setConnectivity(ivector2D * connectivity){
	return getGeometry()->setConnectivity(connectivity);
};

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MimmoGeometry::setWrite(bool write){
	m_write = write;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filenamen Name of output file.
 */
void
MimmoGeometry::setFilename(string filename){
	m_filename = filename;
}

/*!It writes the mesh geometry on output file.
 */
void
MimmoGeometry::write(){
	getGeometry()->write(m_filename);
};

/*!Execution command.
 * It writes the geometry if the condition m_write is true.
 */
void
MimmoGeometry::execute(){
	if (m_write) write();
}

