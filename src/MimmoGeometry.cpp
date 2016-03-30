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
#include "customOperators.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!Default constructor of MimmoGeometry.
 */
MimmoGeometry::MimmoGeometry(){
	m_name 		= "MiMMO.Geometry";
	m_read		= false;
	m_rfilename	= "mimmoGeometry";
	m_write		= false;
	m_wfilename	= "mimmoGeometry";
}

/*!Default destructor of MimmoGeometry.
 */
MimmoGeometry::~MimmoGeometry(){
	if (m_local) delete getGeometry();

};

/*!Copy constructor of MimmoGeometry.
 */
MimmoGeometry::MimmoGeometry(const MimmoGeometry & other){
	*this = other;
};

/*!Assignement operator of MimmoGeometry.
 */
MimmoGeometry & MimmoGeometry::operator=(const MimmoGeometry & other){
	m_type = other.m_type;
	m_read = other.m_read;
	m_rfilename = other.m_rfilename;
	m_write = other.m_write;
	m_wfilename = other.m_wfilename;
	return *this;
};

/*!It gets the type of the geometry Patch.
 * \return Type of geometry mesh (0 = generic (deprecated), 1 = surface, 2 = volume).
 * Return -1 if none geometry linked.
 */
int
MimmoGeometry::getType(){
	if (getGeometry() == NULL) return -1;
	return getGeometry()->getType();
};

/*!It gets the number of vertices of the geometry Patch.
 * \return Number of vertices of geometry mesh.
 * Return -1 if none geometry linked.
 */
long
MimmoGeometry::getNVertex(){
	if (getGeometry() == NULL) return -1;
	return getGeometry()->getNVertex();
};

/*!It gets the number of cells of the geometry Patch.
 * \return Number of cells of geometry mesh.
 * Return -1 if none geometry linked.
 */
long
MimmoGeometry::getNCells(){
	if (getGeometry() == NULL) return -1;
	return getGeometry()->getNCells();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \return Coordinates of vertices of geometry mesh.
 * Return vector of size 0 if none geometry linked.
 */
dvecarr3E
MimmoGeometry::getVertex(){
	if (getGeometry() == NULL) return vector<darray3E>(0);
	return getGeometry()->getVertex();
};

/*!It gets the coordinates of the vertices of the geometry Patch.
 * \param[in] i Index of the vertex of geometry mesh.
 * \return Coordinates of the i-th vertex of geometry mesh.
 * Return -9999 array if none geometry linked.
 */
darray3E
MimmoGeometry::getVertex(long i){
	if (getGeometry() == NULL) return {{-9999, -9999, -9999}};
	return 	getGeometry()->getVertex(i);
};

/*!It gets the connectivity of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \return Connectivity of the i-th cell of geometry mesh.
 * Return vector of size 0 if none geometry linked.
 */
ivector1D
MimmoGeometry::getConnectivity(long i){
	if (getGeometry() == NULL) return vector<int>(0);
	return getGeometry()->getConnectivity(i);
};


/*!It sets the coordinates of the vertices of the geometry Patch.
 * \param[in] vertex Coordinates of vertices of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setVertex(dvecarr3E & vertex){
	if (getGeometry() == NULL) return false;
	return getGeometry()->setVertex(vertex);
};

/*!It adds and it sets the coordinates of one vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex to be added to geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setVertex(darray3E & vertex){
	if (getGeometry() == NULL) return false;
	return getGeometry()->setVertex(vertex);
};

/*!It modifies the coordinates of the vertex of the geometry Patch.
 * \param[in] vertex Coordinates of vertex of geometry mesh.
 * \param[in] id ID of vertex of geometry mesh to modify.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::modifyVertex(darray3E & vertex, long id){
	if (getGeometry() == NULL) return false;
	return getGeometry()->modifyVertex(vertex,id);
};

/*!It sets the connectivity of the cells of the geometry Patch.
 * \param[in] connectivity Connectivity of cells of geometry mesh.
 * \return False if no geometry is linked.
 */
bool
MimmoGeometry::setConnectivity(ivector2D * connectivity){
	if (getGeometry() == NULL) return false;
	return getGeometry()->setConnectivity(connectivity);
};

/*!It sets the type of file to read the geometry during the execution.
 * \param[in] type Extension of file.
 */
void
MimmoGeometry::setFileType(FileType type){
	m_type = type;
}

/*!It sets the type of file to read the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 */
void
MimmoGeometry::setFileType(int type){
	m_type = static_cast<FileType>(type);

}

/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
MimmoGeometry::setRead(bool read){
	m_read = read;
}

/*!It sets the name of file to read the geometry.
 * \param[in] filename Name of input file.
 */
void
MimmoGeometry::setReadFilename(string filename){
	m_rfilename = filename;
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MimmoGeometry::setWrite(bool write){
	m_write = write;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filename Name of output file.
 */
void
MimmoGeometry::setWriteFilename(string filename){
	m_wfilename = filename;
}

/*!It writes the mesh geometry on output .vtu file.
 *\return False if geometry is not linked.
 */
bool
MimmoGeometry::write(){
	if (getGeometry() == NULL) return false;
	getGeometry()->write(m_wfilename);
	return true;
};

/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists.
 */
bool
MimmoGeometry::read(){

	m_local = true;
	//Local Instantiation of mimmo Object.
	MimmoObject* mimmo0 = new MimmoObject();
	setGeometry(mimmo0);

	switch(m_type){

	//Import STL
	case FileType::STL :
		int		np,	nt;
		darray3E point;
		{
			{
				std::ifstream infile(m_rfilename);
				bool check = infile.good();
				if (!check) return false;
			}
			STLObj stl(m_rfilename, true);
			dvector2D V,N;
			ivector2D T;
			stl.load(np, nt, V, N, T);
			for (long ip=0; ip<np; ip++){
				point = conArray<double,3>(V[ip]);
				getGeometry()->setVertex(point);
			}
			getGeometry()->setConnectivity(&T);
			getGeometry()->cleanGeometry();
		}
		break;

	}
	return true;
};


/*!Execution command.
 * It writes the geometry if the condition m_write is true.
 */
void
MimmoGeometry::execute(){
	bool check;
	if (m_read) check = read();
	if (!check){
		std::cout << "MiMMO : ERROR : file not found : "<< m_rfilename << std::endl;
		std::cout << " " << std::endl;
		exit(10);
	}
	if (m_write) check = write();
	if (!check){
		std::cout << "MiMMO : ERROR : write not done : geometry not linked " << std::endl;
		std::cout << " " << std::endl;
		exit(11);
	}
}

