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
	m_dir		= "./";
	m_local 	= false;
	m_wformat	= Short;
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

/*!It gets the PIDs of all the cells of the geometry Patch.
 * \return PIDs of the cells of geometry mesh.
 * Return vector of size 0 if none geometry linked.
 */
ivector1D
MimmoGeometry::getPID(){
	if (getGeometry() == NULL) return ivector1D(0);
	return m_pids;
};

/*!It gets the PID of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \return PID of the i-th cell of geometry mesh.
 * Return -1 if none geometry linked or out of bound of pids stored.
 */
int
MimmoGeometry::getPID(long i){
	if (getGeometry() == NULL || m_pids.size() <= i) return -1;
	return m_pids[i];
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
	//	if (read) m_write = false;
}

/*!It sets the name of directory to read/write the geometry.
 * \param[in] filename Name of directory.
 */
void
MimmoGeometry::setDir(string dir){
	m_dir = dir;
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
	//	if (write)m_read = false;
}

/*!It sets the name of file to write the geometry.
 * \param[in] filename Name of output file.
 */
void
MimmoGeometry::setWriteFilename(string filename){
	m_wfilename = filename;
}

/*!It sets the PIDs of all the cells of the geometry Patch.
 * \param[in] pids PIDs of the cells of geometry mesh.
 */
void
MimmoGeometry::setPID(ivector1D pids){
	m_pids = pids;
};

/*!It sets the PID of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \param[in] PID of the i-th cell of geometry mesh.
 * Return if out of range of number of cells of the geometry linked
 * or of the pids vector stored.
 */
void
MimmoGeometry::setPID(long i, int pid){
	if (m_pids.size() <= i && getNCells() <= i) return;
	long size = m_pids.size();
	m_pids.resize(max(size, getNCells()), 1);
	m_pids[i] = pid;
};

/*!It sets the PID of a cell of the geometry Patch.
 * \param[in] i Index of the cell of geometry mesh.
 * \param[in] PID of the i-th cell of geometry mesh.
 * Resize the pids vector stored if out of bound.
 */
void
MimmoGeometry::setPIDforce(long i, int pid){
	long size = m_pids.size();
	m_pids.resize(max(size, i), 1);
	m_pids[i] = pid;
};

/*!It sets the format to import/export .nas file.
 * \param[in] wform Format of .nas file (Short/Long).
 */
void
MimmoGeometry::setFormatNAS(WFORMAT wform){
	m_wformat = wform;
}


/*!It writes the mesh geometry on output .vtu file.
 *\return False if geometry is not linked.
 */
bool
MimmoGeometry::write(){
	if (getGeometry() == NULL) return false;
	switch(m_type){
	case FileType::STL :
		//Export STL
		//TODO patchkernel writes VTU!!!
	{
		getGeometry()->write(m_dir+"/"+m_wfilename);
		return true;
	}
	break;
	case FileType::SVTU :
		//Export Surface VTU
	{
		dvecarr3E	points = getGeometry()->getVertex();
		ivector2D	connectivity = getGeometry()->getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<3; j++){
				connectivity[i][j] = getGeometry()->getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_dir, m_wfilename, bitpit::VTKElementType::TRIANGLE, points , connectivity );
		vtk.write() ;
		return true;
	}
	break;
	case FileType::VVTU :
		//Export Volume VTU
	{
		dvecarr3E	points = getGeometry()->getVertex();
		ivector2D	connectivity = getGeometry()->getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<4; j++){
				connectivity[i][j] = getGeometry()->getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_dir, m_wfilename, bitpit::VTKElementType::VOXEL, points , connectivity );
		vtk.write() ;
		return true;
	}
	break;
	case FileType::NAS :
		//Export Nastran file
	{
		dvecarr3E	points = getGeometry()->getVertex();
		ivector2D	connectivity = getGeometry()->getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<3; j++){
				connectivity[i][j] = getGeometry()->getMapDataInv(connectivity[i][j]);
			}
		}
		NastranInterface nastran;
		nastran.setWFormat(m_wformat);
		if (m_pids.size() == connectivity.size()){
			nastran.write(m_dir,m_wfilename,points,connectivity, &m_pids);
		}
		else{
			nastran.write(m_dir,m_wfilename,points,connectivity);
		}
	}
	break;
	}
};

/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists.
 */
bool
MimmoGeometry::read(){

	//Local Instantiation of mimmo Object.
	m_local = true;

	switch(m_type){

	//Import STL
	case FileType::STL :
	{
		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);
		string name;

		int		np,	nt;
		darray3E point;
		{
			std::ifstream infile(m_dir+"/"+m_rfilename+".stl");
			bool check = infile.good();
			name = m_dir+m_rfilename+".stl";
			if (!check){
				infile.open(m_dir+m_rfilename+".STL");
				check = infile.good();
				name = m_dir+m_rfilename+".STL";
				if (!check) return false;
			}
		}
		STLObj stl(name, true);
		dvector2D V,N;
		ivector2D T;
		stl.load(np, nt, V, N, T);
		for (long ip=0; ip<np; ip++){
			point = conArray<double,3>(V[ip]);
			getGeometry()->setVertex(point);
		}
		getGeometry()->setConnectivity(&T);
		getGeometry()->cleanGeometry();
		cout << getGeometry()->getGeometry()->getCell(4244).getConnect()[0] << endl;
		cout << getGeometry()->getConnectivity(4244)[0] << endl;
	}
	break;
	case FileType::SVTU :
		//Import Surface VTU
	{
		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_dir, m_rfilename, bitpit::VTKElementType::TRIANGLE, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			getGeometry()->setVertex(Ipoints[ip]);
		}
		getGeometry()->setConnectivity(&Iconnectivity);
		getGeometry()->cleanGeometry();
	}
	break;
	case FileType::VVTU :
		//Import Volume VTU
	{
		MimmoObject* mimmo0 = new MimmoObject(2);
		setGeometry(mimmo0);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_dir, m_rfilename, bitpit::VTKElementType::TETRA, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			getGeometry()->setVertex(Ipoints[ip]);
		}
		getGeometry()->setConnectivity(&Iconnectivity);
		getGeometry()->cleanGeometry();
	}
	break;
	}
	return true;
};


/*!Execution command.
 * It reads the geometry if the condition m_read is true.
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


void NastranInterface::setWFormat(WFORMAT wf){
	wformat = wf;
}

void NastranInterface::writeKeyword(std::string key, std::ofstream& os){

	os.setf(ios_base::left);
	switch (wformat)
	{
	case Short:
	{
		os << setw(8) << key;
		break;
	}
	case Long:
	{
		os << setw(8) << string(key + '*');
		break;
	}
	}
	os.unsetf(ios_base::left);
}

void NastranInterface::writeCoord(darray3E& p, int& pointI, std::ofstream& os){
	// Fixed short/long formats:
	// 1 GRID
	// 2 ID : point ID - requires starting index of 1
	// 3 CP : co-ordinate system ID (blank)
	// 4 X1 : point x cp-ordinate
	// 5 X2 : point x cp-ordinate
	// 6 X3 : point x cp-ordinate
	// 7 CD : co-ordinate system for displacements (blank)
	// 8 PS : single point constraints (blank)
	// 9 SEID : super-element ID

	string separator_("");

	writeKeyword(string("GRID"), os);

	os << separator_;

	os.setf(ios_base::right);

	int pointI1 = pointI+1;
	writeValue(pointI1, os);
	os << separator_;
	writeValue("", os);
	os << separator_;
	writeValue(p[0], os);
	os << separator_;
	writeValue(p[1], os);
	os << separator_;
	switch (wformat)
	{
	case Short:
	{
		writeValue(p[2], os);

		os << nl;
//		os << setw(8) << p[2]
//						   << nl;
		os.unsetf(ios_base::right);

		break;
	}
	case Long:
	{
		os << nl;
		os.unsetf(ios_base::right);
		writeKeyword("", os);
		os.setf(ios_base::right);
		writeValue(p[2], os);
		os << nl;

		break;
	}
	default:
	{
		cout << "Unknown writeFormat enumeration" << endl;
		exit(999);
	}
	}

	os.unsetf(ios_base::right);
}

// NASTRAN INTERFACE //

void NastranInterface::writeFace(string faceType, ivector1D& facePts, int& nFace, ofstream& os, int PID){
	// Only valid surface elements are CTRIA3 and CQUAD4

	// Fixed short/long formats:
	// 1 CQUAD4
	// 2 EID : element ID
	// 3 PID : property element ID; default = EID (blank)
	// 4 G1 : grid point index - requires starting index of 1
	// 5 G2 : grid point index
	// 6 G3 : grid point index
	// 7 G4 : grid point index
	// 8 onwards - not used

	// For CTRIA3 elements, cols 7 onwards are not used

	WFORMAT wformat_ = wformat;
	wformat = Short;

	string separator_("");

	writeKeyword(faceType, os);

	os << separator_;

	os.setf(ios_base::right);

	nFace++;
	writeValue(nFace, os);

	os << separator_;

	writeValue(PID, os);

	int fp1;
	switch (wformat)
	{
	case Short:
	{
		for (int i=0; i< facePts.size(); i++)
		{
			fp1 = facePts[i] + 1;
			writeValue(fp1, os);
		}

		break;
	}
	case Long:
	{
		for (int i=0; i< facePts.size(); i++)
		{

			fp1 = facePts[i] + 1;
			writeValue(fp1, os);
//			if (i == 1)
//			{
//				os << nl;
//				os.unsetf(ios_base::right);
//				writeKeyword("", os);
//				os.setf(ios_base::right);
//			}
		}

		break;
	}
	default:
	{
		cout << "Unknown writeFormat enumeration" << endl;
		exit(999);
	}
	}

	os << nl;
	os.unsetf(ios_base::right);

	wformat = wformat_;

}


void NastranInterface::writeGeometry(dvecarr3E& points, ivector2D& faces, ofstream& os, ivector1D* PIDS){
	// Write points
	os << "$" << nl
			<< "$ Points" << nl
			<< "$" << nl;

	for (int pointI=0; pointI< points.size(); pointI++)
	{
		writeCoord(points[pointI], pointI, os);
	}


	// Write faces
	os << "$" << nl
			<< "$ Faces" << nl
			<< "$" << nl;

	int nFace = 0;
	bool flagpid = (PIDS != NULL);
	for (int faceI=0; faceI< faces.size(); faceI++)
	{
		const ivector1D& f = faces[faceI];
		int PID = 1;
		if (flagpid) PID = (*PIDS)[faceI];

		if (f.size() == 3)
		{
			writeFace("CTRIA3", faces[faceI], nFace, os, PID);
		}
		else if (f.size() == 4)
		{
			writeFace("CQUAD4", faces[faceI], nFace, os, PID);
		}
		else
		{
			cout << "Unknown face format" << endl;
			exit(999);
		}
	}
}


void NastranInterface::writeFooter(ofstream& os){
	int PID = 1;
	string separator_("");

	writeKeyword("PSHELL", os);

	os << separator_;

	writeValue(PID, os);

	for (int i = 0; i < 7; i++)
	{
		// Dummy values
		os << separator_;
		int uno = 1;
		writeValue(uno, os);
	}

	os << nl;
	writeKeyword("MAT1", os);
	os << separator_;

	int MID = 1;

	writeValue(MID, os);

	for (int i = 0; i < 7; i++)
	{
		// Dummy values
		os << separator_;
		writeValue("", os);
	}

	os << nl;
}



void NastranInterface::write(string& outputDir, string& surfaceName, dvecarr3E& points, ivector2D& faces, ivector1D* PIDS){


	ofstream os(outputDir +"/"+surfaceName + ".nas");
//	formatOS(os);

	os << "TITLE=MiMMO " << surfaceName << " mesh" << nl
			<< "$" << nl
			<< "BEGIN BULK" << nl;

	writeGeometry(points, faces, os, PIDS);

	writeFooter(os);

	os << "ENDDATA" << endl;

	return;
}





