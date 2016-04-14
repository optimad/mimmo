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
	m_rtype		= STL;
	m_read		= false;
	m_rfilename	= "mimmoGeometry";
	m_wtype		= STL;
	m_write		= false;
	m_wfilename	= "mimmoGeometry";
	m_rdir		= "./";
	m_wdir		= "./";
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
	m_rtype = other.m_rtype;
	m_wtype = other.m_wtype;
	m_read = other.m_read;
	m_rfilename = other.m_rfilename;
	m_write = other.m_write;
	m_wfilename = other.m_wfilename;
	m_rdir = other.m_rdir;
	m_wdir = other.m_wdir;
	m_local = other.m_local;
	m_wformat = other.m_wformat;
	m_pids = other.m_pids;
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
MimmoGeometry::setReadFileType(FileType type){
	m_rtype = type;
}

/*!It sets the type of file to read the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 */
void
MimmoGeometry::setReadFileType(int type){
	m_rtype = static_cast<FileType>(type);
}

/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
MimmoGeometry::setRead(bool read){
	m_read = read;
}

/*!It sets the name of directory to read the geometry.
 * \param[in] dir Name of directory.
 */
void
MimmoGeometry::setReadDir(string dir){
	m_rdir = dir;
}

/*!It sets the name of file to read the geometry.
 * \param[in] filename Name of input file.
 */
void
MimmoGeometry::setReadFilename(string filename){
	m_rfilename = filename;
}


/*!It sets the type of file to write the geometry during the execution.
 * \param[in] type Extension of file.
 */
void
MimmoGeometry::setWriteFileType(FileType type){
	m_wtype = type;
}

/*!It sets the type of file to write the geometry during the execution.
 * \param[in] type Label of file type (0 = bynary STL).
 */
void
MimmoGeometry::setWriteFileType(int type){
	m_wtype = static_cast<FileType>(type);
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MimmoGeometry::setWrite(bool write){
	m_write = write;
}

/*!It sets the name of directory to write the geometry.
 * \param[in] dir Name of directory.
 */
void
MimmoGeometry::setWriteDir(string dir){
	m_wdir = dir;
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

/*!It sets the format to export .nas file.
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

	switch(m_wtype){
	case FileType::STL :
		//Export STL
		//TODO patchkernel writes VTU update with new bitpit and use writeSTL!!!
	{
		getGeometry()->write(m_wdir+"/"+m_wfilename);
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
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TRIANGLE, points , connectivity );
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
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TETRA, points , connectivity );
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
			nastran.write(m_wdir,m_wfilename,points,connectivity, &m_pids);
		}
		else{
			nastran.write(m_wdir,m_wfilename,points,connectivity);
		}
	}
	break;
<<<<<<< HEAD
	//Export ascii OpenFOAM point cloud
	case FileType::OFP :
	{

		dvecarr3E points = getGeometry()->getVertex();
		writeOFP(m_wdir, m_wfilename, points);

	}
	break;
=======
>>>>>>> be53bd32d1e0802a4bac626aeaeedeaaa46d334c
	}
};

/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists.
 */
bool
MimmoGeometry::read(){

	//Local Instantiation of mimmo Object.
	m_local = true;

	switch(m_rtype){

	//Import STL
	case FileType::STL :
	{
		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);
		string name;

		int		np,	nt;
		darray3E point;
		{
			std::ifstream infile(m_rdir+"/"+m_rfilename+".stl");
			bool check = infile.good();
			name = m_rdir+"/"+m_rfilename+".stl";
			if (!check){
				infile.open(m_rdir+m_rfilename+".STL");
				check = infile.good();
				name = m_rdir+"/"+m_rfilename+".STL";
				if (!check) return false;
			}
		}

		ifstream in(name);
		string	ss, sstype;
		getline(in,ss);
		stringstream ins;
		ins << ss;
		ins >> sstype;
		bool binary = true;
		if (sstype == "solid" || sstype == "SOLID") binary = false;
		in.close();
		STLObj stl(name, binary);

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
	case FileType::SVTU :
		//Import Surface VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TRIANGLE, Ipoints, Iconnectivity );
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

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject* mimmo0 = new MimmoObject(2);
		setGeometry(mimmo0);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		//TODO TETRA OR VOXEL ELEMENTS?!?
		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TETRA, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			getGeometry()->setVertex(Ipoints[ip]);
		}
		getGeometry()->setConnectivity(&Iconnectivity);
		getGeometry()->cleanGeometry();
	}
	break;
	case FileType::NAS :
		//Import Surface NAS
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".nas");
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		NastranInterface nastran;
		nastran.setWFormat(m_wformat);

		nastran.read(m_rdir, m_rfilename, Ipoints, Iconnectivity, m_pids );

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			getGeometry()->setVertex(Ipoints[ip]);
		}
		getGeometry()->setConnectivity(&Iconnectivity);
		getGeometry()->cleanGeometry();
	}
	break;
	//Import ascii OpenFOAM point cloud
	case FileType::OFP :
	{
		//It uses surface type of mimmo object for the moment
		MimmoObject* mimmo0 = new MimmoObject(1);
		setGeometry(mimmo0);

		std::ifstream infile(m_rdir+"/"+m_rfilename);
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		dvecarr3E	Ipoints;
		readOFP(m_rdir, m_rfilename, Ipoints);

		int np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			getGeometry()->setVertex(Ipoints[ip]);
			cout << Ipoints[ip] << endl;
		}
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


//===============================//
//====== OFOAM INTERFACE ======//
//===============================//

void MimmoGeometry::readOFP(string& inputDir, string& surfaceName, dvecarr3E& points){

	ifstream is(inputDir +"/"+surfaceName);

	points.clear();
	int ip = 0;
	int np;
	darray3E point;
	string sread;
	char par;

	for (int i=0; i<18; i++){
		getline(is,sread);
		cout << sread << endl;
	}
	is >> np;
	getline(is,sread);
	getline(is,sread);

	points.resize(np);
	while(!is.eof() && ip<np){
		is.get(par);
		for (int i=0; i<3; i++) is >> point[i];
		is.get(par);
		getline(is,sread);
		points[ip] = point;
		ip++;
	}
	is.close();
	return;

}

void MimmoGeometry::writeOFP(string& outputDir, string& surfaceName, dvecarr3E& points){

	ofstream os(outputDir +"/"+surfaceName);
<<<<<<< HEAD
	char nl = '\n';

	string separator(" ");
	string parl("("), parr(")");
	string hline;

	hline = "/*--------------------------------*- C++ -*----------------------------------*\\" ;
	os << hline << nl;
	hline = "| =========                 |                                                 |";
	os << hline << nl;
	hline = "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |";
	os << hline << nl;
	hline = "|  \\\\    /   O peration     | Version:  2.4.x                                 |";
	os << hline << nl;
	hline = "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |";
	os << hline << nl;
	hline = "|    \\\\/     M anipulation  |                                                 |";
	os << hline << nl;
	hline = "\\*---------------------------------------------------------------------------*/";
	os << hline << nl;
	hline = "FoamFile";
	os << hline << nl;
	hline = "{";
	os << hline << nl;
	hline = "    version     2.0;";
	os << hline << nl;
	hline = "    format      ascii;";
	os << hline << nl;
	hline = "    class       vectorField;";
	os << hline << nl;
	hline = "    location    \"constant/polyMesh\";";
	os << hline << nl;
	hline = "    object      points;";
	os << hline << nl;
	hline = "}";
	os << hline << nl;
	hline = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
	os << hline << nl;
	os << nl;
	os << nl;
	int np = points.size();
	os << np << nl;
	os << parl << nl;
	for (int i=0; i<np; i++){
		os << parl;
		for (int j=0; j<2; j++){
			os << points[i][j] << separator;
		}
		os << points[i][2] << parr << nl;
	}
	os << parr << nl;
	os << nl;
	os << nl;
	hline = "// ************************************************************************* //";
	os << hline << nl;

	os.close();
=======

	string separator(" ");
	string par("(");
	string hline;

	hline = "/*--------------------------------*- C++ -*----------------------------------*\";
			hline = "| =========                 |                                                 |";
	hline = "	| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |";
	hline = "|  \\    /   O peration     | Version:  2.4.x                                 |";
	hline = "|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |";
	hline = "|    \\/     M anipulation  |                                                 |";
	hline = "\*---------------------------------------------------------------------------*/";
	hline = "FoamFile";
	hline = "{";
	hline = "	    version     2.0;";
	hline = "	    format      ascii;";
	hline = "	    class       vectorField;";
	hline = "	    location    "constant/polyMesh";";
	hline = "	    object      points;";
	hline = "	}";
	hline = "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";




	int ip = 0;
	int np = points.size();
	darray3E point;
	string sread;
	char par;




	for (int pointI=0; pointI< points.size(); pointI++)
	{
		writeCoord(points[pointI], pointI, os);
	}



	for (int i=0; i<18; i++){
		getline(is,sread);
		cout << sread << endl;
	}
	is >> np;
	getline(is,sread);
	getline(is,sread);

	points.resize(np);
	while(!is.eof() && ip<np){
		is.get(par);
		for (int i=0; i<3; i++) is >> point[i];
		is.get(par);
		getline(is,sread);
		points[ip] = point;
		ip++;
	}
	is.close();
>>>>>>> be53bd32d1e0802a4bac626aeaeedeaaa46d334c
	return;

}

//===============================//
//====== NASTRAN INTERFACE ======//
//===============================//

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
	os << "TITLE=MiMMO " << surfaceName << " mesh" << nl
			<< "$" << nl
			<< "BEGIN BULK" << nl;
	writeGeometry(points, faces, os, PIDS);
	writeFooter(os);
	os << "ENDDATA" << endl;
	return;
}

//========READ====//

void NastranInterface::read(string& inputDir, string& surfaceName, dvecarr3E& points, ivector2D& faces, ivector1D& PIDS){

	ifstream is(inputDir +"/"+surfaceName + ".nas");

	points.clear();
	faces.clear();
	PIDS.clear();

	int ipoint = 0;
	int iface  = 0;
	darray3E point;
	ivector1D face(3);
	map<int,int> mapnodes;
	int id;
	int pid;
	string sread;
	string ssub = trim(sread.substr(0,8));

	while(!is.eof()){
		while(ssub != "GRID" &&
				ssub != "GRID*" &&
				ssub != "CTRIA3" &&
				ssub != "CTRIA3*" &&
				!is.eof()){

			getline(is,sread);
			ssub = trim(sread.substr(0,8));

		}
		if(ssub == "GRID"){
			id = stoi(sread.substr(8,8));
			mapnodes[id] = ipoint;
			point[0] = stod(convertVertex(trim(sread.substr(24,8))));
			point[1] = stod(convertVertex(trim(sread.substr(32,8))));
			point[2] = stod(convertVertex(trim(sread.substr(40,8))));
			points.push_back(point);
			ipoint++;
			getline(is,sread);
			ssub = trim(sread.substr(0,8));
		}
		else if(ssub == "GRID*"){
			id = stoi(sread.substr(16,16));
			mapnodes[id] = ipoint;
			point[0] = stod(convertVertex(trim(sread.substr(48,16))));
			point[1] = stod(convertVertex(trim(sread.substr(64,16))));
			point[2] = stod(convertVertex(trim(sread.substr(80,16))));
			points.push_back(point);
			ipoint++;
			getline(is,sread);
			ssub = trim(sread.substr(0,8));
		}
		else if(ssub == "CTRIA3"){
			pid = stoi(sread.substr(16,8));
			face[0] = mapnodes[stoi(sread.substr(24,8))];
			face[1] = mapnodes[stoi(sread.substr(32,8))];
			face[2] = mapnodes[stoi(sread.substr(40,8))];
			faces.push_back(face);
			PIDS.push_back(pid);
			iface++;
			getline(is,sread);
			ssub = trim(sread.substr(0,8));
		}
		else if(ssub == "CTRIA3*"){
			pid = stoi(sread.substr(32,16));
			face[0] = mapnodes[stoi(sread.substr(48,16))];
			face[1] = mapnodes[stoi(sread.substr(64,16))];
			face[2] = mapnodes[stoi(sread.substr(80,16))];
			faces.push_back(face);
			PIDS.push_back(pid);
			iface++;
			getline(is,sread);
			ssub = trim(sread.substr(0,8));
		}

	}
	is.close();
	return;

}

string
NastranInterface::trim(string in){

	stringstream out;
	out << in;
	out >> in;
	return in;
}


string
NastranInterface::convertVertex(string in){
	int pos = in.find_last_of("-");
	if (pos<in.size() && pos > 0){
		in.insert(pos, "E");
	}
	else{
		pos = in.find_last_of("+");
		if (pos<in.size() && pos > 0){
				in.insert(pos, "E");
		}
	}
	return in;
}



