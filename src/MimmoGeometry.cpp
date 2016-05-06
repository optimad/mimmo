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
	setDefaults();
}

/*!Default destructor of MimmoGeometry.
 */
MimmoGeometry::~MimmoGeometry(){
	MimmoObject::clear();
};

/*!Copy constructor of MimmoGeometry.
 */
MimmoGeometry::MimmoGeometry(const MimmoGeometry & other){
	*this = other;
};

/*!Assignement operator of MimmoGeometry.
 */
MimmoGeometry & MimmoGeometry::operator=(const MimmoGeometry & other){
	*(static_cast<MimmoObject * >(this)) = *(static_cast<const MimmoObject * >(&other)); //soft copy
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_rtype = other.m_rtype;
	m_wtype = other.m_wtype;
	m_read = other.m_read;
	m_rfilename = other.m_rfilename;
	m_write = other.m_write;
	m_wfilename = other.m_wfilename;
	m_rdir = other.m_rdir;
	m_wdir = other.m_wdir;
	m_wformat = other.m_wformat;
	return *this;
};

/*!
 * Return a const pointer to your current class.
 * Reimplementation of MimmoObject::getCopy;
 */
const MimmoGeometry * MimmoGeometry::getCopy(){
	return this;
}

/*!
 * Return a pointer as MimmoObject to your current class.
 * Static casting is enforced. Reimplementation of BaseManipulation::getGeometry();
 */
MimmoObject * MimmoGeometry::getGeometry(){
	return	static_cast<MimmoObject *>(this);
}
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

/*!Sets your current class as a "soft" copy of the argument.
 * Reimplementation of MimmoObject::setSOFTCopy(), get its documentation
 * for further information.
 * \param[in] other pointer to MimmoGeometry class.
 */
void
MimmoGeometry::setSOFTCopy(const MimmoGeometry * other){
	clear();	
	*this = *other;
}

/*!Sets your current class as a "HARD" copy of the argument.
 * Reimplementation of MimmoObject::setSOFTCopy(), get its documentation
 * for further information.
 * \param[in] other pointer to MimmoGeometry class.
 */
void
MimmoGeometry::setHARDCopy(const MimmoGeometry * other){
	clear();	
	MimmoObject::setHARDCopy(static_cast<const MimmoObject * >(other));
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(other));
	m_rtype = other->m_rtype;
	m_wtype = other->m_wtype;
	m_read = other->m_read;
	m_rfilename = other->m_rfilename;
	m_write = other->m_write;
	m_wfilename = other->m_wfilename;
	m_rdir = other->m_rdir;
	m_wdir = other->m_wdir;
	m_wformat = other->m_wformat;
}

/*!
 * Set geometry. Do Nothing. Reimplementation of BaseManipulation::setGeometry
 */
void
MimmoGeometry::setGeometry(){};

/*!It sets the PIDs of all the cells of the geometry Patch.
 * \param[in] pids PIDs of the cells of geometry mesh.
 */
void
MimmoGeometry::setPID(shivector1D pids){
	if((int)pids.size() != getNCells())	return;
	m_pids.clear();
	m_pids = pids;
	
	//find different type of pids.
	std::unordered_set<short> map;
	map.insert(pids.begin(),pids.end());
	m_pidsType.clear();
	m_pidsType.insert(m_pidsType.end(), map.begin(), map.end());
};

/*!
 * Clear all stuffs in your class
 */

void
MimmoGeometry::clear(){
	MimmoObject::clear();
	BaseManipulation::clear();
	setDefaults();
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
	if (isEmpty()) return false;

	switch(m_wtype){
	case FileType::STL :
		//Export STL
	{
		string name = (m_wdir+"/"+m_wfilename+".stl");
		//forced binary write
		bool binary = true;
		dynamic_cast<SurfUnstructured*>(getPatch())->exportSTL(name, binary);
		return true;
	}
	break;
	case FileType::STVTU :
		//Export Triangulation Surface VTU
	{
		dvecarr3E	points = getVertex();
		ivector2D	connectivity = getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<3; j++){
				connectivity[i][j] = getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TRIANGLE, points , connectivity );
		vtk.write() ;
		return true;
	}
	break;
	case FileType::SQVTU :
		//Export Quadrilateral Surface VTU
	{
		dvecarr3E	points = getVertex();
		ivector2D	connectivity = getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<4; j++){
				connectivity[i][j] = getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::QUAD, points , connectivity );
		vtk.write() ;
		return true;
	}
	break;
	case FileType::VTVTU :
		//Export Tetra Volume VTU
	{
		dvecarr3E	points = getVertex();
		ivector2D	connectivity = getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<4; j++){
				connectivity[i][j] =getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TETRA, points , connectivity );
		vtk.write() ;
		return true;
	}
	case FileType::VHVTU :
		//Export Hexa Volume VTU
	{
		dvecarr3E	points = getVertex();
		ivector2D	connectivity = getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<8; j++){
				connectivity[i][j] = getMapDataInv(connectivity[i][j]);
			}
		}
		bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::HEXAHEDRON, points , connectivity );
		vtk.write() ;
		return true;
	}
	break;
	case FileType::NAS :
		//Export Nastran file
	{
		dvecarr3E	points = getVertex();
		ivector2D	connectivity = getConnectivity();
		for (int i=0; i<connectivity.size(); i++){
			for (int j=0; j<3; j++){
				connectivity[i][j] = getMapDataInv(connectivity[i][j]);
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
		return true;
	}
	break;
	//Export ascii OpenFOAM point cloud
	case FileType::OFP :
	{

		dvecarr3E points = getVertex();
		writeOFP(m_wdir, m_wfilename, points);

		return true;
	}
	break;
	}
};

/*!It reads the mesh geometry from an input file.
 * \return False if file doesn't exists.
 */
bool
MimmoGeometry::read(){

	switch(m_rtype){

	//Import STL
	case FileType::STL :
	{
		MimmoObject::clear();
		m_type = 1;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured * >(m_patch)->setExpert(true);
		string name;

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

		dynamic_cast<SurfUnstructured*>(getPatch())->importSTL(name, binary);
		cleanGeometry();
	}
	break;
	case FileType::STVTU :
		//Import Triangulation Surface VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject::clear();
		m_type = 1;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured * >(m_patch)->setExpert(true);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;
		
		//TODO Import with generic element type ?!?!

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TRIANGLE, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
		}
		
		// convert Iconnectivity to livector1D
		livector2D  Iconn;
		{
			Iconn.resize(Iconnectivity.size());
			livector1D::iterator it;
			int counter=0;
			for(auto && val : Iconnectivity){
				it=Iconn[counter].begin();
				for(auto && element : val){
					*it = long(element);
					++it;
				}
				++counter;
			}
			Iconnectivity.clear();
		}
		
		
		setConnectivity(&Iconn);
		cleanGeometry();
	}
	break;
	case FileType::SQVTU :
		//Import Quadrilateral Surface VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject::clear();
		m_type = 1;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured * >(m_patch)->setExpert(true);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::QUAD, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
		}
		
		// convert Iconnectivity to livector1D
		livector2D  Iconn;
		{
			Iconn.resize(Iconnectivity.size());
			livector1D::iterator it;
			int counter=0;
			for(auto && val : Iconnectivity){
				it=Iconn[counter].begin();
				for(auto && element : val){
					*it = long(element);
					++it;
				}
				++counter;
			}
			Iconnectivity.clear();
		}
		
		
		setConnectivity(&Iconn);
		
		cleanGeometry();
	}
	break;
	case FileType::VTVTU :
		//Import Tetra Volume VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject::clear();
		m_type = 2;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new VolUnstructured(id,3);
		dynamic_cast<VolUnstructured * >(m_patch)->setExpert(true);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		//TODO Import with generic element type ?!?!

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TETRA, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
		}
		
		// convert Iconnectivity to livector1D
		livector2D  Iconn;
		{
			Iconn.resize(Iconnectivity.size());
			livector1D::iterator it;
			int counter=0;
			for(auto && val : Iconnectivity){
				it=Iconn[counter].begin();
				for(auto && element : val){
					*it = long(element);
					++it;
				}
				++counter;
			}
			Iconnectivity.clear();
		}
		
		
		setConnectivity(&Iconn);
		
		cleanGeometry();
	}
	break;
	case FileType::VHVTU :
		//Import Hexa Volume VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		MimmoObject::clear();
		m_type = 2;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new VolUnstructured(id,3);
		dynamic_cast<VolUnstructured * >(m_patch)->setExpert(true);
		
		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::HEXAHEDRON, Ipoints, Iconnectivity );
		vtk.read() ;

		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
		}
		
		// convert Iconnectivity to livector1D
		livector2D  Iconn;
		{
			Iconn.resize(Iconnectivity.size());
			livector1D::iterator it;
			int counter=0;
			for(auto && val : Iconnectivity){
				it=Iconn[counter].begin();
				for(auto && element : val){
					*it = long(element);
					++it;
				}
				++counter;
			}
			Iconnectivity.clear();
		}
		
		
		setConnectivity(&Iconn);
		
		cleanGeometry();
	}
	break;
	case FileType::NAS :
		//Import Surface NAS
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".nas");
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		MimmoObject::clear();
		m_type = 1;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured * >(m_patch)->setExpert(true);

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		NastranInterface nastran;
		nastran.setWFormat(m_wformat);

		shivector1D pids;
		nastran.read(m_rdir, m_rfilename, Ipoints, Iconnectivity, pids );
		setPID(pids);
		
		int	np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
		}
		
		// convert Iconnectivity to livector1D
		livector2D  Iconn;
		{
			Iconn.resize(Iconnectivity.size());
			livector1D::iterator it;
			int counter=0;
			for(auto && val : Iconnectivity){
				it=Iconn[counter].begin();
				for(auto && element : val){
					*it = long(element);
					++it;
				}
				++counter;
			}
			Iconnectivity.clear();
		}
		
		
		setConnectivity(&Iconn);
		
		cleanGeometry();
	}
	break;
	//Import ascii OpenFOAM point cloud
	case FileType::OFP :
	{
		//It uses surface type of mimmo object for the moment
		MimmoObject::clear();
		m_type = 1;
		m_internalPatch = true;
		const int  id =0;
		m_patch = new SurfUnstructured(id);
		dynamic_cast<SurfUnstructured * >(m_patch)->setExpert(true);
		

		std::ifstream infile(m_rdir+"/"+m_rfilename);
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		dvecarr3E	Ipoints;
		readOFP(m_rdir, m_rfilename, Ipoints);

		int np = Ipoints.size();
		for (long ip=0; ip<np; ip++){
			setVertex(Ipoints[ip]);
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
//====== OFOAM INTERFACE ========//
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
			os << setprecision(16) << points[i][j] << separator;
		}
		os << setprecision(16) << points[i][2] << parr << nl;
	}
	os << parr << nl;
	os << nl;
	os << nl;
	hline = "// ************************************************************************* //";
	os << hline << nl;

	os.close();
	return;

}

/*!
 * Set proper member of the class to defaults
 */
void
MimmoGeometry::setDefaults(){
	m_rtype		= STL;
	m_read		= false;
	m_rfilename	= "mimmoGeometry";
	m_wtype		= STL;
	m_write		= false;
	m_wfilename	= "mimmoGeometry";
	m_rdir		= "./";
	m_wdir		= "./";
	m_wformat	= Short;
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

void NastranInterface::writeFace(string faceType, ivector1D& facePts, int& nFace, ofstream& os,int PID){
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

void NastranInterface::writeGeometry(dvecarr3E& points, ivector2D& faces, ofstream& os, shivector1D* PIDS){
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

//TODO we have PIDs, we need to write them here in the footer also. CHECK IT OUT. You can access only
// pidType through the MimmoGeometry member m_pidsType

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

void NastranInterface::write(string& outputDir, string& surfaceName, dvecarr3E& points, ivector2D& faces, shivector1D* PIDS){

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

void NastranInterface::read(string& inputDir, string& surfaceName, dvecarr3E& points, ivector2D& faces, shivector1D& PIDS){

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
			PIDS.push_back((short)pid);
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



