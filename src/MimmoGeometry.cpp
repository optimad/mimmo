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
#include <iostream>

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
	clear();
};

/*!Copy constructor of MimmoGeometry.Soft Copy of MimmoObject;
 */
MimmoGeometry::MimmoGeometry(const MimmoGeometry & other){
	*this = other;
};

/*!
 * Assignement operator of MimmoGeometry. Soft copy of MimmoObject
 */
MimmoGeometry & MimmoGeometry::operator=(const MimmoGeometry & other){
	clear();
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
	m_codex = other.m_codex;
	m_buildBvTree = other.m_buildBvTree;
	m_buildKdTree = other.m_buildKdTree;
	
	if(other.m_isInternal){
		m_geometry = other.m_intgeo.get();
	}	
	m_isInternal = false;
	return *this;
};

/*!
 * Return a const pointer to your current class.
 */
const MimmoGeometry * MimmoGeometry::getCopy(){
	return this;
}

/*!
 * Get current geometry pointer. Reimplementation of BaseManipulation::getGeometry
 */
MimmoObject * MimmoGeometry::getGeometry(){
	if(m_isInternal)	return m_intgeo.get();
	else				return m_geometry;
}

/*!
 * Get current geometry pointer. Reimplementation of BaseManipulation::getGeometry,
 * const overloading
 */
const MimmoObject * MimmoGeometry::getGeometry() const{
	if(m_isInternal)	return m_intgeo.get();
	else				return m_geometry;
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
	m_isInternal  = true;
	m_codex = true;
	m_buildBvTree = false;
	m_buildKdTree = false;
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

/*!
 * set codex ASCII false, BINARY true for writing sessions ONLY.
 * Default is Binary/Appended. Pay attention, binary writing is effective
 * only those file formats which support it.(ex STL, STVTU, SQVTU, VTVTU, VHVTU)
 */
void MimmoGeometry::setCodex(bool binary){
	m_codex = binary;
}

/*!Sets your current class as a "soft" copy of the argument.
 * Soft copy means that only your current geometric object MimmoObject is
 * copied only through its pointer and stored in the internal member m_geometry
 * Other members are exactly copied
 * \param[in] other pointer to MimmoGeometry class.
 */
void
MimmoGeometry::setSOFTCopy(const MimmoGeometry * other){
	clear();	
	*this = *other;
}

/*!Sets your current class as a "HARD" copy of the argument.
 * Hard copy means that only your current geometric object MimmoObject is
 * copied as a stand alone internal member and stored in the unique pointer member m_intgeo. 
 * Other members are exactly copied
 * \param[in] other pointer to MimmoGeometry class.
 */
void
MimmoGeometry::setHARDCopy(const MimmoGeometry * other){

	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(other));
	
	std::unique_ptr<MimmoObject> dum (new MimmoObject());
	dum->setHARDCopy(other->getGeometry());

	m_geometry = NULL;
	m_isInternal = true;
	m_intgeo = std::move(dum);
	
	m_rtype = other->m_rtype;
	m_wtype = other->m_wtype;
	m_read = other->m_read;
	m_rfilename = other->m_rfilename;
	m_write = other->m_write;
	m_wfilename = other->m_wfilename;
	m_rdir = other->m_rdir;
	m_wdir = other->m_wdir;
	m_wformat = other->m_wformat;
	m_codex = other->m_codex;
	m_buildBvTree = other->m_buildBvTree;
	m_buildKdTree = other->m_buildKdTree;
}

/*!
 * Set geometry from an external MimmoObject source, softly linked. 
 * Reimplementation of BaseManipulation::setGeometry
 */
void
MimmoGeometry::setGeometry(MimmoObject * external){
	
	if(getGeometry() == external) return;
	
	m_intgeo.reset(nullptr);
	m_geometry = external;
	m_isInternal = false;
};

/*!
 * Force your class to allocate an internal MimmoObject of type 1-Superficial mesh
 * 2-Volume Mesh. Other internal object allocated or externally linked geometries
 * will be destroyed/unlinked.
 * \param[in] type 1-Surface MimmoObject, 2-Volume MimmoObject. Default is 1, no other type are supported
 */
void
MimmoGeometry::setGeometry(int type){
	int type_ = std::max(std::min(2,type),1);
	m_geometry = NULL;
	m_intgeo.reset(nullptr);
	std::unique_ptr<MimmoObject> dum(new MimmoObject(type_));
	m_intgeo = std::move(dum);
	m_isInternal = true;
};

/*!
 * Return a pointer to the Vertex structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr
 */
bitpit::PiercedVector<bitpit::Vertex> * MimmoGeometry::getVertices(){
	if(isEmpty())	return NULL;
	return &(getGeometry()->getVertices());
};

/*!
 * Return a pointer to the Cell structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr 
 */
bitpit::PiercedVector<bitpit::Cell> * MimmoGeometry::getCells(){
	if(isEmpty())	return NULL;
	return	&(getGeometry()->getCells());
	
};

/*!
 * Wrapping to set vertices to your internal MimmoObject structure ONLY. Vertices in the internal
 * structure will be erased and substituted by the new ones. If an internal object is not 
 * allocated and the class is pointing to an external MimmoObject does nothing and return. Return without 
 * doing anything even in case of argument pointing to nullptr; 
 * \param[in] vertices pointer to a vertex PiercedVector structure 
 */
void
MimmoGeometry::setVertices(bitpit::PiercedVector<bitpit::Vertex> * vertices){
	if(m_intgeo.get() == NULL || vertices == NULL) return ;
	m_intgeo->setVertices(*vertices);
};

/*!
 * Wrapping to set cells to your internal MimmoObject structure ONLY. Cells in the internal
 * structure will be erased and substituted by the new ones. If an internal object is not 
 * allocated and the class is pointing to an external MimmoObject does nothing and return. Return without 
 * doing anything even in case of argument pointing to nullptr; 
 * \param[in] vertices pointer to a cell PiercedVector structure 
 */
void
MimmoGeometry::setCells(bitpit::PiercedVector<bitpit::Cell> * cells){
	if(m_intgeo.get() == NULL || cells == NULL) return;
	m_intgeo->setCells(*cells);
};

/*!It sets the PIDs of all the cells of the geometry Patch.
 * \param[in] pids PIDs of the cells of geometry mesh.
 */
void
MimmoGeometry::setPID(shivector1D pids){
	getGeometry()->setPID(pids);
};

/*!It sets if the BvTree of the patch has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
MimmoGeometry::setBuildBvTree(bool build){
	m_buildBvTree = build;
}

/*!It sets if the KdTree of the patch has to be built during execution.
* \param[in] build If true the KdTree is built in execution and stored in
* the related MimmoObject member.
*/
void
MimmoGeometry::setBuildKdTree(bool build){
	m_buildKdTree = build;
}

/*!
 * Check if geometry is not linked or not locally instantiated in your class.
 * True - no geometry present, False otherwise.
 */
bool MimmoGeometry::isEmpty(){
	return (m_geometry == NULL && (m_intgeo.get() == NULL));
}

/*!
 * Check if geometry is internally instantiated (true) or externally linked(false).
 * Return false if no geometry is checked. Please verify it with isEmpty method first
 */
bool MimmoGeometry::isInternal(){
	return (m_isInternal);
}

/*!
 * Clear all stuffs in your class
 */
void
MimmoGeometry::clear(){
	setDefaults();
	m_intgeo.reset(nullptr);
	BaseManipulation::clear();
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
			dynamic_cast<SurfUnstructured*>(getGeometry()->getPatch())->exportSTL(name, m_codex);
			return true;
			}
		break;
		
		case FileType::STVTU :
			//Export Triangulation Surface VTU
			{
				dvecarr3E	points = getGeometry()->getVertexCoords();
				ivector2D	connectivity = getGeometry()->getCompactConnectivity();
				bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TRIANGLE, points , connectivity );
				if(!m_codex)	vtk.setCodex(bitpit::VTKFormat::ASCII);
				else			vtk.setCodex(bitpit::VTKFormat::APPENDED);
				vtk.write() ;
				return true;
			}
		break;
		
		case FileType::SQVTU :
		//Export Quadrilateral Surface VTU
		{
			dvecarr3E	points = getGeometry()->getVertexCoords();
			ivector2D	connectivity = getGeometry()->getCompactConnectivity();
			bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::QUAD, points , connectivity );
			if(!m_codex)	vtk.setCodex(bitpit::VTKFormat::ASCII);
			else			vtk.setCodex(bitpit::VTKFormat::APPENDED);
			vtk.write() ;
			return true;
		}
		break;
		
		case FileType::VTVTU :
			//Export Tetra Volume VTU
			{
				dvecarr3E	points = getGeometry()->getVertexCoords();
				ivector2D	connectivity = getGeometry()->getCompactConnectivity();
				bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::TETRA, points , connectivity );
				if(!m_codex)	vtk.setCodex(bitpit::VTKFormat::ASCII);
				else			vtk.setCodex(bitpit::VTKFormat::APPENDED);		
				vtk.write() ;
				return true;
			}
		break;
		
		case FileType::VHVTU :
		//Export Hexa Volume VTU
		{
			dvecarr3E	points = getGeometry()->getVertexCoords();
			ivector2D	connectivity = getGeometry()->getCompactConnectivity();
			bitpit::VTKUnstructuredGrid  vtk(m_wdir, m_wfilename, bitpit::VTKElementType::HEXAHEDRON, points , connectivity );
			if(!m_codex)	vtk.setCodex(bitpit::VTKFormat::ASCII);
			else			vtk.setCodex(bitpit::VTKFormat::APPENDED);
			vtk.write() ;
			return true;
		}
		break;
	
		case FileType::NAS :
		//Export Nastran file
		{
			dvecarr3E	points = getGeometry()->getVertexCoords();
			ivector2D	connectivity = getGeometry()->getCompactConnectivity();
			NastranInterface nastran;
			nastran.setWFormat(m_wformat);
			shivector1D & pids = getGeometry()->getPid();
			if (pids.size() == connectivity.size()){
				nastran.write(m_wdir,m_wfilename,points,connectivity, &pids);
			}else{
				nastran.write(m_wdir,m_wfilename,points,connectivity);
			}
			return true;
		}
		break;
	
		case FileType::OFP : //Export ascii OpenFOAM point cloud
		{

			dvecarr3E points = getGeometry()->getVertexCoords();
			writeOFP(m_wdir, m_wfilename, points);
			return true;
		}
		break;
		default: //never been reached
			break;
	}
	
};

/*!It reads the mesh geometry from an input file and reverse it in the internal 
 * MimmoObject container. If an external container is linked skip reading and doing nothing.
 * \return False if file doesn't exists or not found geometry container address.
 */
bool
MimmoGeometry::read(){
	if(!m_isInternal) return false;
	
	switch(m_rtype){

	//Import STL
	case FileType::STL :
	{
		setGeometry(1);
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

		dynamic_cast<SurfUnstructured*>(getGeometry()->getPatch())->importSTL(name, binary);
		getGeometry()->cleanGeometry();
	}
	break;
	
	case FileType::STVTU :
		//Import Triangulation Surface VTU
	{
		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;
		
		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TRIANGLE, Ipoints, Iconnectivity );
		vtk.read() ;
		
		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
		
		setGeometry(1);
		
		int sizeV, sizeC;
		sizeV = Ipoints.size();
		sizeC = Iconnectivity.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		m_intgeo->getPatch()->reserveCells(sizeC);
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		for(auto & cc : Iconnectivity)	{
			livector1D temp(cc.size());
			int counter = 0;
			for(auto && val : cc){
				temp[counter] = val;
				++counter;
			}	
			m_intgeo->addConnectedCell(temp, eltype);
			
		}	
				
		m_intgeo->cleanGeometry();
	}
	break;
	
	case FileType::SQVTU :
		//Import Quadrilateral Surface VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;
	
		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::QUAD, Ipoints, Iconnectivity );
		vtk.read() ;

		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::QUAD;
		
		setGeometry(1);
		
		int sizeV, sizeC;
		sizeV = Ipoints.size();
		sizeC = Iconnectivity.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		m_intgeo->getPatch()->reserveCells(sizeC);
		
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		for(auto & cc : Iconnectivity)	{
			livector1D temp(cc.size());
			int counter = 0;
			for(auto && val : cc){
				temp[counter] = val;
				++counter;
			}	
			m_intgeo->addConnectedCell(temp, eltype);
			
		}	
		
		
		m_intgeo->cleanGeometry();
	}
	break;
	
	case FileType::VTVTU :
		//Import Tetra Volume VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TETRA, Ipoints, Iconnectivity );
		vtk.read() ;

		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TETRA;
		
		setGeometry(2);
		
		int sizeV, sizeC;
		sizeV = Ipoints.size();
		sizeC = Iconnectivity.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		m_intgeo->getPatch()->reserveCells(sizeC);
		
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		for(auto & cc : Iconnectivity)	{
			livector1D temp(cc.size());
			int counter = 0;
			for(auto && val : cc){
				temp[counter] = val;
				++counter;
			}	
			m_intgeo->addConnectedCell(temp, eltype);
			
		}	
		
		
		m_intgeo->cleanGeometry();
	}
	break;
	
	case FileType::VHVTU :
		//Import Hexa Volume VTU
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".vtu");
		bool check = infile.good();
		if (!check) return false;

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::HEXAHEDRON, Ipoints, Iconnectivity );
		vtk.read() ;

		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::HEXAHEDRON;
		
		setGeometry(2);
		
		int sizeV, sizeC;
		sizeV = Ipoints.size();
		sizeC = Iconnectivity.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		m_intgeo->getPatch()->reserveCells(sizeC);
		
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		for(auto & cc : Iconnectivity)	{
			livector1D temp(cc.size());
			int counter = 0;
			for(auto && val : cc){
				temp[counter] = val;
				++counter;
			}	
			m_intgeo->addConnectedCell(temp, eltype);
			
		}	
		
		m_intgeo->cleanGeometry();
	}
	break;
	
	case FileType::NAS :
		//Import Surface NAS
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename+".nas");
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		dvecarr3E	Ipoints ;
		ivector2D	Iconnectivity ;

		NastranInterface nastran;
		nastran.setWFormat(m_wformat);

		shivector1D pids;
		nastran.read(m_rdir, m_rfilename, Ipoints, Iconnectivity, pids );

		bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
		
		setGeometry(1);
		
		int sizeV, sizeC;
		sizeV = Ipoints.size();
		sizeC = Iconnectivity.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		m_intgeo->getPatch()->reserveCells(sizeC);
		
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		for(auto & cc : Iconnectivity)	{
			livector1D temp(cc.size());
			int counter = 0;
			for(auto && val : cc){
				temp[counter] = val;
				++counter;
			}	
			m_intgeo->addConnectedCell(temp, eltype);
			
		}	
		
		m_intgeo->cleanGeometry();
	}
	
	break;
	
	//Import ascii OpenFOAM point cloud
	case FileType::OFP :
	{

		std::ifstream infile(m_rdir+"/"+m_rfilename);
		bool check = infile.good();
		if (!check) return false;
		infile.close();

		dvecarr3E	Ipoints;
		readOFP(m_rdir, m_rfilename, Ipoints);

		setGeometry(3);

		int sizeV = Ipoints.size();
		m_intgeo->getPatch()->reserveVertices(sizeV);
		
		for(auto & vv : Ipoints)		m_intgeo->addVertex(vv);
		
	}
	break;
	
	default: //never been reached
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
	bool check = true;
	if (m_read) check = read();
	if (!check){
		std::cout << "MiMMO : ERROR : file not found : "<< m_rfilename << std::endl;
		std::cout << " " << std::endl;
		exit(10);
	}
	check = true;
	if (m_write) check = write();
	if (!check){
		std::cout << "MiMMO : ERROR : write not done : geometry not linked " << std::endl;
		std::cout << " " << std::endl;
		exit(11);
	}
	if (m_buildBvTree) getGeometry()->buildBvTree();
	if (m_buildKdTree) getGeometry()->buildKdTree();
	return;
}


//===============================//
//====== OFOAM INTERFACE ========//
//===============================//

/*!
 *	Read openFoam format geometry file and absorb it as a point cloud ONLY.
 *\param[in]	inputDir	folder of file
 *\param[in]	surfaceName	name of file
 *\param[out]	points		list of points in the cloud  
 * 
 */
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
/*!
 *	Write geometry file in openFoam format as a point cloud ONLY.
 *\param[in]	inputDir	folder of file
 *\param[in]	surfaceName	name of file
 *\param[out]	points		list of points in the cloud  
 * 
 */
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



