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

#include "MultipleMimmoGeometries.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!Default constructor of MultipleMimmoGeometries.
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(int topo, int formattype){
	initializeClass(topo, formattype);
}

/*!Custom constructor of MultipleMimmoGeometries.
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(int topo, FileType formattype){
	initializeClass(topo, formattype._to_integer());
}

/*!Default destructor of MultipleMimmoGeometries.
 */
MultipleMimmoGeometries::~MultipleMimmoGeometries(){
	clear();
};

/*!Copy constructor of MultipleMimmoGeometries.Soft Copy of MimmoObject;
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(const MultipleMimmoGeometries & other){
	*this = other;
};	

/*!
 * Assignement operator of MultipleMimmoGeometries. Soft copy of MimmoObject
 */
MultipleMimmoGeometries & MultipleMimmoGeometries::operator=(const MultipleMimmoGeometries & other){
	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(&other));
	m_rinfo = other.m_rinfo;
	m_winfo = other.m_winfo;
	m_read = other.m_read;
	m_write = other.m_write;
	m_wformat = other.m_wformat;
	m_codex = other.m_codex;
	m_mapMGeo = other.m_mapMGeo;
	m_buildBvTree = other.m_buildBvTree;
	m_buildKdTree = other.m_buildKdTree;
	m_topo = other.m_topo;
	m_tagtype = other.m_tagtype;
	
	if(other.m_isInternal){
		m_geometry = other.m_intgeo.get();
	}	
	m_isInternal = false;
	return *this;
};

void
MultipleMimmoGeometries::buildPorts(){
	bool built = true;
	built = (built && createPortIn<MimmoObject*, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<std::unordered_map<long,short>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setDivisionMap, PortType::M_MAPGEOM, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGSHORT));
	built = (built && createPortIn<std::vector<FileInfoData>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setReadListAbsolutePathFiles, PortType::M_FINFO, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	built = (built && createPortIn<std::vector<FileInfoData>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setWriteListAbsolutePathFiles, PortType::M_FINFO2, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	
	built = (built && createPortOut<MimmoObject*, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getGeometry, PortType::M_GEOM, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortOut<std::unordered_map<long,short>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getDivisionMap, PortType::M_MAPGEOM, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::LONGSHORT));
	built = (built && createPortOut<std::vector<FileInfoData>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getReadListAbsolutePathFiles, PortType::M_FINFO, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	built = (built && createPortOut<std::vector<FileInfoData>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getWriteListAbsolutePathFiles, PortType::M_FINFO2, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	
	m_arePortsBuilt = built;
}


/*!
 * Return current format file allowed in your class. Refers to enum FileType
 */
int
MultipleMimmoGeometries::getFormatTypeAllowed(){
	return m_tagtype;
}

/*!
 * Return a const pointer to your current class.
 */
const MultipleMimmoGeometries * MultipleMimmoGeometries::getCopy(){
	return this;
}

/*!
 * Get current geometry pointer. Reimplementation of BaseManipulation::getGeometry
 */
MimmoObject * MultipleMimmoGeometries::getGeometry(){
	if(m_isInternal)	return m_intgeo.get();
	else				return m_geometry;
}

/*!
 * Get current geometry pointer. Reimplementation of BaseManipulation::getGeometry,
 * const overloading
 */
const MimmoObject * MultipleMimmoGeometries::getGeometry() const{
	if(m_isInternal)	return m_intgeo.get();
	else				return m_geometry;
}

/*!
 * Get current division map associated to your MimmoObject geometry.
 * (division map is a map which associates for each cell of MimmoObject a tag id
 * of the origin/destination file).
 */
std::unordered_map<long short>	MultipleMimmoGeometries::getDivisionMap(){
	return m_mapMGeo;
}

/*!
 * Get how many sub-files are available for writing, once division map is set. 
 * Meant only for WRITING mode of the class.
 */
int
MultipleMimmoGeometries::getHowManySubFiles(){
	std::unordered_set subtypes;
	for(auto & element: m_mapMGeo){
		subtypes.insert(element.second);
	}
	
	return subtypes.size();
}

/*!
 * Get list of external filenames the class must READ from to get the final MimmoObject.
 */
std::vector<FileInfoData>	MultipleMimmoGeometries::getReadListAbsolutePathFiles(){
	return m_rinfo;
}

/*!
 * Get list of external filenames the class must WRITE to the final MimmoObject.
 */
std::vector<FileInfoData>	MultipleMimmoGeometries::getReadListAbsolutePathFiles(){
	return m_winfo;
}

/*!
 * Get list of PID type actually stored in your MimmoObjects. While reading, PIDs from different source files
 * are all stored in the unique MimmoObject data structure, owned by the class. PIDs are automatically
 * reassigned in case of conflict between PIDs of same value coming from different files. 
 */
std::unordered_set<short>& MultipleMimmoGeometries::getPIDList(){
	getGeometry()->getPIDTypeList();
}

/*!
 * Add an external file path to the list of geometries to read from.
 * Available format type are related to the type of geometry topology 
 * set for your class and the format fixed for your class also. See FileType enum.
 * \param[in] dir  string, absolute path directory
 * \param[in] name string, name of your file, without tag extension
 */
void
MultipleMimmoGeometries::setAddReadAbsolutePathFile(std::string dir, std::string name){
	FileDataInfo temp;
	temp.ftype	= m_tagtype;
	temp.fdir	= dir;
	temp.fname	= name;

	m_rinfo.push_back(temp);
};

/*!
 * Add an external file path to the list of geometries to write to.
 * Available format type are related to the type of geometry topology 
 * set for your class and the format fixed for your class also. See FileType enum.
 * \param[in] dir  string, absolute path directory
 * \param[in] name string, name of your file, without tag extension
 */
void
MultipleMimmoGeometries::setAddWriteAbsolutePathFile(std::string dir, std::string name){
	FileDataInfo temp;
	temp.ftype	= m_tagtype;
	temp.fdir	= dir;
	temp.fname	= name;
	
	m_winfo.push_back(temp);
	
};

/*!
 * Set the whole list of external files data which the geometry is read from.
 */
void
MultipleMimmoGeometries::setReadListAbsolutePathFiles(std::vector<FileInfoData> data){
	m_rinfo.clear();
	m_rinfo.resize(data.size);
	
	int counter = 0;
	for(auto & ele : data){
		if(ele.ftype == m_tagtype){
			m_rinfo[counter] = ele;
			++counter;
		}
	}
	m_rinfo.resize(counter);
};

/*!
 * Set the whole list of external files data which the geometry is written to.
 */
void
MultipleMimmoGeometries::setWriteListAbsolutePathFiles(std::vector<FileInfoData> data){
	m_winfo.clear();
	m_winfo.resize(data.size);
	
	int counter = 0;
	for(auto & ele : data){
		if(ele.ftype == m_tagtype){
			m_winfo[counter] = ele;
			++counter;
		}
	}
	m_winfo.resize(counter);
	
};

/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
MultipleMimmoGeometries::setRead(bool read){
	m_read = read;
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MultipleMimmoGeometries::setWrite(bool write){
	m_write = write;
}

/*!
 * set codex ASCII false, BINARY true for writing sessions ONLY.
 * Default is Binary/Appended. Pay attention, binary writing is effective
 * only those file formats which support it.(ex STL, STVTU, SQVTU, VTVTU, VHVTU)
 */
void MultipleMimmoGeometries::setCodex(bool binary){
	m_codex = binary;
}

/*!Sets your current class as a "soft" copy of the argument.
 * Soft copy means that only your current geometric object MimmoObject is
 * copied only through its pointer and stored in the internal member m_geometry
 * Other members are exactly copied
 * \param[in] other pointer to MultipleMimmoGeometries class.
 */
void
MultipleMimmoGeometries::setSOFTCopy(const MultipleMimmoGeometries * other){
	clear();	
	*this = *other;
}

/*!Sets your current class as a "HARD" copy of the argument.
 * Hard copy means that only your current geometric object MimmoObject is
 * copied as a stand alone internal member and stored in the unique pointer member m_intgeo. 
 * Other members are exactly copied
 * \param[in] other pointer to MultipleMimmoGeometries class.
 */
void
MultipleMimmoGeometries::setHARDCopy(const MultipleMimmoGeometries * other){

	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(other));
	
	std::unique_ptr<MimmoObject> dum (new MimmoObject());
	dum->setHARDCopy(other->getGeometry());

	m_geometry = NULL;
	m_isInternal = true;
	m_intgeo = std::move(dum);
	
	m_rinfo = other->m_rinfo;
	m_winfo = other->m_winfo;
	m_read = other->m_read;
	m_write = other->m_write;
	m_wformat = other->m_wformat;
	m_codex = other->m_codex;
	
	m_topo = other.m_topo;
	m_tagtype = other.m_tagtype;
	
	m_mapMGeo = other->m_mapMGeo;
	
	m_buildBvTree = other->m_buildBvTree;
	m_buildKdTree = other->m_buildKdTree;
}

/*!
 * Set geometry from an external MimmoObject source, softly linked. 
 * Reimplementation of BaseManipulation::setGeometry
 */
void
MultipleMimmoGeometries::setGeometry(MimmoObject * external){
	
	if(getGeometry() == external) return;
	
	m_intgeo.reset(nullptr);
	m_geometry = external;
	m_isInternal = false;
};

/*!
 * Force your class to allocate an internal MimmoObject of m_topo type. 
 * Other internal object allocated or externally linked geometries
 * will be destroyed/unlinked.
 * 
 */
void
MultipleMimmoGeometries::setGeometry(){
	m_geometry = NULL;
	m_intgeo.reset(nullptr);
	std::unique_ptr<MimmoObject> dum(new MimmoObject(m_topo));
	m_intgeo = std::move(dum);
	m_isInternal = true;
};

/*!
 * Return a pointer to the Vertex structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr
 */
bitpit::PiercedVector<bitpit::Vertex> * MultipleMimmoGeometries::getVertices(){
	if(isEmpty())	return NULL;
	return &(getGeometry()->getVertices());
};

/*!
 * Return a pointer to the Cell structure of the MimmoObject geometry actually pointed or allocated by
 * the class. If no geometry is actually available return a nullptr 
 */
bitpit::PiercedVector<bitpit::Cell> * MultipleMimmoGeometries::getCells(){
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
MultipleMimmoGeometries::setVertices(bitpit::PiercedVector<bitpit::Vertex> * vertices){
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
MultipleMimmoGeometries::setCells(bitpit::PiercedVector<bitpit::Cell> * cells){
	if(m_intgeo.get() == NULL || cells == NULL) return;
	m_intgeo->setCells(*cells);
};

/*!It sets the PIDs of all the cells of the geometry Patch currently linked or internal to the class.
 * \param[in] pids PIDs of the cells of geometry mesh, in compact sequential order. If pids size does not match number of current cell does nothing
 */
void
MultipleMimmoGeometries::setPID(shivector1D pids){
	getGeometry()->setPID(pids);
};

/*!It sets the PIDs of part of/all the cells of the geometry Patch currently linked or internal to the class.
 * \param[in] pids PIDs map list w/ id of the cell as first value, and pid as second value.
 */
void
MultipleMimmoGeometries::setPID(std::unordered_map<long, short> pidsMap){
	getGeometry()->setPID(pidsMap);
};

/*!It sets the format to export .nas file.
 * \param[in] wform Format of .nas file (Short/Long).
 */
void
MultipleMimmoGeometries::setFormatNAS(WFORMAT wform){
	m_wformat = wform;
}

/*!It sets if the BvTree of the patch has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
MultipleMimmoGeometries::setBuildBvTree(bool build){
	m_buildBvTree = build;
}

/*!It sets if the KdTree of the patch has to be built during execution.
* \param[in] build If true the KdTree is built in execution and stored in
* the related MimmoObject member.
*/
void
MultipleMimmoGeometries::setBuildKdTree(bool build){
	m_buildKdTree = build;
}

/*!
 * Check if geometry is not linked or not locally instantiated in your class.
 * True - no geometry present, False otherwise.
 */
bool MultipleMimmoGeometries::isEmpty(){
	return (m_geometry == NULL && (m_intgeo.get() == NULL));
}

/*!
 * Check if geometry is internally instantiated (true) or externally linked(false).
 * Return false if no geometry is checked. Please verify it with isEmpty method first
 */
bool MultipleMimmoGeometries::isInternal(){
	return (m_isInternal);
}

/*!
 * Clear all stuffs in your class
 */
void
MultipleMimmoGeometries::clear(){
	clearReadPaths();
	clearWritePaths();
	setDefaults();
	m_intgeo.reset(nullptr);
	BaseManipulation::clear();
};

/*!
 * Clear all reading filenames added to the class 
 */
void
MultipleMimmoGeometries::clearReadPaths(){
	m_rinfo.clear();
}

/*!
 * Clear all writing filenames added to the class 
 */
void
MultipleMimmoGeometries::clearWritePaths(){
	m_winfo.clear();
}

/*!It writes the geometry on multiple output files.
 *\return False if geometry is not linked.
 */
bool
MultipleMimmoGeometries::write(){
	if (isEmpty()) return false;
	
	
	
	
};

/*!It reads the mesh geometry from an input file and reverse it in the internal 
 * MimmoObject container. If an external container is linked skip reading and doing nothing.
 * \return False if file doesn't exists or not found geometry container address.
 */
bool
MultipleMimmoGeometries::read(){
	if(!m_isInternal) return false;
	
	switch(FileType::_from_integral(m_rtype)){

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
		
		//count PID if multi-solid
		auto & map = getGeometry()->getPIDTypeList();
		for(auto & cell : getGeometry()->getCells() ){
			map.insert(cell.getPID());
		}

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
		
		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TRIANGLE);
		vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
		vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
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

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::QUAD);
		vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
		vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
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

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::TETRA);
		vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
		vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
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

		bitpit::VTKUnstructuredGrid  vtk(m_rdir, m_rfilename, bitpit::VTKElementType::HEXAHEDRON);
		vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, Ipoints) ;
		vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, Iconnectivity) ;
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
		m_intgeo->setPID(pids);
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
MultipleMimmoGeometries::execute(){
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
void MultipleMimmoGeometries::readOFP(string& inputDir, string& surfaceName, dvecarr3E& points){

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
void MultipleMimmoGeometries::writeOFP(string& outputDir, string& surfaceName, dvecarr3E& points){

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
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class reads the following parameters:
 * 
 * 1) ReadFlag - activate reading mode boolean
 * 2) ReadDir - reading directory path
 * 3) ReadFileType - file type identifier
 * 4) ReadFilename - name of file for reading
 * 5) WriteFlag - activate writing mode boolean
 * 6) WriteDir - writing directory path
 * 7) WriteFileType - file type identifier
 * 8) WriteFilename - name of file for writing
 * 9) Codex - boolean to write ascii/binary
 * 10) BvTree - evaluate bvTree true/false
 * 11) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void MultipleMimmoGeometries::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	std::string input; 

	if(slotXML.hasOption("ReadFlag")){
		input = slotXML.get("ReadFlag");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setRead(value);
	}; 
	
	if(slotXML.hasOption("ReadDir")){
		input = slotXML.get("ReadDir");
		input = bitpit::utils::trim(input);
		if(input.empty())	input = "./";
		setReadDir(input);
	}; 
	
	if(slotXML.hasOption("ReadFileType")){
		input = slotXML.get("ReadFileType");
		input = bitpit::utils::trim(input);
		if(!input.empty()){
			for(auto c: FileType::_values()){
				if(input == c._to_string()){
					setReadFileType(c);
				}
			}
		}else{
			setReadFileType(0);
		}
	}; 

	if(slotXML.hasOption("ReadFilename")){
		input = slotXML.get("ReadFilename");
		input = bitpit::utils::trim(input);
		if(input.empty())	input = "mimmoGeometry";
		setReadFilename(input);
	}; 
	
	if(slotXML.hasOption("WriteFlag")){
		input = slotXML.get("WriteFlag");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setWrite(value);
	}; 
	
	if(slotXML.hasOption("WriteDir")){
		input = slotXML.get("WriteDir");
		input = bitpit::utils::trim(input);
		if(input.empty())	input = "./";
	   setWriteDir(input);
	}; 
	
	if(slotXML.hasOption("WriteFileType")){
		input = slotXML.get("WriteFileType");
		input = bitpit::utils::trim(input);
		if(!input.empty()){
			for(auto c: FileType::_values()){
				if(input == c._to_string()){
					setWriteFileType(c);
				}
			}
		}else{
			setWriteFileType(0);
		}
	}; 
	
	if(slotXML.hasOption("WriteFilename")){
		input = slotXML.get("WriteFilename");
		input = bitpit::utils::trim(input);
		if(input.empty())	input = "mimmoGeometry";
	   setWriteFilename(input);
	}; 
	
	if(slotXML.hasOption("Codex")){
		input = slotXML.get("Codex");
		bool value = true;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setCodex(value);
	}; 
	
	
	if(slotXML.hasOption("BvTree")){
		input = slotXML.get("BvTree");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setBuildBvTree(value);
	}; 
	
	if(slotXML.hasOption("KdTree")){
		input = slotXML.get("KdTree");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setBuildKdTree(value);
	}; 

return;	
};

/*!
 * Write settings of the class to bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::flushSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class writes the following parameters(if different from default):
 * 
 * 1) ReadFlag - activate reading mode boolean
 * 2) ReadDir - reading directory path
 * 3) ReadFileType - file type identifier
 * 4) ReadFilename - name of file for reading
 * 5) WriteFlag - activate writing mode boolean
 * 6) WriteDir - writing directory path
 * 7) WriteFileType - file type identifier
 * 8) WriteFilename - name of file for writing
 * 9) Codex - boolean to write ascii/binary
 * 10) BvTree - evaluate bvTree true/false
 * 11) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot 
 */
void MultipleMimmoGeometries::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	
	std::string output;
	
	output = std::to_string(m_read);
	slotXML.set("ReadFlag", output);

	if(m_rdir != "./"){
		slotXML.set("ReadDir", m_rdir);
	}	
	
	if(m_rfilename != "mimmoGeometry"){
		slotXML.set("ReadFilename", m_rfilename);
	}	
	
	if(m_rtype != 0){
		std::string temp = (FileType::_from_integral(m_rtype))._to_string();
		slotXML.set("ReadFileType", temp);
	}
	
	output = std::to_string(m_write);
	slotXML.set("WriteFlag", output);
	
	if(m_wdir != "./"){
		slotXML.set("WriteDir", m_wdir);
	}	
	
	if(m_wfilename != "mimmoGeometry"){
		slotXML.set("WriteFilename", m_wfilename);
	}	
	
	if(m_wtype != FileType::STL){
		std::string temp = (FileType::_from_integral(m_wtype))._to_string();
		slotXML.set("WriteFileType", temp);
	}
	
	if(!m_codex){
		output = std::to_string(m_codex);
		slotXML.set("Codex", output);
	}
	
	if(m_buildBvTree){
		output = std::to_string(m_buildBvTree);
		slotXML.set("BvTree", output);
	}
	
	if(m_buildKdTree){
		output = std::to_string(m_buildKdTree);
		slotXML.set("KdTree", output);
	}
	
return;	
};	



/*!
 * Set proper member of the class to defaults
 */
void
MultipleMimmoGeometries::setDefaults(){
	m_rtype		= static_cast<int>(FileType::STL);
	m_read		= false;
	m_rfilename	= "mimmoGeometry";
	m_wtype		= static_cast<int>(FileType::STL);
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


/*!
 * Class initializer, meant to be used in cosntruction.
 */
void
MultipleMimmoGeometries::initializeClass(int topo, int formattype){
	
	m_name 		= "MiMMO.MultipleGeometries";
	m_topo     = std::min(1, topo);
	if(m_topo > 3)	m_topo = 1;
	
	//checking admissible format
	switch(m_topo){
		case 1:
			if(formattype !=0 || formattype !=1 || formattype !=2 ||formattype !=5){
				formattype = 0;
			}
			break;
		case 2:
			if(formattype !=3 || formattype !=4){
				formattype = 3;
			}
			break;
		default:
			if (formattype <0 || formattype > 6){
				formattype= 0;
			}
			break;
	}
	
	m_tagtype = formattype;
	
	setDefaults();
};

/*!
 * Extract MimmoObject cell list from m_mapMGeo, given the id part identifier
 */
livector1D
MultipleMimmoGeometries::cellExtractor(short id){
	livector1D result;
	
	result.resize(m_mapMGeo.size());
	int counter = 0;
	
	for(auto & elem : m_mapMGeo){
		if(elem.second==id){
			result[counter] = elem.first;
			++counter;
		}
	}
	
	result.resize(counter);
	return result;
};

/*!
 * Given a certain sub-group of cells of class linked/internal MimmoObject, fill a stand-alone
 * MimmoObject sub-extraction
 * \param[in] cellList list of cells id of original MimmoObject
 * \param[in,out] subG target extraction, allocated externally and here passed as a pointer.
 */
void
MultipleMimmoGeometries::fillSubStructure(livector1D & cellList, MimmoObject * subG){
	
	if (cellList.empty() || isEmpty() || subG == NULL) return;

	bitpit::PatchKernel * tri = getGeometry()->getPatch();
	bitpit::PiercedVector<bitpit::Vertex> & mapV = subG->getPatch()->getVertices();
	livector1D TT;
	
	long idV;
	short int PID;
	int sizeCC;
	bitpit::ElementInfo::Type eltype;
	
	for(auto & idCell: cellList){
		
		bitpit::Cell & cell = tri->getCell(idCell);
		eltype = cell.getType();
		sizeCC = cell.getVertexCount();
		PID = (short int )cell.getPID();
		TT.resize(sizeCC);
			
			for(int i=0; i<sizeCC; ++i){
				idV = cell.getVertex(i);
				TT[i] = idV;
				
				if(!mapV.exists(idV))	subG->addVertex(tri->getVertexCoords(idV),idV);
			}
		
		subG->addConnectedCell(TT,eltype,PID,idCell);
		TT.clear();
	}
};








