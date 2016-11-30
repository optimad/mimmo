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
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 * \param[in] formattype set tag format type of your geometries for I/O. See FileType enum
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(int topo, int formattype){
	initializeClass(topo, formattype);
}

/*!Custom constructor of MultipleMimmoGeometries.
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 * \param[in] formattype set tag format type of your geometries for I/O. See FileType enum
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
 * Get how many sub-parts are available for writing, once division map is set. 
 * Return each short id identifying the sub-parts in an unordered_set container
 * Meant only for WRITING mode of the class.
 */
std::unordered_set<short> 
MultipleMimmoGeometries::getHowManySubDivisions(){
	std::unordered_set<short> subtypes;
	for(auto & element: m_mapMGeo){
		subtypes.insert(element.second);
	}
	
	return subtypes;
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
	
	if(getGeometry() == external || external.isEmpty()) return;
	if(external.getType() != m_topo)	return;	
	
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
	
	//check format;
	if(!checkCoherentFormat(getGeometry())){
		std::cout<<"Currently linked geometry is not compatible with file format set for writing"<<std::endl;
		return false;
	}
	if((int)m_winfo.size() < 1){
		std::cout<<"Not filenames found for writing"<<std::endl;
		return false;
	}
	
	
	int nFiles = m_winfo.size();
	std::unordered_set<short>  subIds = getHowManySubDivisions();
	std::unordered_set<short>::iterator itS = subIds.begin();
	int nWriters = subIds.size();
	nWriters = std::min(nWriters, nFiles);

	int counter = 0;
	while(counter< nWriters && itS != subIds.end()){
		
		livector1D cells = cellExtractor(*itS);
		if(!cells.empty()){
			
			std::unique_ptr<MimmoObject> geo(new MimmoObject(m_topo));
			std::unique_ptr<MimmoGeometry> writer(new MimmoGeometry());
			fillSubStructure(cells, geo.get());
			
			writer->setWrite(true);
			writer->setWriteDir(m_winfo[counter].fdir);
			writer->setWriteFilename(m_winfo[counter].fname);
			writer->setWriteFileType(m_winfo[counter].ftype);
			writer->setCodex(m_codex);
			writer->setFormatNAS(m_wformat);
			writer->setGeometry(geo.get());
			
			writer->exec();
			counter++;		
		}
		itS++;
	}
	
};

/*!It reads the mesh geometry from an input file and reverse it in the internal 
 * MimmoObject container. If an external container is linked skip reading and doing nothing.
 * \return False if file doesn't exists or not found geometry container address.
 */
bool
MultipleMimmoGeometries::read(){
	if(!m_isInternal) return false;
	
	setGeometry();
	std::unordered_set<short> & pidMotherList = getGeometry()->getPIDTypeList();
	std::vector<std::unique_ptr<MimmoObject> > listRObjects(m_rinfo.size());

	//absorbing Geometries.
	int counter = 0;
	for(auto & data: m_rinfo){
		if(checkReadingFiles(data)){
			
			std::unique_ptr<MimmoGeometry> geo(new MimmoGeometry());
			std::unique_ptr<MimmoObject> subData(new MimmoObject());	
			geo->setRead(true);
			geo->setReadDir(data.fdir);
			geo->setReadFilename(data.fname);
			geo->setReadFileType(data.ftype);
			
			geo->exec();
			
			subData->setHARDCopy(geo->getGeometry());
			listRObjects[counter] = std::move(subData);
			counter++;
		}
	}
	
	listRObjects.resize(counter);
	
	//put all read geometries in the internal unique MimmoObject of the class.
	// eventual PIDding conflicts are solved.
	
	//find cell and vertex sizes for the unique objects
	long sizeCells=0;
	long sizeVerts=0;
	for(auto & obj: listRObjects){
		sizeVerts += obj->getNVerts();
		if(m_topo !=3)	sizeCells += obj->getNCells();
	}

	//reserving space for internal patch structures
	m_intgeo->getPatch()->reserveVertices(sizeVerts);
	m_intgeo->getPatch()->reserveCells(sizeCells);
	
	//let's roll
	
	long offV = 0;
	long offC = 0;
	long offI = 0;
	
	long id = 0;
	darray3E tempV;
	livector1D tempC;
	short tempPID;
	bitpit::ElementInfo::Type type;
	
	for(auto & obj: listRObjects){
		
		//renumbering vertices and cells of the local patch
		obj->getPatch()->consecutiveRenumber(offV, offC, offI);;
		
		//get a map for eventual renumbering of PIDs
		std::unordered_map<short,short> renumPIDfromlocal;
		
		for(auto & ee: obj->getPIDTypeList()){
			short checker = 0;
			if(pidMotherList.count(*ee) > 0){
				while(pidMotherList.count(checker)> 0 && checker < 32000){checker++;}
				renumPIDfromlocal[*ee] = checker;
			}else{
				renumPIDfromlocal[*ee] = *ee;
			}
		}
		
		//copy vertices
		for(auto & vv : obj->getVertices()){
			id = vv.getId();
			tempV = vv.getCoords();
			m_intgeo->addVertex(temp,id);
		}
		
		//copy connected cells
		for(auto & cc : obj->getCells()){
			id = cc.getId();
			tempC = obj->getCellConnectivity(id);
			type = cc.getType();
			tempPID = cc.getPID();
			
			m_intgeo->addConnectedCell(tempC,type, renumPIDfromlocal[PID], id);
		}
		
		//update renumbering thresholds
		offV += obj->getPatch()->getVertexCount();
		offC += obj->getPatch()->getCellCount();
	}
	
	m_intgeo->cleanGeometry();
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

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class reads the following parameters:
 * 
 * 1) ReadFlag - activate reading mode boolean
 * 2) Read Info data - reading files data
 * 3) WriteFlag - activate writing mode boolean
 * 4) Write Info data - writing files data * 
 * 5) Codex - boolean to write ascii/binary
 * 6) BvTree - evaluate bvTree true/false
 * 7) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void MultipleMimmoGeometries::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	std::string input; 
	std::vector<FileInfoData> temp;
	int counter;
	
	//checking topology
	if(slotXML.hasOption("Topology")){
		std::string input = slotXML.get("Topology");
		input = bitpit::utils::trim(input);
		int temp = -1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		if(m_topo != temp)	return;
	}	

	//checking format type
	if(slotXML.hasOption("Format")){
		std::string input = slotXML.get("Format");
		input = bitpit::utils::trim(input);
		int tempformat = -1;
		if(!input.empty()){
			for(auto c: FileType::_values()){
				if(input == c._to_string()){
					tempformat = FileType::_to_integer(c);
				}
			}
		}
		if(m_tagtype != temp)	return;
	}	
	
	if(slotXML.hasOption("ReadFlag")){
		input = slotXML.get("ReadFlag");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setRead(value);
	}; 
	
	if(slotXML.hasSection("ReadFiles")){
		bitpit::Config::Section &reading = slotXML.getSection("ReadFiles");
		int nfile = reading.getSectionCount();
		temp.resize(nfile);
		
		std::string strdum, root="File";
		
		counter=0;
		
		for(int i=0; i<nfile; ++i){
			strdum = root+std::to_string(i);
			if(reading.hasSection(strdum)){
				bitpit::Config::Section &file = reading.getSection(strdum);
			
				temp[counter].ftype =m_tagtype;

				input = file.get("Dir");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "./";
				temp[counter].fdir = input;

				input = file.get("Name");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "MultipleMimmoGeometries";
				temp[counter].fname = input;
				counter++;
			}
		}
		temp.resize(counter);
		setReadListAbsolutePathFiles();
	}

	
	if(slotXML.hasOption("WriteFlag")){
		input = slotXML.get("WriteFlag");
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss >> value;
		}
		setWrite(value);
	}; 
	
	if(slotXML.hasSection("WriteFiles")){
		bitpit::Config::Section & writing = slotXML.getSection("WriteFiles");
		int nfile = reading.getSectionCount();
		temp.resize(nfile);
		
		std::string strdum, root="File";
		
		counter = 0;
		
		for(int i=0; i<nfile; ++i){
			strdum = root+std::to_string(i);
			
			if(writing.hasSection(strdum)){
				bitpit::Config::Section & file = writing.getSection(strdum);
			
				temp[counter].ftype =m_tagtype;
			
				input = file.get("Dir");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "./";
				temp[counter].fdir = input;
	   
				input = file.get("Name");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "MultipleMimmoGeometries";
				temp[counter].fname = input;
				counter++;
			}
		}
		temp.resize(counter);
		setWriteListAbsolutePathFiles();
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
 * 2) Read Info data - reading files data
 * 3) WriteFlag - activate writing mode boolean
 * 4) Write Info data - writing files data * 
 * 5) Codex - boolean to write ascii/binary
 * 6) BvTree - evaluate bvTree true/false
 * 7) KdTree - evaluate kdTree ture/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot 
 */
void MultipleMimmoGeometries::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	slotXML.set("Topology", m_topo);
	
	std::string typetag = (FileType::_from_integral(m_tagtype))._to_string();
	slotXML.set("Format", typetag);
	
	
	std::string output;
	
	output = std::to_string(m_read);
	slotXML.set("ReadFlag", output);

	if (!m_rinfo.empty()){
		
		bitpit::Config::Section & local = slotXML.addSection("ReadFiles");
		int size = m_rinfo.size();
		std::string strdum, root="File";
		
		for(int i=0; i<size; ++i){
			strdum = root+std::to_string(i);
			bitpit::Config::Section & file = slotXML.addSection(strdum);
			file.set("Dir", m_rinfo.fdir);
			file.set("Name", m_rinfo.fname);
		}
	}
	
	output = std::to_string(m_write);
	slotXML.set("WriteFlag", output);
	
	if (!m_winfo.empty()){
		
		bitpit::Config::Section & local = slotXML.addSection("WriteFiles");
		int size = m_winfo.size();
		std::string strdum, root="File";
		
		for(int i=0; i<size; ++i){
			strdum = root+std::to_string(i);
			bitpit::Config::Section & file = slotXML.addSection(strdum);
			file.set("Dir", m_winfo.fdir);
			file.set("Name", m_winfo.fname);
		}
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
 * Class initializer, meant to be used in construction.
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
			if (formattype !=6 || formattype !=7){
				formattype= 7;
			}
			break;
	}
	
	m_tagtype = formattype;
	
	setDefaults();
};

/*!
 * Extract MimmoObject cell list from m_mapMGeo, given the id sub-part identifier.
 * If m_mapGeo retains cell ids which do not exist in MimmoObject, stash them from
 * the extracted list.
 * \param[in] id  identifier of sub-part in m_mapMGeo
 * \result	list of extracted cells
 */
livector1D
MultipleMimmoGeometries::cellExtractor(short id){
	livector1D result;
	flag = true;
	result.resize(m_mapMGeo.size());
	int counter = 0;
	
	auto & check = getGeometry()->getCells();
	
	for(auto & elem : m_mapMGeo){
		if(elem.second==id && check.exists(elem.first)){
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


/*!
 * Check if current linked MimmoObject has a topology/element type coherent with class tag format m_tagtype
 * \param[in]	geo linked MimmoObject geometry
 * \result		boolean true, if requirements are met.
 */
bool 
MultipleMimmoGeometries::checkCoherentFormat(MimmoObject*){
	
	bool check = true;
	
	//point cloud has no restriction, everything can be written in pointcloud file formats OFP and PCVTU
	int   elesize;
	switch(m_tagtype){
		
		case 0:		elesize = 3;
			break;
			
		case 2:		elesize = 4;
			break;
			
		case 3:		elesize = 5;			
			break;
			
		case 4:		elesize = 8;
			break;
		
		case 1:		elesize = 3;
			break;
		
		case 5:		elesize = 3;
			break;
		
		default: 	elesize = -1;
			break;
	}
	
	if(elesize > 0){
		//check if all elements are triangles
		conns = getGeometry()->getConnectivity();
		for(auto &ele : conns){
			check = check && ((int)ele.size() == elesize);
		}	
	}
	
	return check;
}


/*!
 * Check if filenames for READING, stocked in m_rinfo list are available or not on your system.
 * \param[in] filedata 	data of external files to read from
 * \return boolean true/false for valid file check or not.
 */
bool 
MultipleMimmoGeometries::checkReadingFiles(mimmo::FileInfoData & filedata){
	
	switch(FileType::_from_integral(filedata.ftype)){
		case FileType::STL :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".stl");
			bool check = infile.good();
			if (!check){
				infile.open(filedata.fdir+filedata.fname+".STL");
				check = infile.good();
				if (!check) return false;
			}
		}
		break;
		
		case FileType::STVTU :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".vtu");
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::SQVTU :
		{
			
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".vtu");
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::VTVTU :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".vtu");
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::VHVTU :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".vtu");
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::NAS :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".nas");
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::OFP :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname);
			bool check = infile.good();
			if (!check) return false;
		}
		break;
		
		case FileType::PCVTU :
		{
			std::ifstream infile(filedata.fdir+"/"+filedata.fname+".vtu");
			bool check = infile.good();
			if (!check) return false;
		}
		
		default: //never been reached
			break;
			
	}
	
	return true;
};


