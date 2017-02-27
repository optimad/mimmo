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

REGISTER_MANIPULATOR("MiMMO.MultipleGeometries", "multiplegeometries");

/*!Default constructor of MultipleMimmoGeometries.
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(int topo){
	initializeClass(topo, false);
}

/*!Custom constructor of MultipleMimmoGeometries.
 * Format admissible are linked to your choice of topology. See FileType enum
 * \param[in] topo	set topology of your geometries. 1-surface, 2-volume, 3-pointcloud
 * \param[in] IOMode set boolean to activate reading(false) and writing(true) mode of the class
 */
MultipleMimmoGeometries::MultipleMimmoGeometries(int topo, bool IOMode){
	initializeClass(topo, IOMode);
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
	m_buildBvTree = other.m_buildBvTree;
	m_buildKdTree = other.m_buildKdTree;
	m_topo = other.m_topo;
	m_ftype_allow = other.m_ftype_allow;
	
	m_isInternal = other.m_isInternal;
	m_extgeo = other.m_extgeo;
	
	//warning the internal data structure is not copied. Relaunch the execution eventually to fill it.
	return *this;
};

void
MultipleMimmoGeometries::buildPorts(){
	bool built = true;
	built = (built && createPortIn<std::vector<MimmoObject*>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setGeometry, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortIn<std::unordered_map<std::string,std::pair<int, MimmoObject*> >, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setObjMAP, PortType::M_MAPGEOM, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::STRINGPAIRINTMIMMO_));
	built = (built && createPortIn<std::vector<FileDataInfo>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setReadListFDI, PortType::M_FINFO, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	built = (built && createPortIn<std::vector<FileDataInfo>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::setWriteListFDI, PortType::M_FINFO2, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	
	built = (built && createPortOut<std::vector<MimmoObject*>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getGeometry, PortType::M_VECGEOM, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::MIMMO_));
	built = (built && createPortOut<std::unordered_map<std::string,std::pair<int, MimmoObject*> >, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getObjMAP, PortType::M_MAPGEOM, mimmo::pin::containerTAG::UN_MAP, mimmo::pin::dataTAG::STRINGPAIRINTMIMMO_));
	built = (built && createPortOut<std::vector<FileDataInfo>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getReadListFDI, PortType::M_FINFO, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	built = (built && createPortOut<std::vector<FileDataInfo>, MultipleMimmoGeometries>(this, &mimmo::MultipleMimmoGeometries::getWriteListFDI, PortType::M_FINFO2, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FILEINFODATA));
	
	m_arePortsBuilt = built;
}


/*!
 * Return current format files allowed in your class. Refers to enum FileType
 */
std::vector<int>
MultipleMimmoGeometries::getFileTypeAllowed(){
	ivector1D results(m_ftype_allow.size());
	
	int counter=0;
	for( auto &ee : m_ftype_allow){
		results[counter] = ee;
		++counter;
	}
	
	return results;
}

/*!
 * Return a const pointer to your current class.
 */
const MultipleMimmoGeometries * MultipleMimmoGeometries::getCopy(){
	return this;
}

/*!
 * Get current list of geometry pointers. Reimplementation of BaseManipulation::getGeometry
 */
std::vector<MimmoObject *> 
MultipleMimmoGeometries::getGeometry(){
	
	if(m_isInternal)	{
		std::vector<MimmoObject*> results(m_intgeo.size());
		int counter = 0;
		for(auto & val : m_intgeo){
			results[counter] = val.get();
			++counter;
		}
		return results;
	}
	
	else				return m_extgeo;
}

/*!
 * Get current geometry pointer. Reimplementation of BaseManipulation::getGeometry,
 * const overloading
 */
const std::vector<MimmoObject *>
MultipleMimmoGeometries::getGeometry() const{
	if(m_isInternal)	{
		std::vector<MimmoObject*> results(m_intgeo.size());
		int counter = 0;
		for(auto & val : m_intgeo){
			results[counter] = val.get();
			++counter;
		}
		const std::vector<MimmoObject *> temp(results);
		return temp;
	}
	else return m_extgeo;
}
/*!
 * Get list of external filenames the class must READ from to get the final MimmoObject.
 */
std::vector<FileDataInfo>	
MultipleMimmoGeometries::getReadListFDI(){
	return m_rinfo;
}

/*!
 * Get list of external filenames the class must WRITE to the final MimmoObject.
 */
std::vector<FileDataInfo>	
MultipleMimmoGeometries::getWriteListFDI(){
	return m_winfo;
}


/*!
 * Get the complete map of the name of external file vs relative geometry object. 
 * If the class is in reading mode, return reading pathfiles and geometries
 * (these are usually available after the class execution). In writing mode, return writing pathfiles
 * and geometry pointers, that usually are linked by the User as inputs.
 */
std::unordered_map<std::string, std::pair<int,MimmoObject*> >
MultipleMimmoGeometries::getObjMAP(){
	
	std::unordered_map<std::string, std::pair<int, MimmoObject*> > objMap;
	
	std::string fullpath;
	std::vector<FileDataInfo> info;
	if(m_read){
		info = m_rinfo;
	}else{
		info = m_winfo;
	}	
	
	int totVal = info.size();
	int totGeo = m_extgeo.size();
	if(m_isInternal)	totGeo = m_intgeo.size();

	totVal = std::min(totVal,totGeo);
	
	for(int i=0; i<totVal; ++i){
		
		std:string tag = "vtu";
		switch(info[i].ftype){
			case 0: tag="stl"; break;
			case 5: tag="nas"; break;
			case 6: tag="";	   break;
			default:	break;
		}
		
		fullpath = info[i].fdir+"/"+info[i].fname+"."+tag;
		
		if(m_isInternal)	objMap[fullpath] = std::make_pair(info[i].ftype, m_intgeo[i].get());
		else				objMap[fullpath] = std::make_pair(info[i].ftype, m_extgeo[i]);
	}
	return objMap;
}


/*!
 * Add an external file path to the list of geometries to read from.
 * Available format type are related to the type of geometry topology 
 * set for your class and the format fixed for your class also. See FileType enum.
 * \param[in] dir  string, absolute path directory
 * \param[in] name string, name of your file, without tag extension
 * \param[in] ftype enum FileType, extension tag of your geometry file
 */
void
MultipleMimmoGeometries::setAddReadFile(std::string dir, std::string name, FileType ftype){
	
	int fty = ftype._to_integral();
	if(m_ftype_allow.count(fty)== 0) return;
	
	FileDataInfo temp;
	temp.ftype	= fty;
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
 * \param[in] ftype enum FileType, extension tag of your geometry file
 */
void
MultipleMimmoGeometries::setAddWriteFile(std::string dir, std::string name, FileType ftype){

	int fty = ftype._to_integral();
	if(m_ftype_allow.count(fty)== 0) return;
	
	FileDataInfo temp;

	temp.ftype	= fty;
	temp.fdir	= dir;
	temp.fname	= name;
	
	m_winfo.push_back(temp);
};

/*!
 * Set the whole list of external files data which the geometry is read from.
 */
void
MultipleMimmoGeometries::setReadListFDI(std::vector<FileDataInfo> data){
	m_rinfo.clear();
	m_rinfo.resize(data.size());

	int counter = 0;
	for(auto & ele : data){
		
		if(m_ftype_allow.count(ele.ftype) > 0){
			m_rinfo[counter] = ele;
			counter++;
		}	
	}
};

/*!
 * Set the whole list of external files data which the geometry is written to.
 */
void
MultipleMimmoGeometries::setWriteListFDI(std::vector<FileDataInfo> data){
	m_winfo.clear();
	m_winfo.resize(data.size());
	
	int counter = 0;
	for(auto & ele : data){
		if(m_ftype_allow.count(ele.ftype) > 0){
			m_winfo[counter] = ele;
			counter++;
		}	
	}
	
};

/*!
 * Set directly the object map of path to external files and pointers of geometry.
 * Whenever you have ready a list of external files and a list of pointers to MimmoObject 
 * which contains geometrical data of these files, you can use uniquely the current method instead of 
 * setWriteListFDI() and setGeometry(). This method is meant for writing mode only. In reading mode
 * this method do nothing. 
 */
void
MultipleMimmoGeometries::setObjMAP(std::unordered_map<std::string, std::pair<int, MimmoObject*> > map){
	
	if(m_read)	return;
	
	m_winfo.resize(map.size());
	m_extgeo.resize(map.size());
	
	int counter = 0;
	std::string key1="/\\", key2=".";
	std::string name, tag;
	std::size_t found;
	
	for(auto & val : map){
	
		if(m_ftype_allow.count(val.second.first) > 0){
		
			m_winfo[counter].ftype = val.second.first;
			
			found = val.first.find_last_of(key1);
			m_winfo[counter].fdir = val.first.substr(0,found);

			name = val.first.substr(found+1);
			found = name.find_last_of(key2);
			m_winfo[counter].fname = name.substr(0, found);
			
			m_extgeo[counter] = val.second.second;
			++counter;
		}
	}
	
	m_winfo.resize(counter);
	m_extgeo.resize(counter);
};


/*!It sets the condition to read the geometry on file during the execution.
 * \param[in] read Does it read the geometry in execution?
 */
void
MultipleMimmoGeometries::setRead(bool read){
	m_read = read;
	m_write = !read;
}

/*!It sets the condition to write the geometry on file during the execution.
 * \param[in] write Does it write the geometry in execution?
 */
void
MultipleMimmoGeometries::setWrite(bool write){
	m_write = write;
	m_read = !write;
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
 * Soft copy means that the internal member m_intgeo is not copied, but to get it
 * you need to relaunch the class execution at a certain point.
 * On the contrary, other members are exactly copied
 * \param[in] other pointer to MultipleMimmoGeometries class.
 */
void
MultipleMimmoGeometries::setSOFTCopy(const MultipleMimmoGeometries * other){
	clear();	
	*this = *other;
}

/*!Sets your current class as a "HARD" copy of the argument.
 * Hard copy means that the current geometric object MimmoObject is
 * copied as a stand alonelist of  internal members and stored in the unique pointer list m_intgeo. 
 * Other members are exactly copied
 * \param[in] other pointer to MultipleMimmoGeometries class.
 */
void
MultipleMimmoGeometries::setHARDCopy(const MultipleMimmoGeometries * other){

	clear();
	*(static_cast<BaseManipulation * >(this)) = *(static_cast<const BaseManipulation * >(other));

	m_rinfo = other->m_rinfo;
	m_winfo = other->m_winfo;
	m_read = other->m_read;
	m_write = other->m_write;
	m_wformat = other->m_wformat;
	m_codex = other->m_codex;
	
	m_topo = other->m_topo;
	m_ftype_allow = other->m_ftype_allow;
	
	m_buildBvTree = other->m_buildBvTree;
	m_buildKdTree = other->m_buildKdTree;
	m_extgeo = other->m_extgeo;

	m_isInternal = other->m_isInternal;
	
	if(m_isInternal){
		int size = other->m_intgeo.size();
		m_intgeo.resize(size);
		
		for(int i=0; i<size; ++i){
			std::unique_ptr<MimmoObject> dum (new MimmoObject());
			dum->setHARDCopy(other->m_intgeo[i].get());
			m_intgeo[i] = std::move(dum);
		}
	}
}

/*!
 * Set geometry from an external MimmoObject source, softly linked. 
 * Reimplementation of BaseManipulation::setGeometry
 */
void
MultipleMimmoGeometries::setGeometry(std::vector<MimmoObject *> external){
	
	if(external.empty()) return;
	
	m_extgeo.resize(external.size());
	int counter = 0;
	for(auto & obj: external){
		if(!obj->isEmpty()){
			m_extgeo[counter]= obj;
			++counter;
		}
	}
	m_extgeo.resize(counter);	
	m_intgeo.clear();
	m_isInternal = false;
};

/*!
 * Force your class to allocate a list of internal MimmoObject of m_topo type. 
 * Other internal object allocated or externally linked geometries
 * will be destroyed/unlinked. This option is meant for reading mode of the class only
 * 
 */
void
MultipleMimmoGeometries::setGeometry(){
	if(m_write)	return;
	
	m_extgeo.clear();
	m_intgeo.clear();
	int size = m_rinfo.size();
	m_intgeo.resize(size);
	
	for(int i=0; i<size; ++i){
		std::unique_ptr<MimmoObject> dum(new MimmoObject(m_topo));
		m_intgeo[i] = std::move(dum);
	}	
	m_isInternal = true;
};

/*!It sets the format to export .nas files if any in WRITING mode.
 * \param[in] wform Format of .nas file (Short/Long).
 */
void
MultipleMimmoGeometries::setFormatNAS(WFORMAT wform){
	m_wformat = wform;
}

/*!It sets if the BvTree of all the patch geometries has to be built during execution.
 * \param[in] build If true the BvTree is built in execution and stored in
 * the related MimmoObject member.
 */
void
MultipleMimmoGeometries::setBuildBvTree(bool build){
	m_buildBvTree = build;
}

/*!It sets if the KdTree of all the patch geometries has to be built during execution.
* \param[in] build If true the KdTree is built in execution and stored in
* the related MimmoObject member.
*/
void
MultipleMimmoGeometries::setBuildKdTree(bool build){
	m_buildKdTree = build;
}

/*!
 * Check if geometries are not linked or not locally instantiated in your class.
 * True - no geometries present, False otherwise.
 */
bool MultipleMimmoGeometries::isEmpty(){
	return (m_extgeo.empty() && m_intgeo.empty());
}

/*!
 * Check if geometries are internally instantiated (true) or externally linked(false).
 * Return false if no geometry is checked. Please verify it with isEmpty method first
 */
bool MultipleMimmoGeometries::isInternal(){
	return (m_isInternal);
}

/*!
 * Check which current mode of the class is active, Reading-False, Writing-True 
 */
bool
MultipleMimmoGeometries::whichMode(){
	return m_write;
}

/*!
 * Clear all stuffs in your class
 */
void
MultipleMimmoGeometries::clear(){
	clearReadPaths();
	clearWritePaths();
	setDefaults();
	m_intgeo.clear();
	m_extgeo.clear();
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

/*!It writes the listed geometries on multiple output files.
 *\return False if no geometries or relative output info are set.
 */
bool
MultipleMimmoGeometries::write(){
	if (m_extgeo.empty() || m_winfo.empty()) return false;
	
	int totGeo = m_extgeo.size();
	int totInfo = m_winfo.size();
	
	if(totGeo != totInfo)	{
		std::cout<<"WARNING:Not enough output info or some invalid geoemetries found writing in class "<<m_name<<std::endl;
	}
	
	totGeo = std::min(totGeo, totInfo);
	
	for(int k=0; k<totGeo; ++k){
		
		std::unique_ptr<MimmoGeometry> writer(new MimmoGeometry());
		writer->setWrite(true);
		writer->setWriteDir(m_winfo[k].fdir);
		writer->setWriteFilename(m_winfo[k].fname);
		writer->setWriteFileType(m_winfo[k].ftype);
		writer->setCodex(m_codex);
		writer->setFormatNAS(m_wformat);
		writer->setGeometry(m_extgeo[k]);
		writer->execute();
	}
	return true;
};

/*!It reads the mesh geometries from a list of input files and put them in the internal 
 * MimmoObject list container. If an external container is linked, skip reading and do nothing.
 * \return False if files do not exist or not found geometry container to address to.
 */
bool
MultipleMimmoGeometries::read(){
	if(!m_isInternal || m_rinfo.empty()) return false;
	
	setGeometry();

	//absorbing Geometries.
	int counter = 0;
	m_intgeo.resize(m_rinfo.size());
	for(auto & data: m_rinfo){
		std::unique_ptr<MimmoGeometry> geo(new MimmoGeometry());
		std::unique_ptr<MimmoObject> subData(new MimmoObject());	
		geo->setRead(true);
		geo->setReadDir(data.fdir);
		geo->setReadFilename(data.fname);
		geo->setReadFileType(data.ftype);
		geo->execute();
			
		subData->setHARDCopy(geo->getGeometry());
		m_intgeo[counter] = std::move(subData);
		counter++;
	}
	if(counter == 0) return false;
	m_intgeo.resize(counter);
	
	for(auto & val: m_intgeo){
		val->cleanGeometry();
		if(m_buildBvTree)	val->buildBvTree();
		if(m_buildKdTree)	val->buildKdTree();
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
		std::cout << "MiMMO : ERROR : no files to read found : "<< std::endl;
		std::cout << " " << std::endl;
		exit(10);
	}
	check = true;
	if (m_write) check = write();
	if (!check){
		std::cout << "MiMMO : ERROR : write not done : geometry not linked/ write-paths not specified/or incompatible formats " << std::endl;
		std::cout << " " << std::endl;
		exit(11);
	}
	return;
}

/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML. Except of geometry parameter (which is instantiated internally
 * or passed by port linking), the class reads the following parameters:
 * 
 * --> Absorbing data:
 * ReadFlag : activate reading mode boolean
 * Read Info data : reading files data
 * WriteFlag : activate writing mode boolean
 * Write Info data : writing files data  
 * Codex : boolean to write ascii/binary
 * BvTree : evaluate bvTree true/false
 * KdTree : evaluate kdTree true/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void MultipleMimmoGeometries::absorbSectionXML(bitpit::Config::Section & slotXML, std::string name){

	std::string input; 
	std::vector<FileDataInfo> temp;
	int counter;
	
	//checking topology
	if(slotXML.hasOption("Topology")){
		std::string input = slotXML.get("Topology");
		input = bitpit::utils::trim(input);
		int temptop = -1;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temptop;
		}
		if(m_topo != temptop)	return;
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
	
	if(slotXML.hasSection("ReadInfoData")){
		bitpit::Config::Section &reading = slotXML.getSection("ReadInfoData");
		int nfile = reading.getSectionCount();
		temp.resize(nfile);
		counter=0;
		
		for(auto & file : reading.getSections()){
			
			if(!file.second->hasOption("tag") || !file.second->hasOption("dir") || !file.second->hasOption("name")) continue;
			input = file.second->get("tag");
			input = bitpit::utils::trim(input);
			auto maybe_tag = FileType::_from_string_nothrow(input.c_str());
				
			if(!maybe_tag)	continue;
			if(m_ftype_allow.count(maybe_tag->_to_integral()) > 0){
					
				temp[counter].ftype = maybe_tag->_to_integral();
				
				input = file.second->get("dir");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "./";
				temp[counter].fdir = input;

				input = file.second->get("name");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "MultipleMimmoGeometries";
				temp[counter].fname = input;
				
				counter++;
			}	
		}
		temp.resize(counter);
		setReadListFDI(temp);
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
	
	if(slotXML.hasSection("WriteInfoData")){
		bitpit::Config::Section & writing = slotXML.getSection("WriteInfoData");
		int nfile = writing.getSectionCount();
		temp.resize(nfile);
		counter = 0;
		
		for(auto & file : writing.getSections()){
			
			if(!file.second->hasOption("tag") || !file.second->hasOption("dir") || !file.second->hasOption("name")) continue;
				
			input = file.second->get("tag");
			input = bitpit::utils::trim(input);
			auto maybe_tag = FileType::_from_string_nothrow(input.c_str());
				
			if(!maybe_tag)	continue;
			if(m_ftype_allow.count(maybe_tag->_to_integral()) > 0){
					
				temp[counter].ftype = maybe_tag->_to_integral();
					
				input = file.second->get("dir");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "./";
				temp[counter].fdir = input;
	   
				input = file.second->get("name");
				input = bitpit::utils::trim(input);
				if(input.empty())	input = "MultipleMimmoGeometries";
				temp[counter].fname = input;
				counter++;
			}
		}
		temp.resize(counter);
		setWriteListFDI(temp);
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
 * --> Flushing data// how to write it on XML:
 * ClassName : name of the class as "MiMMO.MultipleGeometries"
 * ClassID	  : integer identifier of the class	
 * ReadFlag : activate reading mode boolean
 * ReadInfodata : reading files data
 * 		<ReadInfoData>
 * 			<File0>
 * 				<tag> format tag extension </tag>
 * 				<dir> path to work directory </dir>
 * 				<name> name of file to read </name>
 * 			</File0>
 * 			<File1>
 * 				<tag> format tag extension </tag>
 * 				<dir> path to work directory </dir>
 * 				<name> name of file to read </name>
 * 			</File1>
 * 			...
 * 			...
 * 		</ReadInfoData>
 * WriteFlag : activate writing mode boolean
 * WriteInfodata : writing files data  
 * 		<WriteInfoData>
 * 			<File0>
 * 				<tag> format tag extension </tag>
 * 				<dir> path to work directory </dir>
 * 				<name> name of file to write </name>
 * 			</File0>
 * 			<File1>
 * 				<tag> format tag extension </tag>
 * 				<dir> path to work directory </dir>
 * 				<name> name of file to write </name>
 * 			</File1>
 * 			...
 * 			...
 * 		</WriteInfoData>
 * 
 * Codex : boolean to write ascii/binary
 * BvTree : evaluate bvTree true/false
 * KdTree : evaluate kdTree true/false
 * 
 * \param[in]	slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot 
 */
void MultipleMimmoGeometries::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	slotXML.set("Topology", m_topo);
	
	std::string output;
	
	output = std::to_string(m_read);
	slotXML.set("ReadFlag", output);

	if (!m_rinfo.empty()){
		
		bitpit::Config::Section & local = slotXML.addSection("ReadInfoData");
		int size = m_rinfo.size();
		std::string strdum, root="File";
		
		for(int i=0; i<size; ++i){
			strdum = root+std::to_string(i);
			bitpit::Config::Section & file = slotXML.addSection(strdum);
			file.set("dir", m_rinfo[i].fdir);
			file.set("name", m_rinfo[i].fname);
			file.set("tag", (FileType::_from_integral(m_rinfo[i].ftype))._to_string());
		}
	}
	
	output = std::to_string(m_write);
	slotXML.set("WriteFlag", output);
	
	if (!m_winfo.empty()){
		
		bitpit::Config::Section & local = slotXML.addSection("WriteInfoData");
		int size = m_winfo.size();
		std::string strdum, root="File";
		
		for(int i=0; i<size; ++i){
			strdum = root+std::to_string(i);
			bitpit::Config::Section & file = slotXML.addSection(strdum);
			file.set("dir", m_winfo[i].fdir);
			file.set("name", m_winfo[i].fname);
			file.set("tag", (FileType::_from_integral(m_winfo[i].ftype))._to_string());
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
	m_read		= true;
	clearReadPaths();
	m_write		= false;
	clearWritePaths();
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
MultipleMimmoGeometries::initializeClass(int topo, bool IOMode){
	
	m_name 		= "MiMMO.MultipleGeometries";
	m_read = IOMode; m_write = !IOMode;
	
	m_topo     = std::min(1, topo);
	if(m_topo > 3)	m_topo = 1;
	
	//checking admissible format
	m_ftype_allow.clear();
	switch(m_topo){
		case 1:
			m_ftype_allow.insert(0);
			m_ftype_allow.insert(1);
			m_ftype_allow.insert(2);
			m_ftype_allow.insert(5);
			break;
		case 2:
			m_ftype_allow.insert(3);
			m_ftype_allow.insert(4);
			break;
		default:
			m_ftype_allow.insert(6);
			m_ftype_allow.insert(7);
			break;
	}
	
	setDefaults();
};



