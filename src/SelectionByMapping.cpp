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
 \ *---------------------------------------------------------------------------*/

#include "MeshSelection.hpp"
#include "levelSet.hpp"
#include <cstddef>
using namespace mimmo;


//------------------------------------------------------------------------
//SELECTION	BY MAPPING class 	******************************************
//------------------------------------------------------------------------

REGISTER(BaseManipulation, SelectionByMapping, "MiMMO.SelectionByMapping");

/*!
 * Basic Constructor. Need to know kind of topology chosen 1-3D surface, 2-VolumeMesh
 * . Other options are not available, if forced trigger default value of 1.
 */
SelectionByMapping::SelectionByMapping(int topo){
	m_name = "MiMMO.SelectionByMapping";
	m_type = SelectionType::MAPPING;
	m_tolerance = 1.E-08;
	
	topo = std::min(1,topo);
	if(topo > 2)	topo=1;
	m_topo = topo;
	
	
	m_allowedType.resize(3);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	m_allowedType[2].insert(FileType::VTVTU);
	m_allowedType[2].insert(FileType::VHVTU);
	
	buildPorts();
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
SelectionByMapping::SelectionByMapping(const bitpit::Config::Section & rootXML){
	
	m_name = "MiMMO.SelectionByMapping";
	m_type = SelectionType::MAPPING;
	m_tolerance = 1.E-08;

	std::string fallback_name = "ClassNONE";	
	std::string fallback_topo = "-1";	
	std::string input_name = rootXML.get("ClassName", fallback_name);
	input_name = bitpit::utils::trim(input_name);
	
	std::string input_topo = rootXML.get("Topology", fallback_topo);
	input_topo = bitpit::utils::trim(input_topo);
	
	int topo = std::atoi(input_topo.c_str());
	topo = std::min(1,topo);
	if(topo > 2)	topo=1;
	m_topo = topo;
	
	m_allowedType.resize(3);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	m_allowedType[2].insert(FileType::VTVTU);
	m_allowedType[2].insert(FileType::VHVTU);
	
	buildPorts();
	
	if(input_name == "MiMMO.SelectionByMapping"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml MiMMO::SelectionByMapping constructor. No valid xml data found"<<std::endl;
	};
}

/*!
 * Custom Constructor. Non null geometry set also the topology allowed bu the class.
 * \param[in]	geolist	list of external files to read from, for comparison purposes 
 * \param[in]	target	pointer to target geometry
 * \param[in]	tolerance	proximity criterium tolerance 
 */
SelectionByMapping::SelectionByMapping(std::unordered_map<std::string, int> & geolist, MimmoObject * target, double tolerance){
	m_name = "MiMMO.SelectionByMapping";
	m_type = SelectionType::MAPPING;
	m_tolerance = 1.E-08;
	
	m_allowedType.resize(3);
	m_allowedType[1].insert(FileType::STL);
	m_allowedType[1].insert(FileType::STVTU);
	m_allowedType[1].insert(FileType::SQVTU);
	m_allowedType[1].insert(FileType::NAS);
	m_allowedType[2].insert(FileType::VTVTU);
	m_allowedType[2].insert(FileType::VHVTU);
	
	if(target != NULL){
		m_topo = target->getType();
		m_topo = std::min(1, m_topo);
		if(m_topo > 2)	m_topo = 1;
		setGeometry(target);
		setFiles(geolist);
		m_tolerance = tolerance;
	}
	
	buildPorts();
};

/*!
 * Destructor
 */
SelectionByMapping::~SelectionByMapping(){};

/*!
 * copy Constructor
 */
SelectionByMapping::SelectionByMapping(const SelectionByMapping & other){
	*this = other;
};

/*!
 * Copy Operator
 */
SelectionByMapping & SelectionByMapping::operator=(const SelectionByMapping & other){
	*(static_cast<GenericSelection * >(this)) = *(static_cast<const GenericSelection * >(&other));
	m_tolerance = other.m_tolerance;
	m_geolist = other.m_geolist;
	m_allowedType = other.m_allowedType;
	return *this;
};

/*!
 * It builds the input/output ports of the object
 */
void SelectionByMapping::buildPorts(){

	bool built = true;

	//inheritance
	GenericSelection::buildPorts();

	m_arePortsBuilt = built;
};

/*!
 * Return tolerance actually set in your class
 */
double	SelectionByMapping::getTolerance(){
	return	m_tolerance;
};

/*!
 * Set Proximity tolerance. Under this threshold compared geometries are considered coincident
 */
void 	SelectionByMapping::setTolerance(double tol){
	if(tol == 0.0){
		tol = 1.e-8;
	}
	m_tolerance = tol;
};

/*!
 * Set link to target geometry for your selection. Reimplementation of 
 * GenericSelection::setGeometry();
 */
void	SelectionByMapping::setGeometry( MimmoObject * target){
	
	if(target->getType() != m_topo){
		std::cout<<"SelectionMapping Cannot support current geometry. Topology not supported."<<std::endl;
		return;
	}
	m_geometry = target;
	
};

/*!
 * Return the actual list of external geometry files to read from and compare your target geometry for mapping purposes
 */
const std::unordered_map<std::string, int> & SelectionByMapping::getFiles() const{
	return	m_geolist;
}

/*!
 * Set a list of external geometry files to read from and compare your target geometry for mapping purposes
 * \param[in] files list external geometries to be read
 */
void	SelectionByMapping::setFiles(std::unordered_map<std::string, int>  files){
	for(auto && val : files){
		addFile(val);
	}	
};

/*!
 * Add a new file to a list of external geometry files
 *\param[in] file of external geometry to be read
 */
void 	SelectionByMapping::addFile(std::pair<std::string, int> file){
	int type = m_topo;
	if(m_allowedType[type].find(file.second) != m_allowedType[type].end()){									
		m_geolist.insert(file);
	}	
};

/*!
 * Remove an existent file to a list of external geometry files. If not in the list, do nothing 
 * \param[in] file to be removed from the list 
 */
void 	SelectionByMapping::removeFile(std::string file){
	 if(m_geolist.find(file) != m_geolist.end())	m_geolist.erase(file);
};

/*!
 * Empty your list of file for mapping 
 */
void 	SelectionByMapping::removeFiles(){
	m_geolist.clear();
};

/*!
 * Clear your class
 */
void SelectionByMapping::clear(){
	m_subpatch.reset(nullptr);
	removeFiles();
	m_topo = 0;
	BaseManipulation::clear();
};

/*!
 * Extract portion of target geometry which are enough near to the external geometry provided
 * \return ids of cell of target tesselation extracted 
 */
livector1D SelectionByMapping::extractSelection(){
	
	if(!(getGeometry()->isBvTreeBuilt()))	getGeometry()->buildBvTree();
	std::set<long> cellList;
	
	for (auto && file : m_geolist){
		livector1D list = getProximity(file);
		cellList.insert(list.begin(), list.end());
	}
	
	//check if dual selection is triggered
	livector1D result;
	if(m_dual){
		livector1D totID = getGeometry()->getMapCell();
		result.resize(totID.size() - cellList.size());
		if(result.size() == 0) return result;
		
		std::sort(totID.begin(), totID.end());
		int counter = 0;
		auto tot_it  = totID.begin();
		auto cand_it = cellList.begin();
		while(tot_it != totID.end()){
			long val = *tot_it;
			if (cand_it == cellList.end() || val != *cand_it) {
				result[counter] = val;
				++counter;
			} else {
				++cand_it;
			}
			++tot_it;
		}
	}else{
		result.insert(result.end(), cellList.begin(), cellList.end());
	}	
	return	result;
};

/*!
 * Return portion of target geometry near to an external geometry
 * \param[in] file of external geometry to be compared.
 * \param[in] nbins	number of total raw points for level set evaluation
 */
livector1D SelectionByMapping::getProximity(std::pair<std::string, int> val){
	
	svector1D info = extractInfo(val.first);
	
	MimmoGeometry * geo = new MimmoGeometry();
	geo->setRead(true);
	geo->setWrite(false);
	geo->setReadDir(info[0]);
	geo->setReadFilename(info[1]);
	geo->setReadFileType(val.second);
	geo->setBuildBvTree(true);
	geo->execute();

	if(geo->getGeometry()->getNVertex() == 0 || geo->getGeometry()->getNCells() == 0 ){ 
		std::cout<<"failed to read geometry in SelectionByMapping::getProximity"<<std::endl;
		return livector1D();
	}	
	livector1D result = mimmo::bvTreeUtils::selectByPatch(geo->getGeometry()->getBvTree(), getGeometry()->getBvTree(), m_tolerance);
	delete geo; 
	geo=NULL;
	
	return	result;
};

/*!
 * Extract root dir/filename/tag from an absolute file pattern
 */
svector1D SelectionByMapping::extractInfo(std::string file){
	
	std::string root, name, tag,temp;
	std::string key1=".", key2="/\\";
	
	std::size_t found = file.find_last_of(key2); 
	root = file.substr(0, found);
	temp = file.substr(found+1);
	
	found = temp.find_last_of(key1);
	name = temp.substr(0,found);
	tag = temp.substr(found+1);
	
	svector1D result(3);
	result[0] = root;
	result[1] = name;
	result[2] = tag;
	
	return 	result;
}


/*!
 * Get infos from a XML bitpit::Config::section. The parameters available are
 * 
 * --> Absorbing data:
 * Priority  : uint marking priority in multi-chain execution;
 * Dual       : boolean to get straight what given by selection method or its exact dual
 * Topology   : set Topology for your mapping 1- 3D surface, 2- Volume mesh, 0 none;
 * Tolerance  : tolerance for detect proximity volume in which perform mapping 
 * Files	  : list by their filepaths of external geometries to be mapped on target one 
 * PlotInExecution : boolean 0/1 print optional results of the class.
 * OutputPlot : target directory for optional results writing. 
 * 
 * 
 * Geometry is mandatorily passed through ports. 
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByMapping::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

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
	
	if(slotXML.hasOption("Priority")){
		std::string input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
	//start absorbing
	if(slotXML.hasOption("Dual")){
		std::string input = slotXML.get("Dual");
		input = bitpit::utils::trim(input);
		bool value = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss >> value;
		}
		setDual(value);
	}

	if(slotXML.hasOption("Tolerance")){
		std::string input = slotXML.get("Tolerance");
		input = bitpit::utils::trim(input);
		double temp = 1.0E-8;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
			setTolerance(temp);
			
	}
	
	std::unordered_map<std::string, int> mapp;
	if(slotXML.hasSection("Files")){
		
		const bitpit::Config::Section & filesXML = slotXML.getSection("Files");
		
		for(auto & subfile : filesXML.getSections()){
			std::string path;
			std::string tag;
			
			if(subfile.second->hasOption("fullpath"))	{
				path = subfile.second->get("fullpath");
				path = bitpit::utils::trim(path);
			}
		    if(subfile.second->hasOption("tag")){
				tag = subfile.second->get("tag");
				tag = bitpit::utils::trim(tag);
				//check tag;
				auto maybe_tag = FileType::_from_string_nothrow(tag.c_str());
				if(!maybe_tag)	tag.clear();
				else	tag = maybe_tag->_to_string();
			}	
			
			if(!path.empty() && !tag.empty()){
				mapp[path] = (int) FileType::_from_string(tag.c_str());
			}	
		}

		setFiles(mapp);
		
	}	
	
		if(slotXML.hasOption("PlotInExecution")){
			std::string input = slotXML.get("PlotInExecution");
			input = bitpit::utils::trim(input);
			bool value = false;
			if(!input.empty()){
				std::stringstream ss(input);
				ss >> value;
			}
			setPlotInExecution(value);
		}
		
		if(slotXML.hasOption("OutputPlot")){
			std::string input = slotXML.get("OutputPlot");
			input = bitpit::utils::trim(input);
			std::string temp = ".";
			if(!input.empty())	setOutputPlot(input);
			else			  	setOutputPlot(temp);
		}
	
	return;	
};

/*!
 * Plot infos from a XML bitpit::Config::section. The parameters available are
 * 
 * --> Flushing data// how to write it on XML:
 * ClassName : name of the class as "MiMMO.SelectionByMapping"
 * Priority  : uint marking priority in multi-chain execution; 
 * Dual       : boolean to get straight what given by selection method or its exact dual
 * Topology   : set Topology for your mapping 1- 3D surface, 2- Volume mesh, 0 none;
 * Tolerance  : tolerance for detect proximity volume in which perform mapping 
 * Files	  : list by their filepaths of external geometries to be mapped on target one 
 * 				<Files>
 * 					<file0>
 * 						<fullpath> path to your file </fullpath>
 * 						<tag>	extension tag for your file </tag>
 * 					</file0>
 * 					<file1>
 * 						<fullpath> path to your file </fullpath>
 * 						<tag>	extension tag for your file </tag>
 * 					</file1>
 *					...
 *					... 
 * 				</Files>
 * PlotInExecution : boolean 0/1 print optional results of the class.
 * OutputPlot : target directory for optional results writing. 
 *  
 * Geometry is mandatorily passed through ports. 
 *  
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void SelectionByMapping::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));

	slotXML.set("Topology", m_topo);
	
	int value = m_dual;
	slotXML.set("Dual", std::to_string(value));
	
	
	if(m_tolerance != 1.E-08){
		std::stringstream ss;
		ss<<std::scientific<<m_tolerance;
		slotXML.set("Tolerance", ss.str());
	}
	
	bitpit::Config::Section & filesXML = slotXML.addSection("Files");
	
	int counter = 0;
	for(auto & file : m_geolist){
		std::string name = "file"+std::to_string(counter);
		bitpit::Config::Section & local = filesXML.addSection(name);
		local.set("fullpath", file.first);
		std::string typetag = (FileType::_from_integral(file.second))._to_string(); 
		local.set("tag", typetag);
		++counter;
	}
	
		if(isPlotInExecution()){
			slotXML.set("PlotInExecution", std::to_string(1));
		}
		
		if(m_outputPlot != "."){
			slotXML.set("OutputPlot", m_outputPlot);
		}
	
	return;
};


