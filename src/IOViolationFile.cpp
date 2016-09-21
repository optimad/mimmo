/*---------------------------------------------------------------------------*\
 * 
 *  CAMILO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License Commercial (//TODO Temporary header of license)
 *  This file is part of CAMILO.
 *
 *  CAMILO is a commercial software: you do not own rights to redistribute it 
 * 	and/or modify it both in source or pre-build formats
 *  Please contact Optimad offices for any further informations				
 *
 *  You should have received a copy of the Camilo Commercial License
 *  along with CAMILO, as well as the key to unlock the software.
 *
 \ *---------------------------------------------------------------------------*/

#include "IOViolationFile.hpp"
#include <iostream>
#include <fstream>

using namespace mimmo;

//PUBLIC METHODS
/*!
 * Constructor. All default parameters deactivate refinement and coarsening,
 * and the option to write DoF according to surface normals
 * 
 * \param[in] id set id counter of the class. See BaseManipulation::set/getClassCounter method.
 */
IOViolationFile::IOViolationFile(){
	m_name = "MiMMO.IOViolationFile";
	m_dir = ".";
	m_namefile = "Deformation_ViolationReport";
};

/*!
 * Destructor;
 */
IOViolationFile::~IOViolationFile(){
	clear();
};

/*!
 * Copy constructor
 */
IOViolationFile::IOViolationFile(const IOViolationFile & other){
	*this = other;
};

/*!
 * Copy operator of the class
 */
IOViolationFile & IOViolationFile::operator=(const IOViolationFile & other){
	
	*(static_cast<BaseManipulation *> (this)) = *(static_cast<const BaseManipulation *> (&other));
	m_dir = other.m_dir;
	m_namefile = other.m_namefile;;
	return(*this);
};


/*!
 * It builds the input/output ports of the object
 */
void IOViolationFile::buildPorts(){

	bool built = true;

	//input
	built = (built && createPortIn<std::pair<BaseManipulation*, double>, IOViolationFile>(this, &mimmo::IOViolationFile::addViolationData, PortType::M_VIOLATION, mimmo::pin::containerTAG::PAIR, mimmo::pin::dataTAG::PAIRMIMMO_OBJFLOAT_));
	m_arePortsBuilt = built;
};

/*!
 * Get the maximum violation value, amongst those linked to the class
 */
double IOViolationFile::getMaximumViolation(){
	
	double result = -1.E+18;
	
	for(auto &p : m_list){
		result = std::fmax(result, p.second);
	}
	
	return result;
}

/*!
 * Return current stored list of violation values
 */
std::vector<std::pair<std::string, double> > IOViolationFile::getViolationMap(){
	return m_list;
};

/*!
 * Return directory for writing output
 */
std::string	IOViolationFile::getDir(){
	return m_dir;
};

/*!
 * Return name of output file. Tag extension is implicitly defined as .dat..
 */
std::string	IOViolationFile::getFileName(){
	return m_namefile;
};


/*!
 * Set Output directory
 * \param[in] dir 
 */
void	IOViolationFile::setDir(std::string dir){
		m_dir = dir;
}

/*!
 * Set Output filename
 * \param[in] name
 */
void	IOViolationFile::setFileName(std::string name){
	m_namefile = name;
}


/*!
 * Set list of violation values.  
 * \param[in] list violation values with their etiquettes
 */
void	IOViolationFile::setViolationMap(std::vector<std::pair<std::string, double> > list){
		m_list = list;
}

/*!
 * Add a violation value to the general list.  
 * \param[in] data violation value with its etiquette
 */
void	IOViolationFile::addViolationDataString(std::pair<std::string, double> data){
	m_list.push_back(data);
}


/*!
 * Add a violation value to the general list, but using a general BaseManipulation
 * sender object as key. Value will be stored in list with name of BaseManipulation 
 * object as etiquette.  
 * \param[in] data violation value with BaseManipulation object pointer as key 
 */
void	IOViolationFile::addViolationData(std::pair<BaseManipulation *, double> data){
	
	std::pair<std::string, double> result;
	
	std::string tempName ="NONE";
	if(data.first != NULL)	tempName = data.first->getName()+std::to_string(data.first->getClassCounter());
	result = std::make_pair(tempName, data.second);
	m_list.push_back(result);
}



/*!
 * Clear all list of violation values.  
 * \param[in] list violation values with their etiquettes
 */
void	IOViolationFile::clearViolationMap(){
	m_list.clear();
}

/*!
 * Clear any contents of the class and restore defaults
 * 
 */
void 	IOViolationFile::clear(){
	clearViolationMap();
	m_dir = ".";
	m_namefile = "Deformation_ViolationReport";
	
}

/*!
 * Execution command of the class. 
 */
void 	IOViolationFile::execute(){
	
	std::string filename = m_dir + "/"+ m_namefile + ".dat";
	std::ofstream out(filename.c_str());
	
	//check stream if its good
	if(out.is_open()){
	
		out<<""<<std::endl;	
		out<<"$MAXIMUM_VIOLATION"<<'\t'<<std::scientific<<getMaximumViolation()<<std::endl;	
		out<<""<<std::endl;	
		out<<""<<std::endl;	
		
		for(auto &val : m_list){
			out<<"$SLOT"<<'\t'<< val.first<<'\t'<<'\t'<<std::scientific<<val.second<<std::endl;
		}	
	}else{
		std::cout<<"Error of "<<m_name<<" : cannot open "<<filename<< " requested. Exiting... "<<std::endl;
		exit(1);
	}
	
	out.close();

	return;
}

/*!
 * Plot class infos to a XML bitpit::Config::section. The parameters that can be flushed are
 * 
 *
 *  ViolationData are meant to be passed through ports
 * 
 * \param[in] slotXML 	bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void IOViolationFile::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	slotXML.set("ClassName", m_name);
	slotXML.set("ClassID", std::to_string(getClassCounter()));
	
	return;
};

