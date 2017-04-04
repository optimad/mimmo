/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "GenericInput.hpp"
#include "Operators.hpp"
#include <fstream>

using namespace std;
namespace mimmo {

/*!
 * Default constructor of GenericInput.
 * \param[in] readFromFile True if the object reads the values from file (default value false).
 */
GenericInput::GenericInput(bool readFromFile){
	m_readFromFile = readFromFile;
	m_portsType 	= BaseManipulation::ConnectionType::FORWARD;
	m_name 			= "mimmo.GenericInput";
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
GenericInput::GenericInput(const bitpit::Config::Section & rootXML){
	
	m_readFromFile  = false;
	m_portsType 	= BaseManipulation::ConnectionType::FORWARD;
	m_name 			= "mimmo.GenericInput";
	
	std::string fallback_name = "ClassNONE";	
	std::string input = rootXML.get("ClassName", fallback_name);
	input = bitpit::utils::trim(input);
	if(input == "mimmo.GenericInput"){
		absorbSectionXML(rootXML);
	}else{	
		std::cout<<"Warning in custom xml mimmo::GenericInput constructor. No valid xml data found"<<std::endl;
	};
}

/*!Custom constructor of GenericInput.
 * \param[in] filename Name of the input file.
 * The m_readFromFile flag is set to true.
 */
GenericInput::GenericInput(std::string filename){
	m_readFromFile 	= true;
	m_filename 		= filename;
	m_portsType 	= BaseManipulation::ConnectionType::FORWARD;
};

GenericInput::~GenericInput(){};

/*!Copy constructor of GenericInput.
 */
GenericInput::GenericInput(const GenericInput & other):BaseManipulation(other){
	m_readFromFile 	= other.m_readFromFile;
	m_filename 		= other.m_filename;
};

/*!Assignement operator of GenericInput.
 */
GenericInput & GenericInput::operator=(const GenericInput & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_readFromFile 	= other.m_readFromFile;
	m_filename 		= other.m_filename;
	return *this;
};

/*!It sets if the object imports the displacements from an input file.
 * \param[in] readFromFile True if the object reads the values from file.
 */
void
GenericInput::setReadFromFile(bool readFromFile){
	m_readFromFile = readFromFile;
};

/*!It sets the name of the input file.
 * \param[in] filename Name of the input file.
 */
void
GenericInput::setFilename(std::string filename){
	m_filename = filename;
};

/*! It builds the input/output ports of the object
 */
void
GenericInput::buildPorts(){
	bool built = true;
	built = (built && createPortOut<dvecarr3E, GenericInput>(this, &mimmo::GenericInput::getResult<dvecarr3E>, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvecarr3E, GenericInput>(this, &mimmo::GenericInput::getResult<dvecarr3E>, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvector1D, GenericInput>(this, &mimmo::GenericInput::getResult<dvector1D>, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<dvector1D, GenericInput>(this, &mimmo::GenericInput::getResult<dvector1D>, PortType::M_SCALARFIELD, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, GenericInput>(this, &mimmo::GenericInput::getResult<darray3E>, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<darray3E, GenericInput>(this, &mimmo::GenericInput::getResult<darray3E>, PortType::M_SPAN, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<iarray3E, GenericInput>(this, &mimmo::GenericInput::getResult<iarray3E>, PortType::M_DIMENSION, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));;
	built = (built && createPortOut<double, GenericInput>(this, &mimmo::GenericInput::getResult<double>, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortOut<int, GenericInput>(this, &mimmo::GenericInput::getResult<int>, PortType::M_VALUEI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortOut<bool, GenericInput>(this, &mimmo::GenericInput::getResult<bool>, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortOut<iarray3E, GenericInput>(this, &mimmo::GenericInput::getResult<iarray3E>, PortType::M_DEG, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));
	built = (built && createPortOut<string, GenericInput>(this, &mimmo::GenericInput::getResult<string>, PortType::M_FILENAME, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING));
	built = (built && createPortOut<string, GenericInput>(this, &mimmo::GenericInput::getResult<string>, PortType::M_DIR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING));
// 	built = (built && createPortOut<string*, GenericInput>(this, &mimmo::GenericInput::getResult<string*>, PortType::M_FILENAMEPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
// 	built = (built && createPortOut<string*, GenericInput>(this, &mimmo::GenericInput::getResult<string*>, PortType::M_DIRPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
	
	m_arePortsBuilt = built;
}

/*!It clear the input member of the object
 */
void
GenericInput::clearInput(){
	m_input.reset(nullptr);
}

/*!It clear the result member of the object
 */
void
GenericInput::clearResult(){
	m_result.reset(nullptr);
}

/*!Execution command. It does nothing, the real execution of the object
 * happens in setInput/getResult.
 */
void
GenericInput::execute(){};


/*!
 * Get settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read only data in admissible format (see ports)
 * from a given file (all values in a row).
 * 
 * --> Absorbing data:
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>ReadFromFile</B>: 0/1 set class to read from a file	
 * - <B>Filename</B>: path to your current file data
 * 
 * \param[in] slotXML bitpit::Config::Section which reads from
 * \param[in] name   name associated to the slot
 */
void GenericInput::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	std::string input; 

	if(slotXML.hasOption("Priority")){
		input = slotXML.get("Priority");
		int value =0;
		if(!input.empty()){
			std::stringstream ss(bitpit::utils::trim(input));
			ss>>value;
		}
		setPriority(value);
	}; 
	
	if(slotXML.hasOption("ReadFromFile")){
		std::string input = slotXML.get("ReadFromFile");
		input = bitpit::utils::trim(input);
		bool temp = false;
		if(!input.empty()){
			std::stringstream ss(input);
			ss>>temp;
		}
		setReadFromFile(temp);
	}; 
	
	if(slotXML.hasOption("Filename") && m_readFromFile){
		std::string input = slotXML.get("Filename");
		input = bitpit::utils::trim(input);
		setFilename(input);
	}; 
	
}

/*!
 * Write settings of the class from bitpit::Config::Section slot. Reimplemented from
 * BaseManipulation::absorbSectionXML.The class read only data in admissible format (see ports)
 * from a given file (all values in a row).
 * 
 * --> Flushing data// how to write it on XML:
 * - <B>ClassName</B>: name of the class as "mimmo.GenericInput"
 * - <B>Priority</B>: uint marking priority in multi-chain execution; 
 * - <B>ReadFromFile</B>: 0/1 set class to read from a file	
 * - <B>Filename</B>: path to your current file data
 * 
 * \param[in] slotXML bitpit::Config::Section which writes to
 * \param[in] name   name associated to the slot
 */
void GenericInput::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){
	
	BITPIT_UNUSED(name);
	
	slotXML.set("ClassName", m_name);
	slotXML.set("Priority", std::to_string(getPriority()));
	
	slotXML.set("ReadFromFile", std::to_string((int)m_readFromFile));
	slotXML.set("Filename", m_filename);
};	

}

