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

#include "GenericOutput.hpp"

using namespace std;
using namespace mimmo;

/*!Default constructor of GenericOutput
 * \param[in] filename Name of the output file (default value = "output.txt").
 */
GenericOutput::GenericOutput(std::string filename){
	m_filename	= filename;
	m_portsType	= ConnectionType::BACKWARD;
	m_name 		= "MiMMO.GenericOutput";

};

/*!Default destructor of GenericOutput.
 */
GenericOutput::~GenericOutput(){
	m_portsType = ConnectionType::BACKWARD;
};

/*!Copy constructor of GenericOutput.
 */
GenericOutput::GenericOutput(const GenericOutput & other):BaseManipulation(other){
	m_filename 		= other.m_filename;
};

/*!Assignement operator of GenericOutput.
 */
GenericOutput & GenericOutput::operator=(const GenericOutput & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_filename 		= other.m_filename;
	return *this;
};

/*! It builds the input/output ports of the object
 */
void
GenericOutput::buildPorts(){
	bool built = true;
	built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, PortType::M_COORDS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, PortType::M_DISPLS, mimmo::pin::containerTAG::VECARR3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<dvector1D, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvector1D>, PortType::M_FILTER, mimmo::pin::containerTAG::VECTOR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<darray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<darray3E>, PortType::M_POINT, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, PortType::M_DIMENSION, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<double, GenericOutput>(this, &mimmo::GenericOutput::setInput<double>, PortType::M_VALUED, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::FLOAT));
	built = (built && createPortIn<int, GenericOutput>(this, &mimmo::GenericOutput::setInput<int>, PortType::M_VALUEI, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<bool, GenericOutput>(this, &mimmo::GenericOutput::setInput<bool>, PortType::M_VALUEB, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::BOOL));
	built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, PortType::M_DEG, mimmo::pin::containerTAG::ARRAY3, mimmo::pin::dataTAG::INT));
	built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, PortType::M_FILENAME, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING));
	built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, PortType::M_DIR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING));
// 	built = (built && createPortIn<string*, GenericOutput>(this, &mimmo::GenericOutput::setInput<string*>, PortType::M_FILENAMEPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
// 	built = (built && createPortIn<string*, GenericOutput>(this, &mimmo::GenericOutput::setInput<string*>, PortType::M_DIRPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
	
	m_arePortsBuilt = built;
}

/*!It sets the name of the output file.
 * \param[in] filename Name of the output file.
 */
void
GenericOutput::setFilename(std::string filename){
	m_filename = filename;
};

/*!It clear the input member of the object
 */
void
GenericOutput::clearInput(){
	m_input.reset(nullptr);
}

/*!It clear the result member of the object
 */
void
GenericOutput::clearResult(){
	m_result.reset(nullptr);
}

/*!Execution command. It does nothing, the real execution of the object
 * happens in setInput.
 */
void
GenericOutput::execute(){
	rename("output.txt", m_filename.c_str());
};

