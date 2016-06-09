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
	built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, M_COORDS, {M_GDISPLS, M_DISPLS, M_GLOBAL, M_LOCAL}));
	built = (built && createPortIn<dvecarr3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvecarr3E>, M_DISPLS, {M_GDISPLS, M_COORDS}));
	built = (built && createPortIn<dvector1D, GenericOutput>(this, &mimmo::GenericOutput::setInput<dvector1D>, M_FILTER));
	built = (built && createPortIn<darray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<darray3E>, M_POINT, {M_POINT2}));
	built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, M_DIMENSION));
	built = (built && createPortIn<double, GenericOutput>(this, &mimmo::GenericOutput::setInput<double>, M_VALUED, {M_VALUED2}));
	built = (built && createPortIn<int, GenericOutput>(this, &mimmo::GenericOutput::setInput<int>, M_VALUEI));
	built = (built && createPortIn<bool, GenericOutput>(this, &mimmo::GenericOutput::setInput<bool>, M_VALUEB));
	built = (built && createPortIn<iarray3E, GenericOutput>(this, &mimmo::GenericOutput::setInput<iarray3E>, M_DEG));
	built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, M_FILENAME));
	built = (built && createPortIn<string, GenericOutput>(this, &mimmo::GenericOutput::setInput<string>, M_DIR));
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

