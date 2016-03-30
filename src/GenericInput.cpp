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

#include "GenericInput.hpp"
#include "Operators.hpp"
#include <fstream>

using namespace std;
using namespace mimmo;

/*!Default constructor of GenericInput.
 * \param[in] readFromFile True if the object reads the values from file (default value false).
 */
GenericInput::GenericInput(bool readFromFile){
	m_readFromFile = readFromFile;
	m_pinType 		= BaseManipulation::PinsType::FORWARD;
	m_name 			= "MiMMO.GenericInput";
};

/*!Custom constructor of GenericInput.
 * \param[in] filename Name of the input file.
 * The m_readFromFile flag is set to true.
 */
GenericInput::GenericInput(string filename){
	m_readFromFile 	= true;
	m_filename 		= filename;
	m_pinType 		= BaseManipulation::PinsType::FORWARD;
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

/*!Execution command. It does nothing, the real execution of the object
 * happens in setInput/getResult.
 */
void
GenericInput::execute(){};




