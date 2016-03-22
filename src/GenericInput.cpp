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
#include <fstream>

using namespace std;

/*!Default constructor of OutputDoF.
 * \param[in] readFromFile True if the object reads the values from file.
 */
GenericInput::GenericInput(bool readFromFile){
	m_readFromFile = readFromFile;
	m_pinType 		= BaseManipulation::PinsType::FORWARD;
};

/*!Custom constructor of OutputDoF.
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

/*!Execution command. It reads by file the input or it uses the input of base class
 * set by user.
 */
void
GenericInput::execute(){
	if (m_readFromFile){
		ifstream file;
		file.open(m_filename);
		if (file.is_open()){
			darray3E row;
			dvecarr3E input;
			while(!file.eof()){
				for (int i=0; i<3; i++){
					file >> row[i];
				}
				input.push_back(row);
			}
			file.close();
			setInput(input);
		}else{
			clearInput();
			cout << "file not open --> exit" << endl;
			exit(1);
		}
	}
	return;
};




