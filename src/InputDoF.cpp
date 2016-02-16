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

#include "InputDoF.hpp"
#include <fstream>

using namespace std;

/*!Default constructor of OutputDoF.
 * \param[in] readFromFile True if the object reads the values from file.
 */
InputDoF::InputDoF(bool readFromFile){
	m_readFromFile = readFromFile;
};

/*!Custom constructor of OutputDoF.
 * \param[in] filename Name of the input file.
 * The m_readFromFile flag is set to true.
 */
InputDoF::InputDoF(string filename){
	m_readFromFile 	= true;
	m_filename 		= filename;
};

/*!Custom constructor of OutputDoF.
 * \param[in] ndeg #Degrees of freedom.
 * \param[in] displacements Displacements of the degrees of freedom.
 * The m_readFromFile flag is set to false.
 */
InputDoF::InputDoF(uint32_t ndeg, dvecarr3E & displacements){
	m_readFromFile = false;
	setNDeg(ndeg);
	setDisplacements(displacements);
};

InputDoF::~InputDoF(){};

/*!Copy constructor of InputDoF.
 */
InputDoF::InputDoF(const InputDoF & other):BaseManipulation(other){
	m_readFromFile 	= other.m_readFromFile;
	m_filename 		= other.m_filename;
};

/*!Assignement operator of InputDoF.
 */
InputDoF & InputDoF::operator=(const InputDoF & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_readFromFile 	= other.m_readFromFile;
	m_filename 		= other.m_filename;
	return *this;
};

/*!It sets if the object imports the displacements from an input file.
 * \param[in] readFromFile True if the object reads the values from file.
 */
void
InputDoF::setReadFromFile(bool readFromFile){
	m_readFromFile = readFromFile;
};

/*!It sets the name of the input file.
 * \param[in] filename Name of the input file.
 */
void
InputDoF::setFilename(std::string filename){
	m_filename = filename;
};

/*!Execution command. It reads by file the displacements of the degrees of freedom.
 */
void
InputDoF::exec(){
	if (m_readFromFile){
		ifstream file;
		file.open(m_filename);
		if (file.is_open()){
			file >> BaseManipulation::m_ndeg;
			darray3E displ;
			while(!file.eof()){
				for (int i=0; i<3; i++){
					file >> displ[i];
				}
				BaseManipulation::m_displ.push_back(displ);
			}
			file.close();
		}else{
			BaseManipulation::m_ndeg = 0;
			BaseManipulation::m_displ.clear();
		}
	}
};

