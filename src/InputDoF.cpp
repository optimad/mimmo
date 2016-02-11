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
#ifndef __INPUTDOF_HPP__
#define __INPUTDOF_HPP__

#include "InputDoF.hpp"
#include <ifstream>

using namespace std;

InputDoF::InputDoF(bool readFromFile){
	m_readFromFile = readFromFile;
};

InputDoF::InputDoF(string filename){
	m_readFromFile 	= true;
	m_filename 		= filename;
};

InputDoF::InputDoF(uint32_t ndeg, dvecarr3E & displacements):BaseManipulation(ndeg, displacements){
	m_readFromFile = false;
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
InputDoF & InputDoF::operator=(const InputDoF & other):BaseManipulation(other){
	m_readFromFile 	= other.m_readFromFile;
	m_filename 		= other.m_filename;
};

void
InputDoF::setReadFromFile(bool readFromFile){
	m_readFromFile = readFromFile;
};

void
InputDoF::setFilename(std::string filename){
	m_filename = filename;
};


void
InputDoF::recoverDisplacements(){
	if (readFromFile){
		ifstream file;
		file.open(m_filename);
		if (file.is_open()){
			file >> BaseManipulation::m_ndeg;
			darray3E displ;
			while(!file.eof()){
				for (int i=0; i<3; i++){
					file >> displ[i];
				}
				BaseManipulation::m_displacements.push_back(displ);
			}
			file.close();
		}else{
			BaseManipulation::m_ndeg = 0;
			BaseManipulation::m_displacements.clear();
		}
	}
};

void
InputDoF::exec(){
	recoverDisplacements();
};

