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
#ifndef __OUTPUTDOF_HPP__
#define __OUTPUTDOF_HPP__

#include "OutputDoF.hpp"
#include <ofstream>

using namespace std;

OutputDoF::OutputDoF(BaseManipulation* parent);
OutputDoF::OutputDoF(std::string filename, BaseManipulation* parent){
	m_filename 		= filename;

};

OutputDoF::~OutputDoF(){};

/*!Copy constructor of OutputDoF.
 */
OutputDoF::OutputDoF(const OutputDoF & other):BaseManipulation(other){
	m_filename 		= other.m_filename;
};

/*!Assignement operator of OutputDoF.
 */
OutputDoF & OutputDoF::operator=(const OutputDoF & other):BaseManipulation(other){
	m_filename 		= other.m_filename;
};

void
OutputDoF::setFilename(std::string filename){
	m_filename = filename;
};

void
OutputDoF::recoverDisplacements(){
	if (m_parent == NULL) return;
	setNDeg(m_parent->getNDeg());
	setDisplacements(m_parent->getDisplacements());
};

void
OutputDoF::exec(){
	recoverDisplacements();
	ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << BaseManipulation::m_ndeg;
		for (int iv=0; iv<m_ndeg; iv++){
			for (int i=0; i<3; i++){
				file << BaseManipulation::m_displacements[iv][i];
			}
		}
		file.close();
	}
};

