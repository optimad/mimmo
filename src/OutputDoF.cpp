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

#include "OutputDoF.hpp"
#include <fstream>

using namespace std;

/*!Default constructor of OutputDoF
 * \param[in] filename Name of the output file (default value = "output.txt").
 */
OutputDoF::OutputDoF(std::string filename){
	m_filename.push_back(filename);
};

/*!Default destructor of OutputDoF.
 */
OutputDoF::~OutputDoF(){};

/*!Copy constructor of OutputDoF.
 */
OutputDoF::OutputDoF(const OutputDoF & other):BaseManipulation(other){
	m_filename 		= other.m_filename;
};

/*!Assignement operator of OutputDoF.
 */
OutputDoF & OutputDoF::operator=(const OutputDoF & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_filename 		= other.m_filename;
	return *this;
};

/*!It adds a name of the output files.
 * \param[in] filename Name of the output file.
 */
void
OutputDoF::addFilename(std::string filename){
	m_filename.push_back(filename);
};

/*!Execution command. It writes on file the displacements of the degrees of freedom
 *of the parent manipulation.
 */
void
OutputDoF::execute(){
	for (int i=0; i<getNParent(); i++){
		if (m_parent[i] != NULL && i<m_filename.size()){
			m_ndeg = m_parent[i]->getNDeg();
			m_displ = m_parent[i]->getDisplacements();
			ofstream file;
			file.open(m_filename[i]);
			if (file.is_open()){
				file << BaseManipulation::m_ndeg << "\n";
				for (int iv=0; iv<m_ndeg; iv++){
					for (int i=0; i<3; i++){
						file << BaseManipulation::m_displ[iv][i] << "\t";
					}
					file << "\n";
				}
				file.close();
			}
		}
	}
};

