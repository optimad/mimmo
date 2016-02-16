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

/*!Default constructor of OutputDoF.
 * \param[in] parent Pointer to reference manipulator object to be linked (default value = NULL).
 */
OutputDoF::OutputDoF(BaseManipulation* parent):BaseManipulation(parent){};

/*!Custom constructor of OutputDoF
 * \param[in] filename Name of the output file.
 * \param[in] parent Pointer to reference manipulator object to be linked (default value = NULL).
 */
OutputDoF::OutputDoF(std::string filename, BaseManipulation* parent):BaseManipulation(parent){
	m_filename 		= filename;
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

/*!It sets the name of the output file.
 * \param[in] filename Name of the output file.
 */
void
OutputDoF::setFilename(std::string filename){
	m_filename = filename;
};

/*!Execution command. It writes on file the displacements of the degrees of freedom
 * given by the parent manipulation.
 */
void
OutputDoF::exec(){
	ofstream file;
	file.open(m_filename);
	if (file.is_open()){
		file << BaseManipulation::m_ndeg;
		for (int iv=0; iv<m_ndeg; iv++){
			for (int i=0; i<3; i++){
				file << BaseManipulation::m_displ[iv][i];
			}
		}
		file.close();
	}
};

