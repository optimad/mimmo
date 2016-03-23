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

/*!Default constructor of GenericOutput
 * \param[in] filename Name of the output file (default value = "output.txt").
 */
GenericOutput::GenericOutput(std::string filename){
	m_filename	= filename;
	m_pinType	= PinsType::BACKWARD;
	m_name 		= "MiMMO.GenericOutput";

};

/*!Default destructor of GenericOutput.
 */
GenericOutput::~GenericOutput(){
	m_pinType 	= PinsType::BACKWARD;
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

/*!It adds a name of the output files.
 * \param[in] filename Name of the output file.
 */
void
GenericOutput::setFilename(std::string filename){
	m_filename = filename;
};

/*!Execution command. VOID (setInput does the work).
 */
void
GenericOutput::execute(){};

