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

#include "BaseManipulation.hpp"
#include <string>

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief GenericOutput is the class that write the degrees of freedom.
 *
 *	GenericOutput is derived from BaseManipulation class.
 *	It uses and it write the base members m_ndeg and m_displacements.
 *	After the execution of an object GenericOutput, the number of degrees of freedom and their initial
 *	displacements are write on a text file (ascii).
 */
class GenericOutput: public BaseManipulation{
private:
	//members
	std::string	m_filename;		/**<Name of the output file. The file will be an ascii text file.*/

public:
	GenericOutput(std::string filename = "output.txt");
	~GenericOutput();

	GenericOutput(const GenericOutput & other);
	GenericOutput & operator=(const GenericOutput & other);

	void setFilename(std::string filename);

	template<typename T>
	void 	setInput(T* data);

	template<typename T>
	void 	setInput(T& data);

	void 	execute();

};

#include "GenericOutput.tpp"

#endif /* __OUTPUTDOF_HPP__ */
