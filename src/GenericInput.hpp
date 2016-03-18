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

#include <string>
#include "BaseManipulation.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief GenericInput is the class that set the initialization of a generic input data.
 *
 *	GenericInput is derived from BaseManipulation class. 
 *  Data passed as input are retained and always available in BaseManipulation::Result member.
 *
 */
class GenericInput: public BaseManipulation{
private:
	//members
	bool			m_readFromFile;	/**<True if the object reads the values from file.*/
	std::string		m_filename;		/**<Name of the input file. The file has to be an ascii text file.*/
	
public:
	GenericInput(bool readFromFile = false);
	GenericInput(std::string filename);

	template<typename T>
	GenericInput(T data){
		setInput<T>(data);
	}
	~GenericInput();

	GenericInput(const GenericInput & other);
	GenericInput & operator=(const GenericInput & other);

	void setReadFromFile(bool readFromFile);
	void setFilename(std::string filename);

	void 	execute();

};

#endif /* __INPUTDOF_HPP__ */
