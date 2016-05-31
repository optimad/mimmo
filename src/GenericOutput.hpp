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

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief GenericOutput is the class that write the an output on file.
 *
 *	GenericOutput is derived from BaseManipulation class.
 *	It uses and it writes the base members input.
 *
 *	=========================================================
 *	PORT TABLE
 *
 *	- Input  Ports-
 *
 *	PortID	PortType	variable/function	compatibilities
 *
 *	0 		DISPLS		getResult			GDISPLS COORDS
 *	1 		COORDS		getResult			GDISPLS DISPLS
 *	2 		FILTER		getResult			-
 *
 *
 *	- Output Ports -
 *
 *	PortID	PortType	variable/function
 *
 *	=========================================================
 *
 */
class GenericOutput: public BaseManipulation{
private:
	//members
	std::string				m_filename;		/**<Name of the output file.
											The file will be an ascii text file.*/

	std::unique_ptr<IOData>	m_input;		/**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/
	std::unique_ptr<IOData>	m_result;		/**<Pointer to a base class object Result (derived class is template).*/

public:
	GenericOutput(std::string filename = "output.txt");
	~GenericOutput();

	GenericOutput(const GenericOutput & other);
	GenericOutput & operator=(const GenericOutput & other);

	void buildPorts();

	template<typename T>
	T*					getInput();

	template<typename T>
	T* 					getResult();

	template<typename T>
	void 	setInput(T* data);

	template<typename T>
	void 	setInput(T data);

	template<typename T>
	void 				setResult(T* data);
	template<typename T>
	void 				setResult(T data);

	template<typename T>
	void 				_setInput(T* data);
	template<typename T>
	void 				_setInput(T data);

	template<typename T>
	void 				_setResult(T* data);
	template<typename T>
	void 				_setResult(T data);

	template<typename T>
	T*					_getInput();

	template<typename T>
	T*					_getResult();

	void	clearInput();
	void	clearResult();

	void setFilename(std::string filename);

	void 	execute();

};

}

#include "GenericOutput.tpp"

#endif /* __OUTPUTDOF_HPP__ */
