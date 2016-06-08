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

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief GenericInput is the class that set the initialization of a generic input data.
 *
 *	GenericInput is derived from BaseManipulation class. 
 *  Data passed as input are retained and always available
 *  in BaseManipulation::Result member.
 *  GenericInput can be read the input from a file or
 *   it can be set by using setInput methods.
 *
 *	=========================================================
 * ~~~
 *	|--------------------------------------------------------------|
 *	|                 Port Input                                   |
 *	|-------|----------|-------------------|-----------------------|
 *	|PortID | PortType | variable/function | compatibilities       |
 *	|-------|----------|-------------------|-----------------------|
 *	|-------|----------|-------------------|-----------------------|
 *
 *
 *	|--------------------------------------|
 *	|            Port Output               |
 *	|-------|-----------|-------------------|
 *	|PortID | PortType  | variable/function |
 *	|-------|-----------|-------------------|
 *	| 0     | COORDS    | getResult         |
 *	| 10    | DISPLS    | getResult		    |
 *	| 12    | FILTER    | getResult    	    |
 *	| 20    | POINT     | getResult    	    |
 *	| 24    | DIMENSION | getResult    	    |
 *	|-------|-----------|-------------------|
 * ~~~
 *	=========================================================
 *
 */
class GenericInput: public BaseManipulation{
private:
	//members
	bool			m_readFromFile;	/**<True if the object reads the values from file.*/
	std::string		m_filename;		/**<Name of the input file. The file has to be an ascii text file.*/
	
	std::unique_ptr<IOData>				m_input;		/**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/
	std::unique_ptr<IOData>				m_result;		/**<Pointer to a base class object Result (derived class is template).*/

public:
	GenericInput(bool readFromFile = false);
	GenericInput(std::string filename);

	/*!Custom template constructor of GenericInput.
	 * It sets the base class input with data passed as argument.
	 * \param[in] data Data used to set the input.
	 */
	template<typename T>
	GenericInput(T data){
		setInput<T>(data);
	}
	~GenericInput();

	GenericInput(const GenericInput & other);
	GenericInput & operator=(const GenericInput & other);

	void buildPorts();

	template<typename T>
	T*					getInput();

	template<typename T>
	T 					getResult();

	void setReadFromFile(bool readFromFile);
	void setFilename(std::string filename);

	template<typename T>
	void 				setInput(T* data);
	template<typename T>
	void 				setInput(T& data);

	template<typename T>
	void 				setResult(T* data);
	template<typename T>
	void 				setResult(T& data);

	template<typename T>
	void 				_setInput(T* data);
	template<typename T>
	void 				_setInput(T& data);

	template<typename T>
	void 				_setResult(T* data);
	template<typename T>
	void 				_setResult(T& data);

	template<typename T>
	T*					_getInput();

	template<typename T>
	T*					_getResult();

	void	clearInput();
	void	clearResult();

	void execute();

};

}

#include "GenericInput.tpp"

#endif /* __INPUTDOF_HPP__ */
