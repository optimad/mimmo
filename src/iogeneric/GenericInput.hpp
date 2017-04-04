/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#ifndef __INPUTDOF_HPP__
#define __INPUTDOF_HPP__

#include <string>
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\class GenericInput
 *	\brief GenericInput is the class that set the initialization of a generic input data.
 *
 * GenericInput is derived from BaseManipulation class. 
 * Data passed as input are retained and always available
 * in BaseManipulation::Result member.
 * GenericInput can be read the input from a file or
 * it can be set by using setInput methods.
 * 
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
 *	|-------------------------------------------|-----------------------|
 *	|              Port Output                  |                       |
 *	|-------|---------------|-------------------|-----------------------|
 *	|PortID | PortType      | variable/function | DataType              |
 *	|-------|---------------|-------------------|-----------------------|
 *	| 0     | M_COORDS      | getResult         | (VECARR3, FLOAT)      |
 *	| 10    | M_DISPLS      | getResult         | (VECARR3, FLOAT)      |
 *	| 12    | M_FILTER      | getResult         | (VECTOR, FLOAT)       |
 *	| 19    | M_SCALARFIELD | getResult         | (VECTOR, FLOAT)       |
 *	| 20    | M_POINT       | getResult         | (ARRAY3, FLOAT)       |
 *	| 23    | M_SPAN        | getResult         | (ARRAY3, FLOAT)       |
 *	| 24    | M_DIMENSION   | getResult         | (ARRAY3, INT)         |
 *	| 30    | M_VALUED      | getResult         | (SCALAR, FLOAT)       |
 *	| 31    | M_VALUEI      | getResult         | (SCALAR, INT)         |
 *	| 32    | M_VALUEB      | getResult         | (SCALAR, BOOL)        |
 *	| 40    | M_DEG         | getResult         | (ARRAY3, INT)         |
 *	| 50    | M_FILENAME    | getResult         | (SCALAR, STRING)      |
 *	| 51    | M_DIR         | getResult         | (SCALAR, STRING)      |
 *	|-------|---------------|-------------------|-----------------------|
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
	GenericInput(const bitpit::Config::Section & rootXML);
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
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

};

REGISTER(BaseManipulation, GenericInput, "mimmo.GenericInput")
}

#include "GenericInput.tpp"

#endif /* __INPUTDOF_HPP__ */
