/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
#ifndef __OUTPUTDOF_HPP__
#define __OUTPUTDOF_HPP__

#include "BaseManipulation.hpp"
#include "IOData.hpp"
#include <string>

namespace mimmo{

/*!
 * \class GenericOutput
 * \ingroup iogeneric
 * \brief GenericOutput is the class that write data in an output on file.
 *
 * GenericOutput is derived from BaseManipulation class.
 * GenericOutput can write the output in a data file (unformatted/csv).
 * When writing, the GenericOutput object recognizes the type of data and
 * it adapts the writing method in function of the input port used
 * to build link (pin) with an other object.
 *
 * \n
 * Ports available in GenericOutput Class :
 *
 *	=========================================================

 	|                    Port Input  |||                                   |
 	|-------|-------------|-------------------|-------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 	| 0     | M_COORDS    | setResult         | (VECARR3, FLOAT)        |
 	| 10    | M_DISPLS    | setResult         | (VECARR3, FLOAT)        |
    | 14    | M_DATAFIELD | setResult         | (VECTOR, FLOAT)     |
    | 18    | M_SCALARFIELD | setResult         | (MPVECTOR, FLOAT)     |
    | 19    | M_VECTORFIELD | setResult         | (MPVECARR3, FLOAT)     |
 	| 20    | M_POINT     | setResult         | (ARRAY3, FLOAT)         |
 	| 24    | M_DIMENSION | setResult         | (ARRAY3, INT)           |
 	| 30    | M_VALUED    | setResult         | (SCALAR, FLOAT)         |
 	| 31    | M_VALUEI    | setResult         | (SCALAR, INT)           |
 	| 32    | M_VALUEB    | setResult         | (SCALAR, BOOL)          |
 	| 40    | M_DEG       | setResult         | (ARRAY3, INT)           |


 	|              Port Output  |||              |
 	|-------|-------------|---------|----------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |

 *	=========================================================
 *
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericOutput</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>Filename</B>: name of file to write data;
 * - <B>WriteDir</B>: name of directory to write data;
 * - <B>CSV</B>: true if write in csv format;
 * - <B>Binary</B>: 0/1 set write a BINARY file (default ASCII);
 *
 */
class GenericOutput: public BaseManipulation{
private:
    std::string             m_dir;          /**<Name of directory to write the output file.*/
    std::string             m_filename;     /**<Name of the output file.
                                            The file will be an ascii text file.*/

	bool                    m_csv;          /**<True if write output file in csv format.*/
	std::unique_ptr<IOData>	m_input;		/**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/
	std::unique_ptr<IOData>	m_result;		/**<Pointer to a base class object Result (derived class is template).*/

    bool            m_binary;       /**<Output binary files (used only for MimmoPiercedVector structures).*/

public:
	GenericOutput(std::string dir = "./", std::string filename = "output.txt", bool csv = false);
	GenericOutput(const bitpit::Config::Section & rootXML);
	~GenericOutput();

	GenericOutput(const GenericOutput & other);
    GenericOutput& operator=(GenericOutput other);
    
	void buildPorts();

	BITPIT_DEPRECATED( template<typename T>
	T*					getInput());

	BITPIT_DEPRECATED( template<typename T>
	T* 					getResult());

	template<typename T>
	void 	setInput(T* data);

	template<typename T>
	void 	setInput(T data);

	BITPIT_DEPRECATED( template<typename T>
	void 				setResult(T* data));
	BITPIT_DEPRECATED( template<typename T>
	void 				setResult(T data));

	void	clearInput();
	void	clearResult();

    void setWriteDir(std::string dir);
    void setFilename(std::string filename);
    void setCSV(bool csv);
    void setBinary(bool binary);

	void 	execute();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
protected:
    void swap(GenericOutput & x) noexcept;
    
private:

	template<typename T>
    void                _setInput(T* data);
    template<typename T>
    void                _setInput(T data);

    template<typename T>
    void                _setResult(T* data);
    template<typename T>
    void                _setResult(T data);

    template<typename T>
    T*                  _getInput();

    template<typename T>
    T*                  _getResult();

    template<typename T>
    std::ofstream&  ofstreamcsv(std::ofstream &in, const T &x);
    template<typename T>
    std::ofstream&  ofstreamcsv(std::ofstream &in, const std::vector< T > &x);
    template<typename T, size_t d>
    std::ofstream&  ofstreamcsv(std::ofstream &in, const std::array< T,d > &x);
    template<typename T>
    std::ofstream&  ofstreamcsvend(std::ofstream &in, const T &x);
    template<typename T>
    std::ofstream&  ofstreamcsvend(std::ofstream &in, const std::vector< T > &x);
    template<typename T, size_t d>
    std::ofstream&  ofstreamcsvend(std::ofstream &in, const std::array< T,d > &x);

};

template <>
void  GenericOutput::setResult(dmpvector1D data);

template <>
void  GenericOutput::setResult(dmpvecarr3E data);

REGISTER(BaseManipulation, GenericOutput, "mimmo.GenericOutput")
}

#include "GenericOutput.tpp"

#endif /* __OUTPUTDOF_HPP__ */
