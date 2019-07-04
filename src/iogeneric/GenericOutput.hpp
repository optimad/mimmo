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
 *\brief Utilities to write CSV data to an output stream
 *\ingroup iogeneric
 */
namespace outputCSVStream{
    template<typename T>
    std::fstream&  ofstreamcsv(std::fstream &in, const T &x);
    template<typename T>
    std::fstream&  ofstreamcsvend(std::fstream &in, const T &x);
    template<typename T>
    std::fstream&  ofstreamcsv(std::fstream &in, const std::vector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ofstreamcsv(std::fstream &in, const std::array< T,d > &x);
    template<typename T>
    std::fstream&  ofstreamcsvend(std::fstream &in, const std::vector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ofstreamcsvend(std::fstream &in, const std::array< T,d > &x);
    template<typename T>
    std::fstream&  ofstreamcsv(std::fstream &in, const MimmoPiercedVector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ofstreamcsvend(std::fstream &in, const MimmoPiercedVector< T > &x);
}
/*!
 * \class GenericOutput
 * \ingroup iogeneric
 * \brief GenericOutput is the class that write generic data in a file output.
 *
 * GenericOutput is derived from BaseManipulation class.
 * GenericOutput can write the output in a data file (unformatted/csv).
 * When writing, the GenericOutput object recognizes the type of data and
 * it adapts the writing method in function of the input port used
 * to build link (pin) with an other object.
 *
 * On distributed archs, only the 0 rank procs is deputed to writing.
 * Since the class does not hold data structure that can be unambigously recollected
 * once partitioned, it is supposed that all procs hold the same identical data.
 *
 * \n
 * Ports available in GenericOutput Class :
 *
 *	=========================================================

 	|                    Port Input  ||                                |
 	|-------------|-------------------|-------------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_COORDS    | setResult         | (MC_VECARR3, MD_FLOAT)     |
 	| M_DISPLS    | setResult         | (MC_VECARR3, MD_FLOAT)     |
 	| M_DATAFIELD | setResult         | (MC_VECTOR, MD_FLOAT)      |
    | M_POINT     | setResult         | (MC_ARRAY3, MD_FLOAT)         |
    | M_SPAN      | setResult         | (MC_ARRAY3, MD_FLOAT)         |
 	| M_DIMENSION | setResult         | (MC_ARRAY3, MD_INT)           |
 	| M_VALUED    | setResult         | (MC_SCALAR, MD_FLOAT)         |
 	| M_VALUEI    | setResult         | (MC_SCALAR, MD_INT)           |
 	| M_VALUEB    | setResult         | (MC_SCALAR, MD_BOOL)          |
 	| M_DEG       | setResult         | (MC_ARRAY3, MD_INT)           |


 	|              Port Output  ||              |
 	|-------------|---------|----------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


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
 *
 */
class GenericOutput: public BaseManipulation{
private:
    std::string             m_dir;          /**<Name of directory to write the output file.*/
    std::string             m_filename;     /**<Name of the output file.
                                            The file will be an ascii text file.*/

	bool                    m_csv;          /**<True if write output file in csv format.*/
	std::unique_ptr<IOData>	m_input;		/**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/

public:
	GenericOutput(std::string dir = "./", std::string filename = "output.txt", bool csv = false);
	GenericOutput(const bitpit::Config::Section & rootXML);
	~GenericOutput();

	GenericOutput(const GenericOutput & other);
    GenericOutput& operator=(GenericOutput other);

	void buildPorts();

	BITPIT_DEPRECATED( template<typename T>
	T*      getInput());

	template<typename T>
	void 	setInput(T* data);

	template<typename T>
	void 	setInput(T data);

	void	clearInput();

    void setWriteDir(std::string dir);
    void setFilename(std::string filename);
    void setCSV(bool csv);

	void 	execute();

	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(GenericOutput & x) noexcept;

private:
    template<typename T>
    void _setInput(T & data);
};

/*!
 * \class GenericOutputMPVData
 * \ingroup iogeneric
 * \brief GenericOutputMPVData is the class that write a generic data to file output as
 * mimmo::MimmoPiercedVector.
 *
 * GenericOutputMPVData is derived from BaseManipulation class.
 * GenericOutputMPVData write MimmoPiercedVector data fields to file (unformatted/csv).
 * When writing, the GenericOutputMPVData object recognizes the type of data and
 * it adapts the writing method in function of the input port used
 * to build link (pin) with an other object.
 * Can write binary data in raw format only, activating the binary flag.
 *
 * On distributed archs, only the 0 rank procs is deputed to writing.
 * \n
 * Ports available in GenericOutputMPVData Class :
 *
 *  =========================================================
 *
 *  |                    Port Input  ||                                |
 *  |-------------|-------------------|-------------------------|
 *  | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 *  | M_SCALARFIELD | setInput       | (MC_MPVECTOR, MD_FLOAT)    |
 *  | M_VECTORFIELD | setInput       | (MC_MPVECARR3, MD_FLOAT)   |
 *
 *
 *  |              Port Output  ||              |
 *  |-------------|---------|----------|
 *  | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 *
 *
 *  =========================================================
 *
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericOutputMPVData</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>Filename</B>: name of file to write data;
 * - <B>WriteDir</B>: name of directory to write data;
 * - <B>CSV</B>: true if write in csv format;
 * - <B>Binary</B>: 0/1 set write a BINARY file (default ASCII).If CSV, only ASCII is available;
 *
 */
class GenericOutputMPVData: public BaseManipulation{
private:
    std::string             m_dir;          /**<Name of directory to write the output file.*/
    std::string             m_filename;     /**<Name of the output file.
    The file will be an ascii text file.*/

    bool                    m_csv;          /**<True if write output file in csv format.*/
    std::unique_ptr<IOData> m_input;        /**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/

    bool            m_binary;       /**<Output unformatted binary files.*/

public:
    GenericOutputMPVData(std::string dir = "./", std::string filename = "output.txt", bool csv = false);
    GenericOutputMPVData(const bitpit::Config::Section & rootXML);
    ~GenericOutputMPVData();

    GenericOutputMPVData(const GenericOutputMPVData & other);
    GenericOutputMPVData& operator=(GenericOutputMPVData other);

    void buildPorts();

    BITPIT_DEPRECATED(template<typename T>
    MimmoPiercedVector< T >*    getInput());

    template<typename T>
    void    setInput(MimmoPiercedVector< T > *data);

    template<typename T>
    void    setInput(MimmoPiercedVector< T > data);

    void    clearInput();

    void setWriteDir(std::string dir);
    void setFilename(std::string filename);
    void setCSV(bool csv);
    void setBinary(bool binary);

    void    execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(GenericOutputMPVData & x) noexcept;
#if MIMMO_ENABLE_MPI
    template<typename T>
    void    collectDataFromAllProcs(MimmoPiercedVector<T> & locdata, MimmoPiercedVector<T> * dataglobal);
#endif

private:
    template<typename T>
    void _setInput(MimmoPiercedVector< T > & data);

};

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_DATAFIELD, MC_VECTOR, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_MPVECARR3, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_MPVECTOR, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_DIMENSION, MC_ARRAY3, MD_INT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_VALUEI, MC_SCALAR, MD_INT,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_VALUEB, MC_SCALAR, MD_BOOL,__OUTPUTDOF_HPP__)
REGISTER_PORT(M_DEG, MC_ARRAY3, MD_INT,__OUTPUTDOF_HPP__)

REGISTER(BaseManipulation, GenericOutput, "mimmo.GenericOutput")
REGISTER(BaseManipulation, GenericOutputMPVData, "mimmo.GenericOutputMPVData")

}

#include "GenericOutput.tpp"

#endif /* __OUTPUTDOF_HPP__ */
