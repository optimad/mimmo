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

#ifndef __INPUTDOF_HPP__
#define __INPUTDOF_HPP__

#include <string>
#include "BaseManipulation.hpp"
#include "IOData.hpp"

namespace mimmo{

/*!
 * \brief Utilities to read CSV data from an input stream.
 * \ingroup iogeneric
 */
namespace inputCSVStream{
    template<typename T>
    std::fstream&  ifstreamcsv(std::fstream &in, T &x);
    template<typename T>
    std::fstream&  ifstreamcsvend(std::fstream &in, T &x);
    template<typename T>
    std::fstream&  ifstreamcsv(std::fstream &in, std::vector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ifstreamcsv(std::fstream &in, std::array< T,d > &x);
    template<typename T>
    std::fstream&  ifstreamcsvend(std::fstream &in, std::vector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ifstreamcsvend(std::fstream &in, std::array< T,d > &x);
    template<typename T>
    std::fstream&  ifstreamcsv(std::fstream &in, MimmoPiercedVector< T > &x);
    template<typename T, size_t d>
    std::fstream&  ifstreamcsv(std::fstream &in, MimmoPiercedVector< std::array< T,d > > &x);

}

/*!
 * \class GenericInput
 * \ingroup iogeneric
 * \brief GenericInput is the class that set the initialization of a generic input data.
 *
 * GenericInput is derived from BaseManipulation class.
 * GenericInput can read the input from a data file (unformatted/csv) or
 * it can be set by using setInput methods.
 * When reading, the GenericInput object recognizes the type of data and
 * it adapts the reading method in function of the output port used
 * to build link (pin) with other objects.
 * On MPI versions, GenericInput read data from file with proc rank 0, and send
 * all the data to all other proc.
 *
 * \n
 * Ports available in GenericInput Class :
 *
 *    =========================================================

     |                 Port Input  ||                                 |
     |----------|-------------------|-----------------------|
     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     |              Port Output   |               |                       |
     |---------------|-------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_COORDS      | getResult         | (MC_VECARR3, MD_FLOAT)      |
     | M_DISPLS      | getResult         | (MC_VECARR3, MD_FLOAT)      |
     | M_DATAFIELD   | getResult         | (MC_VECTOR, MD_FLOAT)       |
     | M_POINT       | getResult         | (MC_ARRAY3, MD_FLOAT)       |
     | M_SPAN        | getResult         | (MC_ARRAY3, MD_FLOAT)       |
     | M_DIMENSION   | getResult         | (MC_ARRAY3, MD_INT)         |
     | M_VALUED      | getResult         | (MC_SCALAR, MD_FLOAT)       |
     | M_VALUEI      | getResult         | (MC_SCALAR, MD_INT)         |
     | M_VALUEB      | getResult         | (MC_SCALAR, MD_BOOL)        |
     | M_DEG         | getResult         | (MC_ARRAY3, MD_INT)         |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericInput</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>ReadFromFile</B>: 0/1 set class to read from a file;
 * - <B>CSV</B>: 0/1 set class to read a CSV format;
 * - <B>ReadDir</B>: path to your current file data;
 * - <B>Filename</B>: path to your current file data.
 *
 */
class GenericInput: public BaseManipulation{

// template<typename T>
// friend bitpit::OBinaryStream & (::operator<<)(bitpit::OBinaryStream &buffer, const std::array<T,3> &array3);
// template<typename T>
// friend bitpit::OBinaryStream & (::operator<<)(bitpit::OBinaryStream &buffer, const std::vector<T> &vector);
// template<typename T>
// friend bitpit::OBinaryStream & (::operator<<)(bitpit::OBinaryStream &buffer, const std::vector<std::array<T,3>> &vecarray3);


private:
    bool            m_readFromFile; /**<True if the object reads the values from file.*/
    bool            m_csv;          /**<True if the file is in csv format.*/
    std::string     m_dir;          /**<Name of directory to read the input file.*/
    std::string     m_filename;     /**<Name of the input file. The file has to be an ascii text file.*/

    std::unique_ptr<IOData>                m_input;        /**<Pointer to a base class object Input, meant for input temporary data, cleanable in execution (derived class is template).*/
    std::unique_ptr<IOData>                m_result;        /**<Pointer to a base class object Result (derived class is template).*/

public:
    GenericInput(bool readFromFile = false, bool csv = false);
    GenericInput(const bitpit::Config::Section & rootXML);
    GenericInput(std::string dir, std::string filename, bool csv = false);

    /*!Custom template constructor of GenericInput.
     * It sets the base class input with data passed as argument.
     * \param[in] data Data used to set the input.
     */
    BITPIT_DEPRECATED( template<typename T>
    GenericInput(T data){
        setInput<T>(data);
    });
    ~GenericInput();

    GenericInput(const GenericInput & other);
    GenericInput & operator=(GenericInput other);

    void buildPorts();

    BITPIT_DEPRECATED( template<typename T> T* getInput());

    template<typename T>
    T  getResult();

    void setReadFromFile(bool readFromFile);
    void setCSV(bool csv);
    void setReadDir(std::string dir);
    void setFilename(std::string filename);

    template<typename T>
    void                 setInput(T* data);
    template<typename T>
    void                 setInput(T& data);

    BITPIT_DEPRECATED( template<typename T>
    void                 setResult(T* data));
    BITPIT_DEPRECATED( template<typename T>
    void                 setResult(T& data));

    void    clearInput();
    void    clearResult();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(GenericInput & x) noexcept;
#if MIMMO_ENABLE_MPI
    template<typename T>
    void                sendReadDataToAllProcs(T & data);
#endif

private:

    template<typename T>
    void                 _setInput(T* data);
    template<typename T>
    void                 _setInput(T& data);

    template<typename T>
    void                 _setResult(T* data);
    template<typename T>
    void                 _setResult(T& data);

    template<typename T>
    T*                    _getInput();

    template<typename T>
    T*                    _getResult();

};

/*!
 * \class GenericInputMPVData
 * \ingroup iogeneric
 * \brief GenericInputMPVData is the class that set a generic input data as mimmo::MimmoPiercedVector.
 *
 * GenericInputMPVData is derived from BaseManipulation class.
 * GenericInput read MimmoPiercedVector data fields from file (unformatted/csv)
 * When reading, the GenericInput object recognizes the type of data and
 * it adapts the reading method in function of the output port used
 * to build link (pin) with other objects.
 * If a valid reference geometry is linked through port M_GEOM or through method setGeometry
 * the class checks all MPV entries read  and retains those whose id is contained inside the
 * reference geometry point-data, cell-data or interface-data.
 * If none of the items is referenced by the geometry return an empty MimmoPiercedVector as
 * result.
 * Can read binary data related to PiercedVector, activating the binary flag.
 * On MPI versions, GenericInputMPVData read data from file with proc rank 0, and send
 * all the data to all other procs.
 *
 * \n
 * Ports available in GenericInput Class :
 *
 *    =========================================================
 *
 *     |                 Port Input  ||                                 |
 *     |----------|-------------------|-----------------------|
 *     |<B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 *     |M_GEOM   |  setGeometry | (MC_SCALAR, MD_MIMMO_) |
 *
 *     |              Port Output   |               |                       |
 *     |---------------|-------------------|-----------------------|
 *     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
 *     | M_SCALARFIELD | getResult         | (MC_MPVECTOR, MD_FLOAT)     |
 *     | M_VECTORFIELD | getResult         | (MC_MPVECARR3, MD_FLOAT)    |
 *
 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericInputMPVData</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 *
 * Proper of the class:
 * - <B>CSV</B>: 0/1 set class to read a CSV format;
 * - <B>ReadDir</B>: path to your current file data;
 * - <B>Filename</B>: path to your current file data.
 * - <B>Binary</B>: 0/1 set read a BINARY file (default ASCII);
 *
 */
class GenericInputMPVData: public BaseManipulation{
private:
    bool            m_csv;          /**<True if the file is in csv format.*/
    std::string     m_dir;          /**<Name of directory to read the input file.*/
    std::string     m_filename;     /**<Name of the input file. */

    std::unique_ptr<IOData>                m_result;        /**<Pointer to a base class object Result (derived class is template).*/

    bool            m_binary;       /**<Input binary files (used only for MimmoPiercedVector structures).*/

public:
    GenericInputMPVData(bool csv = false);
    GenericInputMPVData(const bitpit::Config::Section & rootXML);
    GenericInputMPVData(std::string dir, std::string filename, bool csv = false);

    ~GenericInputMPVData();

    GenericInputMPVData(const GenericInputMPVData & other);
    GenericInputMPVData & operator=(GenericInputMPVData other);

    void buildPorts();

    template<typename T>
    MimmoPiercedVector<T>  getResult();

    void setCSV(bool csv);
    void setReadDir(std::string dir);
    void setFilename(std::string filename);
    void setBinary(bool binary);

    void    clearResult();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void swap(GenericInputMPVData & x) noexcept;
#if MIMMO_ENABLE_MPI
    template<typename T>
    void sendReadDataToAllProcs(MimmoPiercedVector<T> & data);
#endif

private:

    template<typename T>
    MimmoPiercedVector< T >*         _getResult();
    template<typename T>
    void                             _setResult(MimmoPiercedVector< T >*);
    template<typename T>
    void                             _setResult(MimmoPiercedVector< T >&);


};

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_DATAFIELD, MC_VECTOR, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_SPAN, MC_ARRAY3, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_DIMENSION, MC_ARRAY3, MD_INT,__INPUTDOF_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_VALUEI, MC_SCALAR, MD_INT,__INPUTDOF_HPP__)
REGISTER_PORT(M_VALUEB, MC_SCALAR, MD_BOOL,__INPUTDOF_HPP__)
REGISTER_PORT(M_DEG, MC_ARRAY3, MD_INT,__INPUTDOF_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_MPVECARR3, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_MPVECTOR, MD_FLOAT,__INPUTDOF_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__INPUTDOF_HPP__)

REGISTER(BaseManipulation, GenericInput, "mimmo.GenericInput")
REGISTER(BaseManipulation, GenericInputMPVData, "mimmo.GenericInputMPVData")

}

#include "GenericInput.tpp"

#endif /* __INPUTDOF_HPP__ */
