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
 * 
 * \n
 * Ports available in GenericInput Class :
 *
 *    =========================================================

     |                 Port Input  |||                                 |
     |-------|----------|-------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |


     |              Port Output   ||               |                       |
     |-------|---------------|-------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 0     | M_COORDS      | getResult         | (VECARR3, FLOAT)      |
     | 10    | M_DISPLS      | getResult         | (VECARR3, FLOAT)      |
     | 12    | M_FILTER      | getResult         | (VECTOR, FLOAT)       |
     | 19    | M_SCALARFIELD | getResult         | (VECTOR, FLOAT)       |
     | 20    | M_POINT       | getResult         | (ARRAY3, FLOAT)       |
     | 23    | M_SPAN        | getResult         | (ARRAY3, FLOAT)       |
     | 24    | M_DIMENSION   | getResult         | (ARRAY3, INT)         |
     | 30    | M_VALUED      | getResult         | (SCALAR, FLOAT)       |
     | 31    | M_VALUEI      | getResult         | (SCALAR, INT)         |
     | 32    | M_VALUEB      | getResult         | (SCALAR, BOOL)        |
     | 40    | M_DEG         | getResult         | (ARRAY3, INT)         |
     | 50    | M_FILENAME    | getResult         | (SCALAR, STRING)      |
     | 51    | M_DIR         | getResult         | (SCALAR, STRING)      |

 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericInput</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>ReadFromFile</B>: 0/1 set class to read from a file;
 * - <B>CSV</B>: 0/1 set class to read a CSV format;
 * - <B>ReadDir</B>: path to your current file data;
 * - <B>Filename</B>: path to your current file data.
 *
 */
class GenericInput: public BaseManipulation{
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
    GenericInput & operator=(const GenericInput & other);

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

    template<typename T>
    std::ifstream&  ifstreamcsv(std::ifstream &in, T &x);
    template<typename T>
    std::ifstream&  ifstreamcsvend(std::ifstream &in, T &x);
    template<typename T>
    std::ifstream&  ifstreamcsv(std::ifstream &in, std::vector< T > &x);
    template<typename T, size_t d>
    std::ifstream&  ifstreamcsv(std::ifstream &in, std::array< T,d > &x);

};

REGISTER(BaseManipulation, GenericInput, "mimmo.GenericInput")
}

#include "GenericInput.tpp"

#endif /* __INPUTDOF_HPP__ */
