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

#include <fstream>
#include "Operators.hpp"

namespace mimmo{

    
namespace inputCSVStream{
    
/*!
 * Recover a data from a stream when import in csv format.
 * \param[in] in import stream.
 * \param[out] x data read.
 * \return Reference to result import stream after data import.
 */
template <typename T>
std::ifstream&
inputCSVStream::ifstreamcsv(std::ifstream &in, T &x){
    T                   dummy{};
    char                delim;
    if ((in.good())) {
        if (in >> dummy && in >> delim) { x = dummy;}
    }
    return(in);
}

/*!
 * Recover a data from a stream when import in csv format when end of line is reached.
 * \param[in] in import stream.
 * \param[out] x data read.
 * \return Reference to result import stream after data import.
 */
template <typename T>
std::ifstream&  inputCSVStream::ifstreamcsvend(std::ifstream &in, T &x){
    T                   dummy{};
    if ((in.good())) {
        if (in >> dummy) { x = dummy;}
    }
    return(in);
}

/*!
 * Recover a data from a stream when import in csv format a vector of data.
 * \param[in] in import stream.
 * \param[out] x vector data read.
 * \return Reference to result import stream after data import.
 */
template <typename T>
std::ifstream&  inputCSVStream::ifstreamcsv(std::ifstream &in, std::vector< T > &x){

    T       dummy;

    while (in.good()) {
        if (ifstreamcsv(in,dummy)) { x.push_back(dummy); }
    }
    return(in);
}

/*!
 * Recover a data from a stream when import in csv format an array of data of dimension d.
 * \param[in] in import stream.
 * \param[out] x array data read.
 * \return Reference to result import stream after data import.
 */
template <typename T, size_t d>
std::ifstream&  inputCSVStream::ifstreamcsv(std::ifstream &in, std::array< T, d> &x){

    T       dummy{};
    int     i;

    i = 0;
    while ((in.good()) && (i < d-1)) {
        if (ifstreamcsv(in,dummy)) { x[i] = dummy;}
        i++;
    } //next i
    if (ifstreamcsvend(in,dummy)) { x[i] = dummy;}
    return(in);
};

/*!
 * Recover a data from a stream when import in csv format a MimmoPiercedVector of data.
 * \param[in] in import stream.
 * \param[out] x MimmoPiercedVector data read.
 * \return Reference to result import stream after data import.
 */
template <typename T>
std::ifstream&  inputCSVStream::ifstreamcsv(std::ifstream &in, MimmoPiercedVector< T > &x){

    T       dummy;
    long    id;
    int     location;

    if(istreamcsv(in, location)){
        if(x.intIsValidLocation(location)) x.setDataLocation(static_cast<MPVLocation>(location));
    }
    while (in.good()) {
        if (ifstreamcsv(in,id) && ifstreamcsv(in,dummy)) { x.insert(id,dummy); }
    }
    return(in);
}

}

///GENERICINPUT////////////////////////////////////////////////////////////////////////////

/*!
 * Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Pointer to data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T* data){
    _setInput(data);
    _setResult(data);
}

/*!
 * Overloaded function of base class setInput.
 * It sets the input of the object, but at the same time it sets even the result.
 * \param[in] data Data to be used to set the input/result.
 */
template<typename T>
void
GenericInput::setInput(T& data){
    _setInput(data);
    _setResult(data);
}

/*!
 * Overloaded function of base class getResult.
 * It gets the result of the object, equal to the input.
 * In the case it reads the input from file before to set and to get the result.
 * \return Pointer to data stored in result member.
 */
template<typename T>
T
GenericInput::getResult(){
    if (m_readFromFile){
        T data;
        std::ifstream file;
        file.open(m_dir+"/"+m_filename);
        if (file.is_open()){
            if (m_csv){
                ifstreamcsv(file, data);
            }
            else{
                file >> data;
            }
            file.close();
        }else{
            (*m_log) << "file not open --> exit" << std::endl;
            throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
        }
        _setResult(data);
    }

    return  *static_cast<IODataT<T>*>(m_result.get())->getData();
}


/*!
 * It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericInput::getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

/*!
 * It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInput::setResult(T* data){
    m_result = std::move(std::unique_ptr<IOData>(new IODataT<T>(*data)));
}

/*!It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInput::setResult(T& data){
    m_result = std::move(std::unique_ptr<IOData>(new IODataT<T>(data)));
}

/*!
 * It sets the input member of the object.
 * \param[in] data Pointer to data to be stored in the input member.
 */
template<typename T>
void
GenericInput::_setInput(T* data){
    m_input = std::move(std::unique_ptr<IOData>(new IODataT<T>(*data)));
}

/*!
 * It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
GenericInput::_setInput(T& data){
    m_input = std::move(std::unique_ptr<IOData>(new IODataT<T>(data)));
}

/*!
 * It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericInput::_getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

/*!
 * It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInput::_setResult(T* data){
    m_result = std::move(std::unique_ptr<IOData> (new IODataT<T>(*data)));
}

/*!
 * It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInput::_setResult(T& data){
    m_result = std::move(std::unique_ptr<IOData> (new IODataT<T>(data)));
}

/*!
 * It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
GenericInput::_getResult(){
    return(static_cast<IODataT<T>*>(m_result.get())->getData());
}

//GENERICINPUTMPVDATA/////////////////////////////////////////////////////

/*!
 * It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericInputMPVData::getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}


/*!
 * It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInputMPVData::_setResult(T* data){
    m_result = std::move(std::unique_ptr<IOData>(new IODataT<T>(*data)));
}

/*!
 * It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInputMPVData::_setResult(T& data){
    m_result = std::move(std::unique_ptr<IOData>(new IODataT<T>(data)));
}

/*!
 * It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
GenericInputMPVData::_getResult(){
    return(static_cast<IODataT<T>*>(m_result.get())->getData());
}





}
