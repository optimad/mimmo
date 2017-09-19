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
std::fstream&
ifstreamcsv(std::fstream &in, T &x){
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
std::fstream&  ifstreamcsvend(std::fstream &in, T &x){
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
std::fstream&  ifstreamcsv(std::fstream &in, std::vector< T > &x){

    T       dummy;

    while (in.good()) {
        if (ifstreamcsv(in,dummy)) { x.push_back(dummy); }
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
std::fstream&  ifstreamcsvend(std::fstream &in, std::vector< T > &x){
        return ifstreamcsv(in, x);
}

/*!
 * Recover a data from a stream when import in csv format an array of data of dimension d.
 * \param[in] in import stream.
 * \param[out] x array data read.
 * \return Reference to result import stream after data import.
 */
template <typename T, size_t d>
std::fstream&  ifstreamcsv(std::fstream &in, std::array< T, d> &x){

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
 * Recover a data from a stream when import in csv format an array of data of dimension d.
 * \param[in] in import stream.
 * \param[out] x array data read.
 * \return Reference to result import stream after data import.
 */
template <typename T, size_t d>
std::fstream&  ifstreamcsvend(std::fstream &in, std::array< T, d> &x){
    return ifstreamcsvend(in, x);
};

/*!
 * Recover a data from a stream when import in csv format a MimmoPiercedVector of data.
 * \param[in] in import stream.
 * \param[out] x MimmoPiercedVector data read.
 * \return Reference to result import stream after data import.
 */
template <typename T>
std::fstream&  ifstreamcsv(std::fstream &in, MimmoPiercedVector< T > &x){

    T       dummy;
    long    id;
    int     location;
    long    sizeData = 0;

    if(ifstreamcsvend(in, location)){
        if(x.intIsValidLocation(location)) x.setDataLocation(static_cast<MPVLocation>(location));
    }
    if(ifstreamcsvend(in, sizeData)){
        x.reserve(sizeData);
    }
    
    for(long count = 0; count<sizeData; ++count){
        ifstreamcsv(in,id) ;
        ifstreamcsvend(in,dummy);
        x.insert(id,dummy); 
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
        std::fstream file;
        file.open(m_dir+"/"+m_filename, std::fstream::in);
        if (file.is_open()){
            if (m_csv){
                inputCSVStream::ifstreamcsv(file, data);
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
 * It gets the result of the object.
 * Result may be empty both for failed reading or invalid read data. In that case, 
 * return a NULL pointer.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
MimmoPiercedVector< T >*
GenericInputMPVData::_getResult(){
    MimmoPiercedVector< T > data;
    if (getGeometry() == NULL) return NULL;

    int n_loc;
    long nSize = 0;
    long id;
    T data_T;

    if(m_csv)   m_binary = false;
    data.setGeometry(getGeometry());

    std::fstream file;
    file.open(m_dir+"/"+m_filename, std::fstream::in);
    if (file.is_open()){
        if (m_binary){
            bitpit::genericIO::absorbBINARY(file, n_loc);
            bitpit::genericIO::absorbBINARY(file, nSize);
            data.reserve(nSize);
            for(int i=0; i<nSize; ++i){
                bitpit::genericIO::absorbBINARY(file, id);
                bitpit::genericIO::absorbBINARY(file, data_T);
                data.insert(id,data_T);
            }
            data.setDataLocation(n_loc);
        }
        else if (m_csv){
            inputCSVStream::ifstreamcsv(file, data);
        }else{
            bitpit::genericIO::absorbASCII(file, n_loc);
            bitpit::genericIO::absorbASCII(file, nSize);
            data.reserve(nSize);
            for(int i=0; i<nSize; ++i){
                bitpit::genericIO::absorbASCII(file, id);
                bitpit::genericIO::absorbASCII(file, data_T);
                data.insert(id,data_T);
            }
            data.setDataLocation(n_loc);
        }
        file.close();
    }else{
        (*m_log) << "file not open --> exit" << std::endl;
        throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
    }

    //mandatory check:if field is id-uncoherent with current geometry, does not set anything
    if(data.checkDataIdsCoherence()) _setResult(data);

    return static_cast<IODataT<MimmoPiercedVector< T > >*>(m_result.get())->getData();
}

/*!
 * It gets the result of the object.
 * Result may be empty both for failed reading or invalid read data. In that case, 
 * return an empty/default copy of the result data.
 * \return copy of the data stored in result member.
 */
template < typename T>
MimmoPiercedVector< T >
GenericInputMPVData::getResult(){
    MimmoPiercedVector< T > * pres = _getResult< T >() ;
    if (pres != NULL){
        return *pres;
    }else{
        return MimmoPiercedVector< T >();
    }
}

/*!
 * It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericInputMPVData::_setResult(MimmoPiercedVector< T > * data){
    m_result = std::move(std::unique_ptr<IOData> (new IODataT<MimmoPiercedVector< T > >(*data)));
}

/*!
 * It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericInputMPVData::_setResult(MimmoPiercedVector< T > & data){
    m_result = std::move(std::unique_ptr<IOData> (new IODataT<MimmoPiercedVector < T > >(data)));
}

}
