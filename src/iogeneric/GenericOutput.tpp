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

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Pointer to data to be written and to be used to set the input/result.
 */
template<typename T>
void
GenericOutput::setInput(T* data){
    _setInput(data);
    _setResult(data);
    std::ofstream file;
    file.open(m_dir+"/"+m_filename);
    if (file.is_open()){
        if (m_csv){
            ofstreamcsv(file, *data);
        }
        else{
            file << *data;
        }
        file.close();
    }
}

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Data to be written and to be used to set the input/result.
 */
template<typename T>
void
GenericOutput::setInput(T data){
    _setInput(data);
    _setResult(data);
    std::ofstream file;
    file.open(m_dir+"/"+m_filename);
    if (file.is_open()){
        if (m_csv){
            ofstreamcsv(file, data);
        }
        else{
            file << data;
        }
        file.close();
    }
}

/*!
 * It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericOutput::getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

/*!
 * It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericOutput::setResult(T* data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_result = std::move(dummy);
}

/*!
 * It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericOutput::setResult(T data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_result = std::move(dummy);
}

/*!
 * It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
GenericOutput::getResult(){
    return(static_cast<IODataT<T>*>(m_result.get())->getData());
}

/*!
 * It sets the input member of the object.
 * \param[in] data Pointer to data to be stored in the input member.
 */
template<typename T>
void
GenericOutput::_setInput(T* data){
    clearInput();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_input = std::move(dummy);
}

/*!It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
GenericOutput::_setInput(T data){
    clearInput();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_input = std::move(dummy);
}

/*!It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
GenericOutput::_getInput(){
    return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

/*!It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
GenericOutput::_setResult(T* data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
    m_result = std::move(dummy);
}

/*!It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
GenericOutput::_setResult(T data){
    clearResult();
    std::unique_ptr<IOData> dummy(new IODataT<T>(data));
    m_result = std::move(dummy);
}

/*!It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
GenericOutput::_getResult(){
    return(static_cast<IODataT<T>*>(m_result.get())->getData());
}

/*!
 * Store a data in a stream when export in csv format.
 * \param[in] out output stream.
 * \param[in] x data to write.
 * \return Reference to result output stream after data export.
 */
template <class T>
std::ofstream& GenericOutput::ofstreamcsv(std::ofstream &out, const T &x)
{
    out << x << ",";
    return(out);
};

/*!
 * Store a data in a stream when export in csv format a vector of data.
 * \param[in] out output stream.
 * \param[in] x vector of data to write.
 * \return Reference to result output stream after data export.
 */
template <class T>
std::ofstream& GenericOutput::ofstreamcsv(std::ofstream &out, const std::vector< T > &x)
{

    size_t n = x.size();
    if (n == 0) {
        return(out);
    }
    for (size_t i = 0; i < n-1; i++) {
        ofstreamcsv(out,x[i]);
    } //next i
    ofstreamcsvend(out,x[n-1]);
    out << "\n";
    return(out);
};

/*!
 * Store a data in a stream when export in csv format an array of data of dimension d.
 * \param[in] out output stream.
 * \param[in] x array of data to write.
 * \return Reference to result output stream after data export.
 */
template <class T, size_t d>
std::ofstream& GenericOutput::ofstreamcsv(std::ofstream &out, const std::array< T,d > &x)
{

    if (d == 0) return(out);
    for (size_t i = 0; i < d-1; i++) {
        ofstreamcsv(out,x[i]);
    } //next i
    ofstreamcsvend(out,x[d-1]);
    out << "\n";
    return(out);
};

/*!
 * Store a data in a stream when export in csv format when end of line is reached.
 * \param[in] out output stream.
 * \param[in] x data to write.
 * \return Reference to result output stream after data export.
 */
template <class T>
std::ofstream& GenericOutput::ofstreamcsvend(std::ofstream &out, const T &x)
{
    out << x;
    return(out);
};

/*!
 * Store a data in a stream when export in csv format a vector of data when end of line is reached.
 * \param[in] out output stream.
 * \param[in] x vector of data to write.
 * \return Reference to result output stream after data export.
 */
template <class T>
std::ofstream& GenericOutput::ofstreamcsvend(std::ofstream &out, const std::vector< T > &x)
{

    ofstreamcsv(out,x);
    return(out);
};

/*!
 * Store a data in a stream when export in csv format an array of data
 * of dimension d when end of line is reached.
 * \param[in] out output stream.
 * \param[in] x array of data to write.
 * \return Reference to result output stream after data export.
 */
template <class T, size_t d>
std::ofstream& GenericOutput::ofstreamcsvend(std::ofstream &out, const std::array< T,d > &x)
{

    ofstreamcsv(out,x);
    return(out);
};

}
