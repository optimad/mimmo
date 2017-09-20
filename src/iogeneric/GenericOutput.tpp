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

namespace outputCSVStream {

    /*!
     * Store a data in a stream when export in csv format.
     * \param[in] out output stream.
     * \param[in] x data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsv(std::fstream &out, const T &x)
    {
        out << x << ",";
        return(out);
    };
    
    /*!
     * Store a data in a stream when export in csv format when end of line is reached.
     * \param[in] out output stream.
     * \param[in] x data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsvend(std::fstream &out, const T &x)
    {
        out << x;
        out<<'\n';
        return(out);
    };
    
    /*!
     * Store a data in a stream when export in csv format a vector of data.
     * \param[in] out output stream.
     * \param[in] x vector of data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsv(std::fstream &out, const std::vector< T > &x)
    {

        size_t n = x.size();
        if (n == 0) {
            return(out);
        }
        for (size_t i = 0; i < n-1; i++) {
            ofstreamcsv(out,x[i]);
        } //next i
        ofstreamcsvend(out,x[n-1]);
        return(out);
    };

    /*!
     * Store a data in a stream when export in csv format a vector of data when end of line is reached.
     * \param[in] out output stream.
     * \param[in] x vector of data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsvend(std::fstream &out, const std::vector< T > &x)
    {
        ofstreamcsv(out,x);
        return(out);
    };

    /*!
     * Store a data in a stream when export in csv format an array of data of dimension d.
     * \param[in] out output stream.
     * \param[in] x array of data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T, size_t d>
    std::fstream& ofstreamcsv(std::fstream &out, const std::array< T,d > &x)
    {
        
        if (d == 0) return(out);
        for (size_t i = 0; i < d-1; i++) {
            ofstreamcsv(out,x[i]);
        } //next i
        ofstreamcsvend(out,x[d-1]);
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
    std::fstream& ofstreamcsvend(std::fstream &out, const std::array< T,d > &x)
    {
        ofstreamcsv(out,x);
        return(out);
    };

    /*!
     * Store a data in a stream when export in csv format a vector of data.
     * \param[in] out output stream.
     * \param[in] x vector of data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsv(std::fstream &out, const MimmoPiercedVector< T > &x)
    {
        size_t n = x.size();
        if (n == 0) {
            return(out);
        }
        int loc = static_cast<int>(x.getConstDataLocation());
        out << loc<< '\n';
        out << long(x.size())<<'\n';
        typename bitpit::PiercedVector<T>::const_iterator itB, itE = x.end();
        for (itB = x.begin(); itB != itE; ++itB) {
            ofstreamcsv(out, itB.getId());
            ofstreamcsvend(out,*itB);
        } //next i
        return(out);
    };

    /*!
     * Store a data in a stream when export in csv format a vector of data when end of line is reached.
     * \param[in] out output stream.
     * \param[in] x vector of data to write.
     * \return Reference to result output stream after data export.
     */
    template <class T>
    std::fstream& ofstreamcsvend(std::fstream &out, const MimmoPiercedVector< T > &x)
    {
        ofstreamcsv(out,x);
        return(out);
    };

}//end of namespace

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Pointer to data to be written and to be used to set the input/result.
 */
template<typename T>
void
GenericOutput::setInput(T* data){
    _setInput(*data);
    std::fstream file;
    file.open(m_dir+"/"+m_filename, std::fstream::out);
    if (file.is_open()){
        if (m_csv){
            outputCSVStream::ofstreamcsv(file, *data);
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
    std::fstream file;
    file.open(m_dir+"/"+m_filename, std::fstream::out);
    if (file.is_open()){
        if (m_csv){
            outputCSVStream::ofstreamcsv(file, data);
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

/*!It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
GenericOutput::_setInput(T & data){
    m_input = std::move(std::unique_ptr<IOData>(new IODataT<T>(data)));
}

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Pointer to data to be written and to be used to set the input/result.
 */
template<typename T>
void
GenericOutputMPVData::setInput(MimmoPiercedVector< T > * data){
    _setInput(*data);
    if(m_csv)   m_binary = false;
    int loc = static_cast<int>(data->getDataLocation());
    std::fstream file;
    file.open(m_dir+"/"+m_filename, std::fstream::out);
    if (file.is_open()){
        if(m_binary){
            bitpit::genericIO::flushBINARY(file, loc);
            bitpit::genericIO::flushBINARY(file, long(data->size()));
            typename bitpit::PiercedVector<T>::const_iterator dataItr, dataEnd = data->end();
            for (dataItr = data->begin(); dataItr != dataEnd; ++dataItr) {
                bitpit::genericIO::flushBINARY(file, dataItr.getId());
                bitpit::genericIO::flushBINARY(file, *dataItr);
            }    
        } else if (m_csv){
            outputCSVStream::ofstreamcsv(file, *data);
        }else{
            bitpit::genericIO::flushASCII(file,loc);
            file<<'\n';
            bitpit::genericIO::flushASCII(file,long(data->size()));
            file<<'\n';
            typename bitpit::PiercedVector<T>::const_iterator dataItr, dataEnd = data->end();
            for (dataItr = data->begin(); dataItr != dataEnd; ++dataItr) {
                bitpit::genericIO::flushASCII(file, dataItr.getId());
                bitpit::genericIO::flushASCII(file, *dataItr);
                file<<'\n';
            }
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
GenericOutputMPVData::setInput(MimmoPiercedVector< T > data){
    setInput(&data);
}

/*!
 * It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
MimmoPiercedVector< T >*
GenericOutputMPVData::getInput(){
    return(static_cast<IODataT<MimmoPiercedVector< T > >*>(m_input.get())->getData());
}


/*!It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
GenericOutputMPVData::_setInput(MimmoPiercedVector< T > & data){
    m_input = std::move(std::unique_ptr<IOData>(new IODataT<MimmoPiercedVector< T > >(data)));
}


}
