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
        if (ifstreamcsv(in,dummy)) {
            x[i] = dummy;
        }
        i++;
    } //next i
    if (ifstreamcsvend(in,dummy)) {
        x[i] = dummy;
    }
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
    return ifstreamcsv(in,x);
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
    long    sizeData = 0, readSizeData;

    if(ifstreamcsvend(in, location)){
        if(x.intIsValidLocation(location)) x.setDataLocation(static_cast<MPVLocation>(location));
    }
    if(ifstreamcsvend(in, readSizeData)){
        sizeData = std::max(sizeData, readSizeData);
    }
    x.reserve(sizeData);
    for(long count = 0; count<sizeData; ++count){
        ifstreamcsv(in,id) ;
        ifstreamcsvend(in,dummy);
        x.insert(id,dummy);
    }
    return(in);
}

/*!
 * Recover a data from a stream when import in csv format a MimmoPiercedVector of data.
 * \param[in] in import stream.
 * \param[out] x MimmoPiercedVector data read.
 * \return Reference to result import stream after data import.
 */
template<typename T, size_t d >
std::fstream&  ifstreamcsv(std::fstream &in, MimmoPiercedVector< std::array< T,d > > &x){

    std::array<T,d>  dummy;
    long    id;
    int     location;
    long    sizeData = 0, readSizeData;

    if(ifstreamcsvend(in, location)){
        if(x.intIsValidLocation(location)) x.setDataLocation(static_cast<MPVLocation>(location));
    }
    if(ifstreamcsvend(in, readSizeData)){
        sizeData = std::max(sizeData, readSizeData);
    }
    x.reserve(sizeData);
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
#if MIMMO_ENABLE_MPI
        if(m_rank == 0)
# endif
        {
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
                (*m_log) << "Error GenericInput: file not open --> exit" << std::endl;
                throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
            }

        }

#if MIMMO_ENABLE_MPI
        // 0 read something, let's share with all its friends.
        if(m_nprocs > 1)    sendReadDataToAllProcs(data);
#endif
        //finally store it in result slot.
        _setResult(data);
    }

    //ready to return the data.
    return  *static_cast<IODataT<T>*>(m_result.get())->getData();
}

#if MIMMO_ENABLE_MPI
/*!
 * For reading part only. Since we assume the 0 rank proc read from file,
 * we need to communicate data to all the other procs.
 * \param[in] dataTC data to be communicated.
 */
template<typename T>
void
GenericInput::sendReadDataToAllProcs(T & dataTC){

    if(m_rank == 0){

        //create char output data buffer and reverse data into it.
        bitpit::OBinaryStream dataBuffer;
        dataBuffer << dataTC;
        long dataBufferSize = dataBuffer.getSize();
        //Send data to all other procs
        for (int sendRank=1; sendRank<m_nprocs; sendRank++){
           MPI_Send(&dataBufferSize, 1, MPI_LONG, sendRank, 100, m_communicator);
           MPI_Send(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, sendRank, 110, m_communicator);
        }
        //hey 0, your job is done.
    }else{

        long dataBufferSize;
        MPI_Recv(&dataBufferSize, 1, MPI_LONG, 0, 100, m_communicator, MPI_STATUS_IGNORE);
        bitpit::IBinaryStream dataBuffer(dataBufferSize);
        MPI_Recv(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, 0, 110, m_communicator, MPI_STATUS_IGNORE);

        dataBuffer >> dataTC;
    }
}

#endif


//GENERICINPUTMPVDATA/////////////////////////////////////////////////////


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
    if (pres != nullptr){
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
    MimmoObject * refgeo = getGeometry();

#if MIMMO_ENABLE_MPI
    if(m_rank == 0)
#endif
    {
        int n_loc;
        long nSize = 0, readNSize;
        long id;
        T data_T;

        if(m_csv)   m_binary = false;

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
            	std::cout << " loc " << n_loc << std::endl;
                bitpit::genericIO::absorbASCII(file, readNSize);
            	std::cout << " readNSize " << readNSize << std::endl;
                nSize = std::max(nSize,readNSize);
                data.reserve(nSize);
                for(int i=0; i<nSize; ++i){
                    bitpit::genericIO::absorbASCII(file, id);
                    bitpit::genericIO::absorbASCII(file, data_T);
                	std::cout << " id " << id << std::endl;
                	std::cout << " data_T " << data_T << std::endl;
                    data.insert(id,data_T);
                }
                data.setDataLocation(n_loc);
            }
            file.close();
        }else{
            (*m_log) << "file not open --> exit" << std::endl;
            throw std::runtime_error (m_name + " : cannot open " + m_filename + " requested");
        }
    }

#if MIMMO_ENABLE_MPI
    //if there are any other procs send data to them.
    if(m_nprocs > 1)    sendReadDataToAllProcs(data);
#endif

    if(refgeo == nullptr){
        _setResult(data);
    }else{
        data.setGeometry(refgeo);
        MimmoPiercedVector<T> data2 = data.resizeToCoherentDataIds();
        _setResult(data2);
    }

    return static_cast<IODataT<MimmoPiercedVector< T > >*>(m_result.get())->getData();
}


#if MIMMO_ENABLE_MPI
/*!
 * For reading part only. Since we assume the 0 rank proc read from file,
 * we need to communicate data to all the other procs.
 * \param[in] dataTC data to be communicated.
 */
template<typename T>
void
GenericInputMPVData::sendReadDataToAllProcs(MimmoPiercedVector<T> & dataTC){

    if(m_rank == 0){

        bitpit::OBinaryStream dataBuffer;
        //create char output data buffer and reverse data into it.
        dataBuffer << dataTC;
        long dataBufferSize = dataBuffer.getSize();

        //Send data to all other procs
        for (int sendRank=1; sendRank<m_nprocs; sendRank++){
           MPI_Send(&dataBufferSize, 1, MPI_LONG, sendRank, 100, m_communicator);
           MPI_Send(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, sendRank, 110, m_communicator);
        }
        //hey 0, your job is done.
    }else{

        long dataBufferSize;
        MPI_Recv(&dataBufferSize, 1, MPI_LONG, 0, 100, m_communicator, MPI_STATUS_IGNORE);

        bitpit::IBinaryStream dataBuffer(dataBufferSize);
        MPI_Recv(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, 0, 110, m_communicator, MPI_STATUS_IGNORE);

        dataBuffer >> dataTC;
    }
}

#endif

}
