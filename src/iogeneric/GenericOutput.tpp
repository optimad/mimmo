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
#if MIMMO_ENABLE_MPI
    if(m_rank == 0)
#endif
    {
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
}

/*!
 * Overloaded function of base class setInput.
 * It sets the input/result and write on file at the same time.
 * \param[in] data Data to be written and to be used to set the input/result.
 */
template<typename T>
void
GenericOutput::setInput(T data){
    setInput(&data);
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
    MimmoPiercedVector<T> * workingptr_ = data;

#if MIMMO_ENABLE_MPI
    std::unique_ptr< MimmoPiercedVector<T> > dataglobal(new MimmoPiercedVector<T>());
    collectDataFromAllProcs(*data, dataglobal.get());
    workingptr_ = dataglobal.get();
    // get only 0 to work;
    if(m_rank == 0)
#endif
    {
        std::fstream file;
        file.open(m_dir+"/"+m_filename, std::fstream::out);
        if (file.is_open()){
            if(m_binary){
                bitpit::genericIO::flushBINARY(file, loc);
                bitpit::genericIO::flushBINARY(file, long(workingptr_->size()));
                for (auto datait = workingptr_->begin(); datait != workingptr_->end(); ++datait) {
                    bitpit::genericIO::flushBINARY(file, datait.getId());
                    bitpit::genericIO::flushBINARY(file, *datait);
                }
            } else if (m_csv){
                outputCSVStream::ofstreamcsv(file, *workingptr_);
            }else{
                bitpit::genericIO::flushASCII(file,loc);
                file<<'\n';
                bitpit::genericIO::flushASCII(file,long(workingptr_->size()));
                file<<'\n';
                for (auto datait = workingptr_->begin(); datait != workingptr_->end(); ++datait) {
                    bitpit::genericIO::flushASCII(file, datait.getId());
                    bitpit::genericIO::flushASCII(file, *datait);
                    file<<'\n';
                }
            }
            file.close();
        }
    }// exiting scope writing.
}

#if MIMMO_ENABLE_MPI
/*!
 * For reading part only. Since we assume the 0 rank proc read from file,
 * we need to communicate data to all the other procs.
 * \param[in] locdata localdata to communicate
 * \param[in] dataTC collecting structure.
 */
template<typename T>
void
GenericOutputMPVData::collectDataFromAllProcs(MimmoPiercedVector<T> & locdata, MimmoPiercedVector<T> * globdata){

    if(m_rank == 0){
        //prefill global data with 0 rank info.
        globdata->clear();
        globdata->setDataLocation(locdata.getDataLocation());
        for(auto it=locdata.begin(); it!=locdata.end(); ++it){
            globdata->insert(it.getId(), *it);
        }


        //Receive data to all other procs and fill globdata
        for (int sendRank=1; sendRank<m_nprocs; sendRank++){
            long dataBufferSize;
            MPI_Recv(&dataBufferSize, 1, MPI_LONG, sendRank, 100, m_communicator, MPI_STATUS_IGNORE);
            bitpit::IBinaryStream dataBuffer(dataBufferSize);
            MPI_Recv(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, sendRank, 110, m_communicator, MPI_STATUS_IGNORE);

            //reverse in temp
            MimmoPiercedVector<T> temp;
            dataBuffer >> temp;

            //check values of temp and insert into globdata. If ID already exists skip it.
            for(auto it=temp.begin(); it!=temp.end(); ++it){
                if(!globdata->exists(it.getId())){
                    globdata->insert(it.getId(), *it);
                }
            }
        }
        //hey 0, your job is done.
    }else{

        bitpit::OBinaryStream dataBuffer;
        dataBuffer << locdata;
        long dataBufferSize = dataBuffer.getSize();
        //Send data to rank 0
        MPI_Send(&dataBufferSize, 1, MPI_LONG, 0, 100, m_communicator);
        MPI_Send(dataBuffer.data(), dataBuffer.getSize(), MPI_CHAR, 0, 110, m_communicator);
    }
}

#endif



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
