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
#if MIMMO_ENABLE_MPI

namespace mimmo{

/*!
    Creates a new scalar streamer

    \param data is the dataset that will be exchanged
*/
template<class mpvt>
MimmoDataBufferStreamer<mpvt>::MimmoDataBufferStreamer(MimmoPiercedVector<mpvt>* data)
    : ExchangeBufferStreamer(sizeof(mpvt))
{
	m_data = data;
}

template<class mpvt>
void MimmoDataBufferStreamer<mpvt>::setData(MimmoPiercedVector<mpvt> *data)
{
	m_data = data;
}

/*!
    Read the dataset from the buffer.

    \param rank is the rank of the processor who sent the data
    \param buffer is the buffer where the data will be read from
    \param list is the list of ids that will be read
*/
template<class mpvt>
void MimmoDataBufferStreamer<mpvt>::read(int const &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list)
{
    BITPIT_UNUSED(rank);

    // Read the dataset
    for (const long id : list) {
    	mpvt val;
    	buffer >> val;
    	m_data->at(id) = val;
    }
}

/*!
    Write the dataset into the buffer.

    \param rank is the rank of the processor who will receive the data
    \param buffer is the buffer where the data will be written to
    \param list is the list of ids that will be written
*/
template<class mpvt>
void MimmoDataBufferStreamer<mpvt>::write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list)
{
    BITPIT_UNUSED(rank);

    // Write the dataset
    for (const long id : list) {
    	mpvt val = m_data->at(id);
    	buffer << val;
    }
}

///*!
//    Creates a new scalar streamer
//
//    \param data is the dataset that will be exchanged
//*/
//template<std::size_t NCOMP>
//MimmoPointDataBufferStreamer<NCOMP>::MimmoPointDataBufferStreamer(MimmoPiercedVector<std::array<double, NCOMP> >* data)
//    : ExchangeBufferStreamer(NCOMP*sizeof(double))
//{
//	m_data = data;
//}
//
//template<std::size_t NCOMP>
//void MimmoPointDataBufferStreamer<NCOMP>::setData(MimmoPiercedVector<std::array<double, NCOMP> > *data)
//{
//	m_data = data;
//}
//
///*!
//    Read the dataset from the buffer.
//
//    \param rank is the rank of the processor who sent the data
//    \param buffer is the buffer where the data will be read from
//    \param list is the list of ids that will be read
//*/
//template<std::size_t NCOMP>
//void MimmoPointDataBufferStreamer<NCOMP>::read(int const &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list)
//{
//    // Read the dataset
//    for (const long id : list) {
//    	std::array<double, NCOMP> value;
//    	for (double & val : value)
//    		buffer >> val;
//    	if (!m_receivedFromRank.count(id)){
//    		m_data->at(id) = value;
//    		m_receivedFromRank[id] = rank;
//    	}
//    	else if(rank <= m_receivedFromRank[id]){
//    		m_data->at(id) = value;
//    		m_receivedFromRank[id] = rank;
//    	}
//    }
//}
//
///*!
//    Write the dataset into the buffer.
//
//    \param rank is the rank of the processor who will receive the data
//    \param buffer is the buffer where the data will be written to
//    \param list is the list of ids that will be written
//*/
//template<std::size_t NCOMP>
//void MimmoPointDataBufferStreamer<NCOMP>::write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list)
//{
//    BITPIT_UNUSED(rank);
//
//    // Write the dataset
//    for (const long id : list) {
//    	std::array<double, NCOMP> value = m_data->at(id);
//    	for (double & val : value)
//    		buffer << val;
//    }
//}


/*!
    Creates a new base streamer

    \param data is the dataset that will be exchanged
*/
template<class mpvt>
MimmoPointDataBufferStreamer<mpvt>::MimmoPointDataBufferStreamer(MimmoPiercedVector<mpvt>* data)
    : ExchangeBufferStreamer(sizeof(mpvt))
{
    m_data = data;
}

template<class mpvt>
void MimmoPointDataBufferStreamer<mpvt>::setData(MimmoPiercedVector<mpvt> *data)
{
    m_data = data;
}

/*!
    Read the dataset from the buffer.

    \param rank is the rank of the processor who sent the data
    \param buffer is the buffer where the data will be read from
    \param list is the list of ids that will be read
*/
template<class mpvt>
void MimmoPointDataBufferStreamer<mpvt>::read(int const &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list)
{
    // Read the dataset
    for (const long id : list) {
        mpvt val;
        buffer >> val;
        if (!m_receivedFromRank.count(id)){
            m_data->at(id) = val;
            m_receivedFromRank[id] = rank;
        }
        else if(rank <= m_receivedFromRank[id]){
            m_data->at(id) = val;
            m_receivedFromRank[id] = rank;
        }
    }
}

/*!
    Write the dataset into the buffer.

    \param rank is the rank of the processor who will receive the data
    \param buffer is the buffer where the data will be written to
    \param list is the list of ids that will be written
*/
template<class mpvt>
void MimmoPointDataBufferStreamer<mpvt>::write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list)
{
    BITPIT_UNUSED(rank);

    // Write the dataset
    for (const long id : list) {
        mpvt val = m_data->at(id);
        buffer << val;
    }
}

}
#endif
