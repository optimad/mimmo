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
#ifndef __MIMMO_STREAMERS_HPP__
#define __MIMMO_STREAMERS_HPP__

#include "communications.hpp"
#include "MimmoPiercedVector.hpp"

namespace mimmo{

/*!
 * \class MimmoDataBufferStreamer
 * \ingroup parallel
 * \brief Specialized buffer streamer to exchange data defined on cells.
 */
template<class mpvt>
class MimmoDataBufferStreamer : public ExchangeBufferStreamer {

public:
	MimmoDataBufferStreamer(MimmoPiercedVector<mpvt > *data);

	void setData(MimmoPiercedVector<mpvt> *data);

    void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>());
    void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>());

private:

    MimmoPiercedVector<mpvt>* m_data;

};

/*!
 * \class MimmoPointDataBufferStreamerBase
 * \ingroup parallel
 * \brief Specialized buffer streamer to exchange data defined on points.
 */
template<class mpvt>
class MimmoPointDataBufferStreamer : public ExchangeBufferStreamer {

public:
    MimmoPointDataBufferStreamer(MimmoPiercedVector<mpvt> *data);

    void setData(MimmoPiercedVector<mpvt> *data);

    void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>());
    void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>());

private:

    MimmoPiercedVector<mpvt>* m_data;
    std::unordered_map<long, int> m_receivedFromRank; //Which rank sent the received id-th point data?

};
}

#include "mimmo_streamers.tpp"

#endif
#endif
