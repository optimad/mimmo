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
template<std::size_t NCOMP>
class MimmoDataBufferStreamer : public ExchangeBufferStreamer {

public:
	MimmoDataBufferStreamer(MimmoPiercedVector<std::array<double, NCOMP> > *data);

	void setData(MimmoPiercedVector<std::array<double, NCOMP> > *data);

    void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>());
    void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>());

private:

    MimmoPiercedVector<std::array<double, NCOMP> >* m_data;

};

/*!
 * \class MimmoPointDataBufferStreamer
 * \ingroup parallel
 * \brief Specialized buffer streamer to exchange data defined on points.
 */
template<std::size_t NCOMP>
class MimmoPointDataBufferStreamer : public ExchangeBufferStreamer {

public:
	MimmoPointDataBufferStreamer(MimmoPiercedVector<std::array<double, NCOMP> > *data);

	void setData(MimmoPiercedVector<std::array<double, NCOMP> > *data);

    void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>());
    void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>());

private:

    MimmoPiercedVector<std::array<double, NCOMP> >* m_data;
    std::unordered_map<long, int> m_receivedFromRank; //Which rank sent the received id-th point data?

};

}

#include "mimmo_streamers.tpp"

#endif
#endif
