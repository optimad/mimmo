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

#ifndef __MIMMO_COMMUNICATIONS_HPP__
#define __MIMMO_COMMUNICATIONS_HPP__

#include <mpi.h>
#include <vector>
#include <unordered_map>

#include <bitpit_containers.hpp>
#include <bitpit_communications.hpp>
#include <bitpit_patchkernel.hpp>

#include <mimmo_core.hpp>

namespace mimmo{

class ExchangeBufferStreamer
{

public:
    ExchangeBufferStreamer(const size_t &itemSize);

    virtual ~ExchangeBufferStreamer();

    size_t getItemSize() const;

    virtual void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>()) = 0;
    virtual void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>()) = 0;

private:
    size_t m_itemSize;

    void setItemSize(const size_t &itemSize);

};

template<typename container_t, typename value_t = typename container_t::value_type>
class ListBufferStreamer : public ExchangeBufferStreamer
{

public:
    typedef value_t value_type;

    ListBufferStreamer(container_t *container);
    ListBufferStreamer(container_t *container, const size_t &itemSize);

    container_t & getContainer();

    void read(const int &rank, bitpit::RecvBuffer &buffer, const std::vector<long> &list = std::vector<long>());
    void write(const int &rank, bitpit::SendBuffer &buffer, const std::vector<long> &list = std::vector<long>());

protected:
    container_t *m_container;

};

class ListCommunicator : public bitpit::DataCommunicator
{

public:
    typedef std::vector<long> RankExchangeList;
    typedef std::unordered_map<int, RankExchangeList> ExchangeList;

    enum ListType {
        LIST_SEND,
        LIST_RECV
    };

    ListCommunicator(const MPI_Comm &communicator);

    virtual ~ListCommunicator();

    size_t getItemSize() const;

    const ExchangeList & getSendList() const;
    const RankExchangeList & getSendList(int rank) const;
    const ExchangeList & getRecvList() const;
    const RankExchangeList & getRecvList(int rank) const;
    virtual void setExchangeLists(const ExchangeList &sendList, const ExchangeList &recvList);
    virtual void setExchangeList(ListType listcontainer_type, const ExchangeList &list);

    bool hasData() const;
    void addData(ExchangeBufferStreamer *streamer);
    void addData(ExchangeBufferStreamer *writer, ExchangeBufferStreamer *reader);

    void startAllExchanges();
    void completeAllExchanges();

    void completeAllRecvs();
    int completeAnyRecv(const std::vector<int> &blacklist = std::vector<int>());

    void completeAllSends();

    void remapExchangeLists(const std::unordered_map<long, long> &mapper);
    void remapExchangeLists(const std::unordered_map<int, std::vector<long>> &mapper);

    void remapSendList(const std::unordered_map<long, long> &mapper);
    void remapSendList(const std::unordered_map<int, std::vector<long>> &mapper);

    void remapRecvList(const std::unordered_map<long, long> &mapper);
    void remapRecvList(const std::unordered_map<int, std::vector<long>> &mapper);

protected:
    size_t m_itemSize;
    ExchangeList m_sendList;
    ExchangeList m_recvList;

    void updateExchangeInfo();

    ExchangeList scatterExchangeList(const ExchangeList &inputList);

    void remapList(ExchangeList &list, const std::unordered_map<long, long> &mapper);
    void remapList(ExchangeList &list, const std::unordered_map<int, std::vector<long>> &mapper);

    virtual const RankExchangeList & getStreamableSendList(int rank, ExchangeBufferStreamer *reader);
    virtual const RankExchangeList & getStreamableRecvList(int rank, ExchangeBufferStreamer *writer);

private:
    std::vector<ExchangeBufferStreamer *> m_writers;
    std::vector<ExchangeBufferStreamer *> m_readers;

};

class GhostCommunicator : public ListCommunicator
{

public:
    GhostCommunicator(const bitpit::PatchKernel *patch);

    void resetExchangeLists();
    void setExchangeLists(const ExchangeList &sendList, const ExchangeList &recvList);
    void setExchangeList(ListType listcontainer_type, const ExchangeList &list);

protected:
    void createStreamableLists();

    const RankExchangeList & getStreamableSendList(int rank, ExchangeBufferStreamer *reader);
    const RankExchangeList & getStreamableRecvList(int rank, ExchangeBufferStreamer *writer);

private:
    const bitpit::PatchKernel *m_patch;

    ExchangeList m_sendListIds;
    ExchangeList m_recvListIds;

    ExchangeList sequentialIndexesConversion(const ExchangeList &list);

};


// POINT DATA COMMUNICATION STRUCTURES

class PointListCommunicator : public bitpit::DataCommunicator
{

public:
    typedef std::vector<long> RankExchangeList;
    typedef std::unordered_map<int, RankExchangeList> ExchangeList;

    enum ListType {
        LIST_SEND,
        LIST_RECV
    };

    PointListCommunicator(const MPI_Comm &communicator);

    virtual ~PointListCommunicator();

    size_t getItemSize() const;

    const ExchangeList & getSendList() const;
    const RankExchangeList & getSendList(int rank) const;
    const ExchangeList & getRecvList() const;
    const RankExchangeList & getRecvList(int rank) const;
    virtual void setExchangeLists(const ExchangeList &sendList, const ExchangeList &recvList);
    virtual void setExchangeList(ListType listcontainer_type, const ExchangeList &list);

    bool hasData() const;
    void addData(ExchangeBufferStreamer *streamer);
    void addData(ExchangeBufferStreamer *writer, ExchangeBufferStreamer *reader);

    void startAllExchanges();
    void completeAllExchanges();

    void completeAllRecvs();
    int completeAnyRecv(const std::vector<int> &blacklist = std::vector<int>());

    void completeAllSends();

    void remapExchangeLists(const std::unordered_map<long, long> &mapper);
    void remapExchangeLists(const std::unordered_map<int, std::vector<long>> &mapper);

    void remapSendList(const std::unordered_map<long, long> &mapper);
    void remapSendList(const std::unordered_map<int, std::vector<long>> &mapper);

    void remapRecvList(const std::unordered_map<long, long> &mapper);
    void remapRecvList(const std::unordered_map<int, std::vector<long>> &mapper);

protected:
    size_t m_itemSize;
    ExchangeList m_sendList;
    ExchangeList m_recvList;

    void updateExchangeInfo();

    ExchangeList scatterExchangeList(const ExchangeList &inputList);

    void remapList(ExchangeList &list, const std::unordered_map<long, long> &mapper);
    void remapList(ExchangeList &list, const std::unordered_map<int, std::vector<long>> &mapper);

    virtual const RankExchangeList & getStreamableSendList(int rank, ExchangeBufferStreamer *reader);
    virtual const RankExchangeList & getStreamableRecvList(int rank, ExchangeBufferStreamer *writer);

private:
    std::vector<ExchangeBufferStreamer *> m_writers;
    std::vector<ExchangeBufferStreamer *> m_readers;

};


class PointGhostCommunicator : public PointListCommunicator
{

public:
	PointGhostCommunicator(const MimmoObject *object);

    void resetExchangeLists();
    void setExchangeLists(const ExchangeList &sendList, const ExchangeList &recvList);
    void setExchangeList(ListType listcontainer_type, const ExchangeList &list);

protected:
    void createStreamableLists();

    const RankExchangeList & getStreamableSendList(int rank, ExchangeBufferStreamer *reader);
    const RankExchangeList & getStreamableRecvList(int rank, ExchangeBufferStreamer *writer);

private:
    const bitpit::PatchKernel *m_patch;
    const MimmoObject *m_object;

    ExchangeList m_sendListIds;
    ExchangeList m_recvListIds;

    ExchangeList sequentialIndexesConversion(const ExchangeList &list);

};

};

#include "communications.tpp"

#endif

#endif
