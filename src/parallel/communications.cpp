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

#include "communications.hpp"

namespace mimmo{

/*!
    \class ExchangeBufferStreamer

    \brief The ExchangeBufferStreamer class provides the base infrastructure
    needed to write and read the buffer of a DataCommunicator.
*/

ExchangeBufferStreamer::ExchangeBufferStreamer(const size_t &itemSize)
{
    setItemSize(itemSize);
}

/*!
    Destructuor.
*/
ExchangeBufferStreamer::~ExchangeBufferStreamer()
{
}

/*!
    Gets the size, expressed in bytes, of the single list item to exchange.

    \result The size, expressed in bytes, of the single list item to exchange.
*/
size_t ExchangeBufferStreamer::getItemSize() const
{
    return m_itemSize;
}

/*!
    Sets the size, expressed in bytes, of the single list item to exchange.

    \param itemSize is the size, expressed in bytes, of the single list
    item to exchange
*/
void ExchangeBufferStreamer::setItemSize(const size_t &itemSize)
{
    m_itemSize = itemSize;
}

/*!
    \class ListCommunicator

    \brief The ListCommunicator class provides the infrastructure needed to
    exchange a list of data among processors.
*/

/*!
    Creates a new communicator for data exchange

    \param patch is the patch
*/
ListCommunicator::ListCommunicator(const MPI_Comm &communicator)
    : DataCommunicator(communicator),
      m_itemSize(0)
{
}

/*!
    Destructor.
*/
ListCommunicator::~ListCommunicator()
{
}

/*!
    Gets the size, expressed in bytes, of the single list item to exchange.

    \result The size, expressed in bytes, of the single list item to exchange.
*/
size_t ListCommunicator::getItemSize() const
{
    return m_itemSize;
}

/*!
    Gets the list of ids for which the current processor will send data.

    \result The list of ids for which the current processor will send data.
*/
const ListCommunicator::ExchangeList & ListCommunicator::getSendList() const
{
    return m_sendList;
}

/*!
    Gets the list of ids for which the current processor will send data
    to the specified rank.

    \result The list of ids for which the current processor will send data
    to the specified rank.
*/
const ListCommunicator::RankExchangeList & ListCommunicator::getSendList(int rank) const
{
    return m_sendList.at(rank);
}

/*!
    Gets the list of ids from which the current processor expects to receive
    data.

    \result The list of ids from which the current processor expects to receive
    data.
*/
const ListCommunicator::ExchangeList & ListCommunicator::getRecvList() const
{
    return m_recvList;
}

/*!
    Gets the list of ids from which the current processor expects to receive
    data from the specified rank.

    \result The list of ids from which the current processor expects to receive
    data from the specified rank.
*/
const ListCommunicator::RankExchangeList & ListCommunicator::getRecvList(int rank) const
{
    return m_recvList.at(rank);
}

/*!
    Sets both the send list and the receive list.

    \param sendList is a list of ids the current processor will send data to
    \param recvList is a list of ids the current processor expects to receive
    data from
*/
void ListCommunicator::setExchangeLists(const ExchangeList &sendList,
                                        const ExchangeList &recvList)
{
    m_sendList = sendList;
    m_recvList = recvList;

    if (hasData()) {
        updateExchangeInfo();
    }
}

/*!
    Sets either the send list or the receive list.

    \param listType controls if the list defines ids to send or to receive
    \param list is a list of ids the current processor expectst to receive
    data from or will send data for
*/
void ListCommunicator::setExchangeList(ListType listType, const ExchangeList &list)
{
    switch (listType) {
        case LIST_SEND:
            m_sendList = list;
            m_recvList = scatterExchangeList(m_sendList);
            break;

        case LIST_RECV:
            m_recvList = list;
            m_sendList = scatterExchangeList(m_recvList);
            break;
    }

    if (hasData()) {
        updateExchangeInfo();
    }
}

/*!
    Communicate an exchange list among the other processors

    \param inputList is the list the will be communicated
    \result The exchange list recevied from the other processors.

*/
ListCommunicator::ExchangeList ListCommunicator::scatterExchangeList(const ExchangeList &inputList)
{
    // Clear all requests
    //
    // Before setting new requests it is important to cancel all the old ones.
    // Otherwise, if the continuous receive feature is enabled, the new sends
    // may start sending data to the old receives. The synchronous version
    // of clear functions will be used.
    //
    // When using the synchronous version of the clear functions, send and
    // receive requests have to match. It's not possible to first clear the
    // sends synchronously and then clear the receives synchronously, when
    // clearing the receives the sends would have already been cleared and
    // thus the send/receive requests would not match. We need to first
    // cancel the requests and then clear them.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    //
    // Send the list of ids
    //

    // Set the sends for exchanging the ids
    for (const auto &entry : inputList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * sizeof(long);

        setSend(rank, bufferSize);
    }

    // Discover the receives for exchanging the ids
    discoverRecvs();

    // Start the receives
    if (!areRecvsContinuous()) {
        startAllRecvs();
    }

    // Send the list of ids
    for (int rank : getSendRanks()) {
    	bitpit::SendBuffer buffer = getSendBuffer(rank);
        for (const long id : inputList.at(rank)) {
            buffer << id;
        }

        startSend(rank);
    }

    //
    // Get the list of ids
    //

    // Receive the ids
    int nPendingRecvs = getRecvCount();

    ExchangeList outputList;
    outputList.reserve(nPendingRecvs);
    while (nPendingRecvs != 0) {
        int rank = waitAnyRecv();
        bitpit::RecvBuffer buffer = getRecvBuffer(rank);
        long rankListSize = buffer.getSize() / sizeof(long);

        RankExchangeList &rankList = outputList[rank];
        rankList.resize(rankListSize);
        for (size_t i = 0; i < rankList.size(); ++i) {
            buffer >> rankList[i];
        }

        --nPendingRecvs;
    }

    // Wait all sends
    waitAllSends();

    // Clear the requests used for scattering the list
    //
    // Also here the synchronous version of clear functions will be used.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    return outputList;
}

/*!
    Updates the information needed for data exchange
*/
void ListCommunicator::updateExchangeInfo()
{
    // Clear all requests
    //
    // Before setting new requests it is important to cancel all the old ones.
    // Otherwise, if the continuous receive feature is enabled, the new sends
    // may start sending data to the old receives. The synchronous version
    // of clear functions will be used.
    //
    // When using the synchronous version of the clear functions, send and
    // receive requests have to match. It's not possible to first clear the
    // sends synchronously and then clear the receives synchronously, when
    // clearing the receives the sends would have already been cleared and
    // thus the send/receive requests would not match. We need to first
    // cancel the requests and then clear them.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    // Update send info
    for (const auto &entry : m_sendList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * m_itemSize;

        setSend(rank, bufferSize);
    }

    // Update receive info
    for (const auto &entry : m_recvList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * m_itemSize;

        setRecv(rank, bufferSize);
    }
}

/*!
    Checks if the communicator has some data to communicate.

    \result Returns true if the communicator has some data to communicate,
    otherwise it returns false.
*/
bool ListCommunicator::hasData() const
{
    return (m_writers.size() > 0 || m_readers.size() > 0);
}

/*!
    Add the data of the streamer to the communicator.

    \param streamer is the streamer of the data
*/
void ListCommunicator::addData(ExchangeBufferStreamer *streamer)
{
    addData(streamer, streamer);
}

/*!
    Add the data of the streamers to the communicator.

    \param writer is the streamer that will write the data to the buffer
    \param reader is the streamer that will read the data from the buffer
*/
void ListCommunicator::addData(ExchangeBufferStreamer *writer, ExchangeBufferStreamer *reader)
{
    // Update the item size
    int writerItemSize = writer->getItemSize();
    int readerItemSize = reader->getItemSize();
    if (readerItemSize != writerItemSize) {
        throw std::runtime_error("The item size of the writer differs from the item size of the reader");
    }

    m_itemSize += writerItemSize;

    // Add the streamers
    m_writers.push_back(writer);
    m_readers.push_back(reader);

    // Update exchange info
    updateExchangeInfo();
}

/*!
    Send ghosts data using non-blocking communications

    \param cellData is the container of the cell data
*/
void ListCommunicator::startAllExchanges()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    // Start the receives
    for (int rank : getRecvRanks()) {
        if (!isRecvActive(rank)) {
            startRecv(rank);
        }
    }

    // Wait previous sends
    waitAllSends();

    // Fill the buffer with the given field and start sending the data
    for (int rank : getSendRanks()) {
        // Get send buffer
        bitpit::SendBuffer &buffer = getSendBuffer(rank);

        // Write the buffer
        for (ExchangeBufferStreamer *streamer : m_writers) {
            streamer->write(rank, buffer, getStreamableSendList(rank, streamer));
        }

        // Start the send
        startSend(rank);
    }
}

/*!
    Receive ghosts data using non-blocking communications
*/
void ListCommunicator::completeAllExchanges()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    completeAllRecvs();

    completeAllSends();
}

/*!
    Receive ghosts data using non-blocking communications
*/
void ListCommunicator::completeAllRecvs()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    size_t nTotalRecvs = getRecvCount();

    std::vector<int> completedRecvs;
    while (completedRecvs.size() != nTotalRecvs) {
        int rank = completeAnyRecv(completedRecvs);
        completedRecvs.push_back(rank);
    }
}

/*!
    Receive ghost data from any receive and return the associated rank.

    \result The rank of the completed receive or MPI_UNDEFINED if there was
    no active receives.
*/
int ListCommunicator::completeAnyRecv(const std::vector<int> &blacklist)
{
    // Wait for a receve to finish

    int rank = waitAnyRecv(blacklist);

    // Get receive buffer
    bitpit::RecvBuffer &buffer = getRecvBuffer(rank);

    // Read the buffer
    for (ExchangeBufferStreamer *streamer : m_readers) {
        streamer->read(rank, buffer, getStreamableRecvList(rank, streamer));
    }

    return rank;
}

/*!
    Receive ghosts data using non-blocking communications
*/
void ListCommunicator::completeAllSends()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    waitAllSends();
}

/*!
    Remap the exchange lists according to the specified mapper.

    \param mapper is the mapper used for mapping the exchange lists
*/
void ListCommunicator::remapExchangeLists(const std::unordered_map<long, long> &mapper)
{
    remapRecvList(mapper);
    remapSendList(mapper);
}

/*!
    Remap the exchange lists according to the specified mapper.

    \param mapper is the mapper used for mapping the exchange lists
*/
void ListCommunicator::remapExchangeLists(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapRecvList(mapper);
    remapSendList(mapper);
}

/*!
    Remap the send list according to the specified mapper.

    \param mapper is the mapper used for mapping the send list
*/
void ListCommunicator::remapSendList(const std::unordered_map<long, long> &mapper)
{
    remapList(m_sendList, mapper);
}

/*!
    Remap the send list according to the specified mapper.

    \param mapper is the mapper used for mapping the send list
*/
void ListCommunicator::remapSendList(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapList(m_sendList, mapper);
}

/*!
    Remap the receive list according to the specified mapper.

    \param mapper is the mapper used for mapping the receive list
*/
void ListCommunicator::remapRecvList(const std::unordered_map<long, long> &mapper)
{
    remapList(m_recvList, mapper);
}

/*!
    Remap the receive list according to the specified mapper.

    \param mapper is the mapper used for mapping the receive list
*/
void ListCommunicator::remapRecvList(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapList(m_recvList, mapper);
}

/*!
    Remap the list according to the specified mapper.

    \param list it the list to remap
    \param mapper is the mapper used for mapping the list
*/
void ListCommunicator::remapList(ExchangeList &list, const std::unordered_map<long, long> &mapper)
{
    for (auto &entry : list) {
        for (long &id : entry.second) {
            id = mapper.at(id);
        }
    }
}

/*!
    Remap the list according to the specified mapper.

    \param list it the list to remap
    \param mapper is the mapper used for mapping the list
*/
void ListCommunicator::remapList(ExchangeList &list, const std::unordered_map<int, std::vector<long>> &mapper)
{
    for (auto &entry : list) {
        const std::vector<long> &rankMapper = mapper.at(entry.first);
        for (long &id : entry.second) {
            id = rankMapper[id];
        }
    }
}

/*!
    Gets a constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be send to
    \param reader is the streamer that will read the buffer
    \result A constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.
*/
const ListCommunicator::RankExchangeList & ListCommunicator::getStreamableSendList(int rank, ExchangeBufferStreamer *reader)
{
    BITPIT_UNUSED(reader);

    return m_sendList.at(rank);
}

/*!
    Gets a constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be received from
    \param writer is the streamer that will write the buffer
    \result A constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.
*/
const ListCommunicator::RankExchangeList & ListCommunicator::getStreamableRecvList(int rank, ExchangeBufferStreamer *writer)
{
    BITPIT_UNUSED(writer);

    return m_recvList.at(rank);
}

/*!
    \class GhostCommunicator

    \brief The GhostCommunicator class provides the infrastructure
    needed to exchange data among ghosts.
*/

/*!
    Creates a new communicator for data exchange

    \param patch is the patch

*/
GhostCommunicator::GhostCommunicator(const bitpit::PatchKernel *patch)
    : ListCommunicator(patch->getCommunicator()),
      m_patch(patch)
{
}

/*!
    Resets the exchange list.

    After the exchange list has been reset, data will be exchanged among all
    ghosts.
*/
void GhostCommunicator::resetExchangeLists()
{
    const ExchangeList &exchangeSources = m_patch->getGhostExchangeSources();
    const ExchangeList &exchangeTargets = m_patch->getGhostExchangeTargets();

    ExchangeList sendList = sequentialIndexesConversion(exchangeSources);
    ExchangeList recvList = sequentialIndexesConversion(exchangeTargets);

    setExchangeLists(sendList, recvList);
}

/*!
    Sets both the send list and the receive list.

    \param sendList is a list of ghosts indexes the current processor will
    send data to
    \param recvList is a list of ghosts indexes the current processor expects
    to receive data from
*/
void GhostCommunicator::setExchangeLists(const ExchangeList &sendList,
                                        const ExchangeList &recvList)
{
    // Set the index lists
    ListCommunicator::setExchangeLists(sendList, recvList);

    // Create lists for the streamers
    createStreamableLists();
}

/*!
    Sets either the send list or the receive list.

    \param listType controls if the list defines indexes to send or to receive
    \param list is a list of ghosts indexes the current processor expectst to
    receive data from or will send data for
*/
void GhostCommunicator::setExchangeList(ListType listType, const ExchangeList &list)
{
    // Set the index lists
    ListCommunicator::setExchangeList(listType, list);

    // Create lists for the streamers
    createStreamableLists();
}

/*!
    Converts a list of ghost ids to a list of sequential ghost indexes.

    \param idList is the list of ghost ids
    \return The list converted to sequential ghost indexes.
*/
GhostCommunicator::ExchangeList GhostCommunicator::sequentialIndexesConversion(const ExchangeList &idList)
{
    ExchangeList indexList;
    indexList.reserve(idList.size());

    for (const auto &entry : idList) {
        short rank = entry.first;
        const RankExchangeList &rankIdList = entry.second;

        RankExchangeList &rankIndexList = indexList[rank];
        rankIndexList.reserve(rankIdList.size());
        for (size_t i = 0; i < rankIdList.size(); ++i) {
            rankIndexList.push_back(i);
        }
    }

    return indexList;
}

/*!
    Gets a constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be send to
    \param reader is the streamer that will read the buffer
    \result A constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.
*/
const GhostCommunicator::RankExchangeList & GhostCommunicator::getStreamableSendList(int rank, ExchangeBufferStreamer *reader)
{
    BITPIT_UNUSED(reader);

    return m_sendListIds.at(rank);
}

/*!
    Creates the lists that will be used by the streamer.

    The local ids of the cells are different on the differen prrocessors,
    what remains constant is the sequential ghost index. The exchange
    lists for the data communicator contains the ghost index, in this way
    it is possible to communicate these lists to all the processors without
    remapping. The lists tht will be passed to the streamers are generated
    from the data communicator lists, rempping ghost index to local cell
    ids.
*/
void GhostCommunicator::createStreamableLists()
{
    m_sendListIds = m_sendList;
    remapList(m_sendListIds, m_patch->getGhostExchangeSources());

    m_recvListIds = m_recvList;
    remapList(m_recvListIds, m_patch->getGhostExchangeTargets());
}

/*!
    Gets a constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be received from
    \param writer is the streamer that will write the buffer
    \result A constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.
*/
const GhostCommunicator::RankExchangeList & GhostCommunicator::getStreamableRecvList(int rank, ExchangeBufferStreamer *writer)
{
    BITPIT_UNUSED(writer);

    return m_recvListIds.at(rank);
}




/*!
    \class PointListCommunicator

    \brief The PointListCommunicator class provides the infrastructure needed to
    exchange a list of data among processors.
*/

/*!
    Creates a new communicator for data exchange

    \param patch is the patch
*/
PointListCommunicator::PointListCommunicator(const MPI_Comm &communicator)
    : DataCommunicator(communicator),
      m_itemSize(0)
{
}

/*!
    Destructor.
*/
PointListCommunicator::~PointListCommunicator()
{
}

/*!
    Gets the size, expressed in bytes, of the single list item to exchange.

    \result The size, expressed in bytes, of the single list item to exchange.
*/
size_t PointListCommunicator::getItemSize() const
{
    return m_itemSize;
}

/*!
    Gets the list of ids for which the current processor will send data.

    \result The list of ids for which the current processor will send data.
*/
const PointListCommunicator::ExchangeList & PointListCommunicator::getSendList() const
{
    return m_sendList;
}

/*!
    Gets the list of ids for which the current processor will send data
    to the specified rank.

    \result The list of ids for which the current processor will send data
    to the specified rank.
*/
const PointListCommunicator::RankExchangeList & PointListCommunicator::getSendList(int rank) const
{
    return m_sendList.at(rank);
}

/*!
    Gets the list of ids from which the current processor expects to receive
    data.

    \result The list of ids from which the current processor expects to receive
    data.
*/
const PointListCommunicator::ExchangeList & PointListCommunicator::getRecvList() const
{
    return m_recvList;
}

/*!
    Gets the list of ids from which the current processor expects to receive
    data from the specified rank.

    \result The list of ids from which the current processor expects to receive
    data from the specified rank.
*/
const PointListCommunicator::RankExchangeList & PointListCommunicator::getRecvList(int rank) const
{
    return m_recvList.at(rank);
}

/*!
    Sets both the send list and the receive list.

    \param sendList is a list of ids the current processor will send data to
    \param recvList is a list of ids the current processor expects to receive
    data from
*/
void PointListCommunicator::setExchangeLists(const ExchangeList &sendList,
                                        const ExchangeList &recvList)
{
    m_sendList = sendList;
    m_recvList = recvList;

    if (hasData()) {
        updateExchangeInfo();
    }
}

/*!
    Sets either the send list or the receive list.

    \param listType controls if the list defines ids to send or to receive
    \param list is a list of ids the current processor expectst to receive
    data from or will send data for
*/
void PointListCommunicator::setExchangeList(ListType listType, const ExchangeList &list)
{
    switch (listType) {
        case LIST_SEND:
            m_sendList = list;
            m_recvList = scatterExchangeList(m_sendList);
            break;

        case LIST_RECV:
            m_recvList = list;
            m_sendList = scatterExchangeList(m_recvList);
            break;
    }

    if (hasData()) {
        updateExchangeInfo();
    }
}

/*!
    Communicate an exchange list among the other processors

    \param inputList is the list the will be communicated
    \result The exchange list recevied from the other processors.

*/
PointListCommunicator::ExchangeList PointListCommunicator::scatterExchangeList(const ExchangeList &inputList)
{
    // Clear all requests
    //
    // Before setting new requests it is important to cancel all the old ones.
    // Otherwise, if the continuous receive feature is enabled, the new sends
    // may start sending data to the old receives. The synchronous version
    // of clear functions will be used.
    //
    // When using the synchronous version of the clear functions, send and
    // receive requests have to match. It's not possible to first clear the
    // sends synchronously and then clear the receives synchronously, when
    // clearing the receives the sends would have already been cleared and
    // thus the send/receive requests would not match. We need to first
    // cancel the requests and then clear them.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    //
    // Send the list of ids
    //

    // Set the sends for exchanging the ids
    for (const auto &entry : inputList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * sizeof(long);

        setSend(rank, bufferSize);
    }

    // Discover the receives for exchanging the ids
    discoverRecvs();

    // Start the receives
    if (!areRecvsContinuous()) {
        startAllRecvs();
    }

    // Send the list of ids
    for (int rank : getSendRanks()) {
    	bitpit::SendBuffer buffer = getSendBuffer(rank);
        for (const long id : inputList.at(rank)) {
            buffer << id;
        }

        startSend(rank);
    }

    //
    // Get the list of ids
    //

    // Receive the ids
    int nPendingRecvs = getRecvCount();

    ExchangeList outputList;
    outputList.reserve(nPendingRecvs);
    while (nPendingRecvs != 0) {
        int rank = waitAnyRecv();
        bitpit::RecvBuffer buffer = getRecvBuffer(rank);
        long rankListSize = buffer.getSize() / sizeof(long);

        RankExchangeList &rankList = outputList[rank];
        rankList.resize(rankListSize);
        for (size_t i = 0; i < rankList.size(); ++i) {
            buffer >> rankList[i];
        }

        --nPendingRecvs;
    }

    // Wait all sends
    waitAllSends();

    // Clear the requests used for scattering the list
    //
    // Also here the synchronous version of clear functions will be used.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    return outputList;
}

/*!
    Updates the information needed for data exchange
*/
void PointListCommunicator::updateExchangeInfo()
{
    // Clear all requests
    //
    // Before setting new requests it is important to cancel all the old ones.
    // Otherwise, if the continuous receive feature is enabled, the new sends
    // may start sending data to the old receives. The synchronous version
    // of clear functions will be used.
    //
    // When using the synchronous version of the clear functions, send and
    // receive requests have to match. It's not possible to first clear the
    // sends synchronously and then clear the receives synchronously, when
    // clearing the receives the sends would have already been cleared and
    // thus the send/receive requests would not match. We need to first
    // cancel the requests and then clear them.
    cancelAllSends(true);
    cancelAllRecvs(true);

    clearAllSends();
    clearAllRecvs();

    // Update send info
    for (const auto &entry : m_sendList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * m_itemSize;

        setSend(rank, bufferSize);
    }

    // Update receive info
    for (const auto &entry : m_recvList) {
        const int rank = entry.first;
        const auto &list = entry.second;

        long bufferSize = list.size() * m_itemSize;

        setRecv(rank, bufferSize);
    }
}

/*!
    Checks if the communicator has some data to communicate.

    \result Returns true if the communicator has some data to communicate,
    otherwise it returns false.
*/
bool PointListCommunicator::hasData() const
{
    return (m_writers.size() > 0 || m_readers.size() > 0);
}

/*!
    Add the data of the streamer to the communicator.

    \param streamer is the streamer of the data
*/
void PointListCommunicator::addData(ExchangeBufferStreamer *streamer)
{
    addData(streamer, streamer);
}

/*!
    Add the data of the streamers to the communicator.

    \param writer is the streamer that will write the data to the buffer
    \param reader is the streamer that will read the data from the buffer
*/
void PointListCommunicator::addData(ExchangeBufferStreamer *writer, ExchangeBufferStreamer *reader)
{
    // Update the item size
    int writerItemSize = writer->getItemSize();
    int readerItemSize = reader->getItemSize();
    if (readerItemSize != writerItemSize) {
        throw std::runtime_error("The item size of the writer differs from the item size of the reader");
    }

    m_itemSize += writerItemSize;

    // Add the streamers
    m_writers.push_back(writer);
    m_readers.push_back(reader);

    // Update exchange info
    updateExchangeInfo();
}

/*!
    Send ghosts data using non-blocking communications

    \param cellData is the container of the cell data
*/
void PointListCommunicator::startAllExchanges()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    // Start the receives
    for (int rank : getRecvRanks()) {
        if (!isRecvActive(rank)) {
            startRecv(rank);
        }
    }

    // Wait previous sends
    waitAllSends();

    // Fill the buffer with the given field and start sending the data
    for (int rank : getSendRanks()) {
        // Get send buffer
    	bitpit::SendBuffer &buffer = getSendBuffer(rank);

        // Write the buffer
        for (ExchangeBufferStreamer *streamer : m_writers) {
            streamer->write(rank, buffer, getStreamableSendList(rank, streamer));
        }

        // Start the send
        startSend(rank);
    }
}

/*!
    Receive ghosts data using non-blocking communications
*/
void PointListCommunicator::completeAllExchanges()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    completeAllRecvs();

    completeAllSends();
}

/*!
    Receive ghosts data using non-blocking communications
*/
void PointListCommunicator::completeAllRecvs()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    size_t nTotalRecvs = getRecvCount();

    std::vector<int> completedRecvs;
    while (completedRecvs.size() != nTotalRecvs) {
        int rank = completeAnyRecv(completedRecvs);
        completedRecvs.push_back(rank);
    }
}

/*!
    Receive ghost data from any receive and return the associated rank.

    \result The rank of the completed receive or MPI_UNDEFINED if there was
    no active receives.
*/
int PointListCommunicator::completeAnyRecv(const std::vector<int> &blacklist)
{
    // Wait for a receve to finish
    int rank = waitAnyRecv(blacklist);

    // Get receive buffer
    bitpit::RecvBuffer &buffer = getRecvBuffer(rank);

    // Read the buffer
    for (ExchangeBufferStreamer *streamer : m_readers) {
        streamer->read(rank, buffer, getStreamableRecvList(rank, streamer));
    }

    return rank;
}

/*!
    Receive ghosts data using non-blocking communications
*/
void PointListCommunicator::completeAllSends()
{
    if (getCommunicator() == MPI_COMM_NULL || !hasData()) {
        return;
    }

    waitAllSends();
}

/*!
    Remap the exchange lists according to the specified mapper.

    \param mapper is the mapper used for mapping the exchange lists
*/
void PointListCommunicator::remapExchangeLists(const std::unordered_map<long, long> &mapper)
{
    remapRecvList(mapper);
    remapSendList(mapper);
}

/*!
    Remap the exchange lists according to the specified mapper.

    \param mapper is the mapper used for mapping the exchange lists
*/
void PointListCommunicator::remapExchangeLists(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapRecvList(mapper);
    remapSendList(mapper);
}

/*!
    Remap the send list according to the specified mapper.

    \param mapper is the mapper used for mapping the send list
*/
void PointListCommunicator::remapSendList(const std::unordered_map<long, long> &mapper)
{
    remapList(m_sendList, mapper);
}

/*!
    Remap the send list according to the specified mapper.

    \param mapper is the mapper used for mapping the send list
*/
void PointListCommunicator::remapSendList(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapList(m_sendList, mapper);
}

/*!
    Remap the receive list according to the specified mapper.

    \param mapper is the mapper used for mapping the receive list
*/
void PointListCommunicator::remapRecvList(const std::unordered_map<long, long> &mapper)
{
    remapList(m_recvList, mapper);
}

/*!
    Remap the receive list according to the specified mapper.

    \param mapper is the mapper used for mapping the receive list
*/
void PointListCommunicator::remapRecvList(const std::unordered_map<int, std::vector<long>> &mapper)
{
    remapList(m_recvList, mapper);
}

/*!
    Remap the list according to the specified mapper.

    \param list it the list to remap
    \param mapper is the mapper used for mapping the list
*/
void PointListCommunicator::remapList(ExchangeList &list, const std::unordered_map<long, long> &mapper)
{
    for (auto &entry : list) {
        for (long &id : entry.second) {
            id = mapper.at(id);
        }
    }
}

/*!
    Remap the list according to the specified mapper.

    \param list it the list to remap
    \param mapper is the mapper used for mapping the list
*/
void PointListCommunicator::remapList(ExchangeList &list, const std::unordered_map<int, std::vector<long>> &mapper)
{
    for (auto &entry : list) {
        const std::vector<long> &rankMapper = mapper.at(entry.first);
        for (long &id : entry.second) {
            id = rankMapper[id];
        }
    }
}

/*!
    Gets a constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be send to
    \param reader is the streamer that will read the buffer
    \result A constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.
*/
const PointListCommunicator::RankExchangeList & PointListCommunicator::getStreamableSendList(int rank, ExchangeBufferStreamer *reader)
{
    BITPIT_UNUSED(reader);

    return m_sendList.at(rank);
}

/*!
    Gets a constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be received from
    \param writer is the streamer that will write the buffer
    \result A constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.
*/
const PointListCommunicator::RankExchangeList & PointListCommunicator::getStreamableRecvList(int rank, ExchangeBufferStreamer *writer)
{
    BITPIT_UNUSED(writer);

    return m_recvList.at(rank);
}

/*!
    \class PointGhostCommunicator

    \brief The PointGhostCommunicator class provides the infrastructure
    needed to exchange data among ghosts.
*/

/*!
    Creates a new communicator for data exchange

    \param object is the MimmoObject

*/
PointGhostCommunicator::PointGhostCommunicator(const MimmoObject *object)
    : PointListCommunicator(object->getCommunicator()),
      m_object(object), m_patch(object->getPatch())
{
}

/*!
    Resets the exchange list.

    After the exchange list has been reset, data will be exchanged among all
    ghosts.
*/
void PointGhostCommunicator::resetExchangeLists()
{

	//Use the getGhostExchangeSources and getGhostExchangeSources defined in MimmoObject

    const ExchangeList &exchangeSources = m_object->getPointGhostExchangeSources();
    const ExchangeList &exchangeTargets = m_object->getPointGhostExchangeTargets();

    ExchangeList sendList = sequentialIndexesConversion(exchangeSources);
    ExchangeList recvList = sequentialIndexesConversion(exchangeTargets);

    setExchangeLists(sendList, recvList);
}

/*!
    Sets both the send list and the receive list.

    \param sendList is a list of ghosts indexes the current processor will
    send data to
    \param recvList is a list of ghosts indexes the current processor expects
    to receive data from
*/
void PointGhostCommunicator::setExchangeLists(const ExchangeList &sendList,
                                        const ExchangeList &recvList)
{
    // Set the index lists
    PointListCommunicator::setExchangeLists(sendList, recvList);

    // Create lists for the streamers
    createStreamableLists();
}

/*!
    Sets either the send list or the receive list.

    \param listType controls if the list defines indexes to send or to receive
    \param list is a list of ghosts indexes the current processor expectst to
    receive data from or will send data for
*/
void PointGhostCommunicator::setExchangeList(ListType listType, const ExchangeList &list)
{
    // Set the index lists
    PointListCommunicator::setExchangeList(listType, list);

    // Create lists for the streamers
    createStreamableLists();
}

/*!
    Converts a list of ghost ids to a list of sequential ghost indexes.

    \param idList is the list of ghost ids
    \return The list converted to sequential ghost indexes.
*/
PointGhostCommunicator::ExchangeList PointGhostCommunicator::sequentialIndexesConversion(const ExchangeList &idList)
{
    ExchangeList indexList;
    indexList.reserve(idList.size());

    for (const auto &entry : idList) {
        short rank = entry.first;
        const RankExchangeList &rankIdList = entry.second;

        RankExchangeList &rankIndexList = indexList[rank];
        rankIndexList.reserve(rankIdList.size());
        for (size_t i = 0; i < rankIdList.size(); ++i) {
            rankIndexList.push_back(i);
        }
    }

    return indexList;
}

/*!
    Gets a constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be send to
    \param reader is the streamer that will read the buffer
    \result A constant reference to a receive list that can be used by the given
    streamer to send data to the specified rank.
*/
const PointGhostCommunicator::RankExchangeList & PointGhostCommunicator::getStreamableSendList(int rank, ExchangeBufferStreamer *reader)
{
    BITPIT_UNUSED(reader);

    return m_sendListIds.at(rank);
}

/*!
    Creates the lists that will be used by the streamer.

    The local ids of the cells are different on the differen prrocessors,
    what remains constant is the sequential ghost index. The exchange
    lists for the data communicator contains the ghost index, in this way
    it is possible to communicate these lists to all the processors without
    remapping. The lists tht will be passed to the streamers are generated
    from the data communicator lists, rempping ghost index to local cell
    ids.
*/
void PointGhostCommunicator::createStreamableLists()
{
	//Use the getGhostExchangeSources and getGhostExchangeSources defined in MimmoObject

    m_sendListIds = m_sendList;
    remapList(m_sendListIds, m_object->getPointGhostExchangeSources());

    m_recvListIds = m_recvList;
    remapList(m_recvListIds, m_object->getPointGhostExchangeTargets());
}

/*!
    Gets a constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.

    \param rank is the rank data will be received from
    \param writer is the streamer that will write the buffer
    \result A constant reference to a send list that can be used by the given
    streamer to send data to the specified rank.
*/
const PointGhostCommunicator::RankExchangeList & PointGhostCommunicator::getStreamableRecvList(int rank, ExchangeBufferStreamer *writer)
{
    BITPIT_UNUSED(writer);

    return m_recvListIds.at(rank);
}


};

#endif
