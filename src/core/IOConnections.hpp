/*---------------------------------------------------------------------------*\
*
*  mimmo
*
*  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#ifndef __IOCONNECTIONS_HPP__
#define __IOCONNECTIONS_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class IOConnections_MIMMO
 * \brief IOConnections_MIMMO is a class to read/write mimmo pin connections from XML IO/parser defined in bitpit::Config.
 * \ingroup core
 *
 * Read and write declared connections from/to an external input/output XML configuration file. Needs to be constructed
 * with an external map reporting a string with "name of the object in XML file" as key, and a "BaseManipulation pointer"
 * of all the object instantiated and connectable as argument. Search the right couple of object to connect and connects
 * them, in input mode, or write existing connections in output mode.
 *
 * In the XML interface a connection between two executable blocks must
   be declared following the nomenclature:

<tt><B>\< name_of_connection\></B> \n
&nbsp;&nbsp;&nbsp;<B>\<sender\></B> name of the sender <B>\</sender\></B> \n
&nbsp;&nbsp;&nbsp;<B>\<senderPort\></B> type of the sender Port <B>\</senderPort\></B> \n
&nbsp;&nbsp;&nbsp;<B>\<receiver\></B> name of the receiver <B>\</receiver\></B> \n
&nbsp;&nbsp;&nbsp;<B>\<receiverPort\></B> name of the receiver Port <B>\</receiverPort\></B> \n
<B>\</name_of_connection\></B> \n</tt>

 * Multiple connections are declared repeating this block format as many times as necessary
 */
class IOConnections_MIMMO{

protected:
    std::unordered_map<std::string, BaseManipulation * > m_mapConn;         /**< direct map of connectable object */
    std::unordered_map<BaseManipulation *, std::string > m_invMapConn;      /**< inverse map of connectable object */

    bitpit::Logger*             m_log;             /**<Pointer to logger.*/

public:
    IOConnections_MIMMO(std::unordered_map<std::string, BaseManipulation * > mapConn);
    virtual ~IOConnections_MIMMO();

    IOConnections_MIMMO(const IOConnections_MIMMO & other);

    //execution utils
    void 	absorbConnections(const bitpit::Config & slotXML, bool debug = false);
    void 	flushConnections(bitpit::Config & slotXML, bool debug = false);
};

};

#endif /* __IOCONNECTIONS_HPP__ */
