/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __IOCONNECTIONS_HPP__
#define __IOCONNECTIONS_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 *	\class IOConnections_MIMMO
 *	\brief IOConnections_MIMMO is a class to read/write mimmo pin connections from XML IO/parser defined in bitpit::Config 
 *
 * Read and write declared connections from/to an external input/output XML configuration file. Needs to be constructed
 * with an external map reporting a string with "name of the object in XML file" as key, and a "BaseManipulation pointer" 
 * of all the object instantiated and connectable as argument. Search the right couple of object to connect and connects them, in input mode, or write 
 * existing connections in output mode.
 */
class IOConnections_MIMMO{

protected:	
	std::unordered_map<std::string, BaseManipulation * > m_mapConn; /**< direct map of connectable object */
	std::unordered_map<BaseManipulation *, std::string > m_invMapConn; /**< inverse map of connectable object */
	std::unordered_map<std::string, short int> m_mapPorts;	/**< map of ports available in mimmo */
	std::unordered_map<short int, std::string> m_invMapPorts;	/**< inverse map of ports available in mimmo */

public:
	IOConnections_MIMMO(std::unordered_map<std::string, BaseManipulation * > mapConn);
	virtual ~IOConnections_MIMMO();

	IOConnections_MIMMO(const IOConnections_MIMMO & other);
	IOConnections_MIMMO & operator=(const IOConnections_MIMMO & other);

	//get methods
	std::unordered_map<std::string, short int> getMapPorts();
	std::unordered_map<short int, std::string> getInvMapPorts();
	
	//execution utils
	void 	absorbConnections(const bitpit::Config & slotXML, bool debug = false);
	void 	flushConnections(bitpit::Config & slotXML, bool debug = false);
};

};

#endif /* __IOCONNECTIONS_HPP__ */
