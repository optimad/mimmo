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
#include "IOConnections.hpp"

namespace mimmo{

/*!
 * Default constructor of IOConnections_MIMMO.
 * Needs a map of executable objects actually instantiated externally, with a string identifying the object.
 * Such strings identifies the object in the XML section from/to which reading/writing also.
 * When builds itself, the class evaluates the inverse map of executable objects and compile the maps of 
 * PortType-s actually available in the API.
 * \param[in]  mapConn map of instantiated connectable objects with their XML names as keys.
 */
IOConnections_MIMMO::IOConnections_MIMMO(std::unordered_map<std::string, BaseManipulation * > mapConn) {

    m_log = &bitpit::log::cout(MIMMO_LOG_FILE);

    for(auto &val: mapConn){
        if(val.second == NULL) continue;
        m_mapConn[val.first] = val.second;
        m_invMapConn[val.second] =  val.first;
    }

    for(auto enumVal : PortType::_values()){
        m_mapPorts[enumVal._to_string()] = (short int) enumVal._to_integral();
        m_invMapPorts[(short int) enumVal._to_integral()] = enumVal._to_string();
    }
};

/*!
 * Default destructor
 */
IOConnections_MIMMO::~IOConnections_MIMMO(){};

/*!
 * Copy constructor
 */
IOConnections_MIMMO::IOConnections_MIMMO(const IOConnections_MIMMO & other){
    *this = other;
};

/*!
 * Assignement operator
 */ 
IOConnections_MIMMO & IOConnections_MIMMO::operator=(const IOConnections_MIMMO & other){
    m_mapConn		= other.m_mapConn;
    m_mapPorts 		= other.m_mapPorts;
    m_invMapConn	= other.m_invMapConn;
    m_invMapPorts 	= other.m_invMapPorts;
    return (*this);
};

/*!
 * \return the map of all ports available in the mimmo API as string key (name of the ports)
 * and corrispective special enum PortType (in its integer form) as argument.
 */
std::unordered_map<std::string, short int>
IOConnections_MIMMO::getMapPorts(){
    return m_mapPorts;
};

/*!
 * \return the inverse map of all ports available in the mimmo API as a PortType enum 
 * (in its integer form) as key and the name (as string) of the port as argument.
 */
std::unordered_map<short int, std::string >
IOConnections_MIMMO::getInvMapPorts(){
    return m_invMapPorts;
};

/*!
 * Reading declared connections from a bitpit XML slot Configuration Section and apply them. Map of 
 * instantiated objects helps class to recognize which objects are already istantiated and ready 
 * to be eventually connected. If the objects map reports a not instantiated object, connection
 * will be refused. The same applies for not existent or not compatible kind of connection ports,
 * declared in the XML section. The method is meant to read connections from xml file.
 *
 * \param[in] slotXML 	Config::Section from which read infos. 
 * \param[in] debug		true, activate verbose debug mode. Default option false.	
 */
void
IOConnections_MIMMO::absorbConnections(const bitpit::Config & slotXML, bool debug){

    if(debug)
        m_log->setPriority(bitpit::log::NORMAL);

    if(slotXML.getSectionCount() == 0)	{
        (*m_log)<<"IOConnections::absorbConnections does not found any connection to read in current XML slot."<<std::endl;
        return;
    }	

    for( auto & sect : slotXML.getSections()){

        std::string snd_str; 
        std::string rcv_str;
        std::string sndP_str;
        std::string rcvP_str;

        if(sect.second->hasOption("sender"))		snd_str = sect.second->get("sender");
        if(sect.second->hasOption("receiver"))		rcv_str = sect.second->get("receiver");
        if(sect.second->hasOption("senderPort"))	sndP_str = sect.second->get("senderPort");
        if(sect.second->hasOption("receiverPort"))	rcvP_str = sect.second->get("receiverPort");

        snd_str = bitpit::utils::trim(snd_str);
        rcv_str = bitpit::utils::trim(rcv_str);
        sndP_str = bitpit::utils::trim(sndP_str);
        rcvP_str = bitpit::utils::trim(rcvP_str);


        auto itSend = m_mapConn.find(snd_str);
        auto itRece = m_mapConn.find(rcv_str);
        auto itSPort = m_mapPorts.find(sndP_str);
        auto itRPort = m_mapPorts.find(rcvP_str);

        if(itSend == m_mapConn.end() || itRece == m_mapConn.end() || itSPort == m_mapPorts.end() || itRPort == m_mapPorts.end()){
            (*m_log)<<"---------------------------------------------"<<std::endl;
            (*m_log)<<"sender: "<<snd_str<<std::endl;
            (*m_log)<<"receiver: "<<rcv_str<<std::endl;
            (*m_log)<<"senderPort: "<<sndP_str<<std::endl;
            (*m_log)<<"receiverPort: "<<rcvP_str<<std::endl;
            (*m_log)<<"---------------------------------------------"<<std::endl;
            (*m_log)<<""<<std::endl;
            (*m_log)<<"IOConnections::absorbConnections does not found any sender/receiver and/or port compatible with API."<<std::endl;
            (*m_log)<<"Please control connection data entries in your XML file"<<std::endl;
            (*m_log)<<""<<std::endl;
            continue;
        }

        bool check = pin::addPin(itSend->second, itRece->second, itSPort->second, itRPort->second);

        if(!check){
            (*m_log)<<"---------------------------------------------"<<std::endl;
            (*m_log)<<"sender: "<<snd_str<<std::endl;
            (*m_log)<<"receiver: "<<rcv_str<<std::endl;
            (*m_log)<<"senderPort: "<<sndP_str<<std::endl;
            (*m_log)<<"receiverPort: "<<rcvP_str<<std::endl;
            (*m_log)<<"---------------------------------------------"<<std::endl;
            (*m_log)<<""<<std::endl;
            (*m_log)<<"IOConnections::absorbConnections failed creating connection."<<std::endl;
            (*m_log)<<""<<std::endl;
        }	
        else{
            (*m_log)<<"IOConnections::absorbConnections successfully created connection."<<std::endl;
        }
    }
    m_log->setPriority(bitpit::log::DEBUG);

    return;
};

/*!
 * Write declared connections to a bitpit XML slot Configuration Section. The map of 
 * instantiated objects helps class to recognize which objects are are connected and compile correctly
 * their xml field. If the objects map reports a not instantiated object, connection
 * will be refused and not written. The same applies for not existent or not compatible kind of connection ports. 
 * The method is meant to write connections to a xml file.
 *
 * \param[in] slotXML   Config::Section to which write infos. 
 * \param[in] debug     true, activate verbose debug mode. Default option false.
 */
void
IOConnections_MIMMO::flushConnections(bitpit::Config & slotXML, bool debug ){

    int counterGlob = 0;
    for( auto object : m_mapConn){

        std::string sender = m_invMapConn[object.second]; //recover sender name
        auto mapPOut = (object.second)->getPortsOut(); //get its slot of output port


        for(auto & pOut : mapPOut ){

            std::string senderPort = m_invMapPorts[pOut.first]; //gets the Port Type of sender

            std::vector<BaseManipulation * > links = (pOut.second)->getLink(); //get all possible links to this specific port
            std::vector<short int>			lids = (pOut.second)->getPortLink();

            int counter = 0;
            for( auto & val : links){
                std::string receiver = m_invMapConn[val];
                std::string receiverPort = m_invMapPorts[lids[counter]];

                std::string title = "conn" + std::to_string(counter+counterGlob);
                bitpit::Config::Section & conn = slotXML.addSection(title);
                conn.set("sender", sender);
                conn.set("senderPort", senderPort);
                conn.set("receiver", receiver);
                conn.set("receiverPort", receiverPort);

                ++counter;
            }//loop on possible links on the same port	
            counterGlob += counter;
        }//loop on all the port of the single object
    }//loop on all instantiated object.

    if(debug)
        m_log->setPriority(bitpit::log::NORMAL);
    (*m_log)<<"IOConnections::flushConnections wrote "<< slotXML.getSectionCount()<< " active connections."<<std::endl;
    m_log->setPriority(bitpit::log::DEBUG);
    return;
};

}
