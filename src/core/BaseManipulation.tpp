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
namespace mimmo{
	
/*!It adds an output port to the object.The input ports with compatible tags can receive data 
 * form this output port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] var_ Pointer to T member to communicate.
 * \param[in] portS marker of the port; the port will be created at portS-th slot of the output ports of the object.
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(T* var_, PortID portS){
	bool check = false;
    
    //checking if portS containerTAG and dataTAG are registered in mimmo::PortManager
    
    if(!mimmo::PortManager::instance().containsPort(portS) || m_portOut.count(portS) != 0 ){
        (*m_log)<<"Unable to physically create the Output Port."<<std::endl;
        (*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
        return check;
    }
    
    auto info = mimmo::PortManager::instance().getPortData(portS);
	DataType datat(info.container, info.datatype);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(var_, datat);
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

/*!It adds an output port to the object.The input ports with compatible tags can receive data 
 * form this output port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] getVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] portS marker of the port; the port will be created at portS-th slot of the output ports of the object.
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(O* obj_, T (O::*getVar_)(), PortID portS){
	bool check = false;

    //checking if portS containerTAG and dataTAG are registered in mimmo::PortManager
  
    if(!mimmo::PortManager::instance().containsPort(portS) || m_portOut.count(portS) != 0 ){
        (*m_log)<<"Unable to physically create the Output Port."<<std::endl;
        (*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
        return check;
    }
    
    auto info = mimmo::PortManager::instance().getPortData(portS);    
    DataType datat(info.container, info.datatype);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(obj_, getVar_, datat);
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

/*!It adds an input port to the object. The output ports with compatible container/data tags can send data 
 * to this input port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] var_ Pointer to member to fill.
 * \param[in] portR marker of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] mandatory does the port have to be mandatorily linked?
 * \param[in] family tag of family of mandatory ports (family=0 independent ports with no alternative).
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(T* var_, PortID portR, bool mandatory, int family ){
	bool check = false;
	
    //checking if portR containerTAG and dataTAG are registered in mimmo::PortManager
    
    if(!mimmo::PortManager::instance().containsPort(portR) || m_portOut.count(portR) != 0 ){
        (*m_log)<<"Unable to physically create the Input Port."<<std::endl;
        (*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
        return check;
    }
    
    auto info = mimmo::PortManager::instance().getPortData(portR);
	DataType datat(info.container, info.datatype);
	PortInT<T, O>* portIn = new PortInT<T, O>(var_, datat, mandatory, family);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

/*!It adds an input port to the object.The output ports compatible tags can send data 
 * to this input port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] setVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] portR marker of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] mandatory does the port have to be mandatorily linked?
 * \param[in] family tag of family of mandatory ports (family=0 independent ports with no alternative).
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(O* obj_, void (O::*setVar_)(T), PortID portR, bool mandatory, int family ){
	bool check = false;
    
    //checking if portR containerTAG and dataTAG are registered in mimmo::PortManager
    if(!mimmo::PortManager::instance().containsPort(portR) || m_portOut.count(portR) != 0 ){
        (*m_log)<<"Unable to physically create the Input Port."<<std::endl;
        (*m_log)<<"Port name was not regularly registered or a port with the same name was already instantiated in the current class."<<std::endl;
        return check;
    }
    
    auto info = mimmo::PortManager::instance().getPortData(portR);
	DataType datat(info.container, info.datatype);
	PortInT<T, O>* portIn = new PortInT<T, O>(obj_, setVar_, datat, mandatory, family);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

};
