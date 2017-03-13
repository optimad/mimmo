//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PORTS METHODS	//
//==================================================//

namespace mimmo{
	
/*!It adds an output port to the object.The input ports with compatible tags can receive data 
 * form this output port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] var_ Pointer to T member to communicate.
 * \param[in] portS ID of the port; the port will be created at portS-th slot of the output ports of the object.
 * \param[in] conType port container type identified through containerTAG enum
 * \param[in] dataType port data type identified through dataTAG enum
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(T* var_, PortID portS, containerTAG conType, dataTAG dataType ){
	bool check = false;
	if (m_portOut.count(portS) != 0 ) return (check);
	DataType datat(conType, dataType);
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
 * \param[in] portS ID of the port; the port will be created at portS-th slot of the output ports of the object.
 * \param[in] conType port container type identified through containerTAG enum
 * \param[in] dataType port data type identified through dataTAG enum
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(O* obj_, T (O::*getVar_)(), PortID portS, containerTAG conType, dataTAG dataType ){
	bool check = false;
	if (m_portOut.count(portS) != 0 ) return (check);
	DataType datat(conType, dataType);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(obj_, getVar_, datat);
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

/*!It adds an input port to the object. The output ports with compatible tags can send data 
 * to this input port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] var_ Pointer to member to fill.
 * \param[in] portR ID of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] conType port container type identified through containerTAG enum
 * \param[in] dataType port data type identified through dataTAG enum
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(T* var_, PortID portR, containerTAG conType, dataTAG dataType ){
	bool check = false;
	if (m_portIn.count(portR) != 0 ) return (check);
	DataType datat(conType, dataType);
	PortInT<T, O>* portIn = new PortInT<T, O>(var_, datat);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

/*!It adds an input port to the object.The output ports compatible tags can send data 
 * to this input port. The compatibility between ports with same tag is automatic.
 * 
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] setVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] portR ID of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] conType port container type identified through containerTAG enum
 * \param[in] dataType port data type identified through dataTAG enum
 * \return true/false, if port is created
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(O* obj_, void (O::*setVar_)(T), PortID portR, containerTAG conType, dataTAG dataType ){
	bool check = false;
	if (m_portIn.count(portR) != 0 ) return (check);
	DataType datat(conType, dataType);
	PortInT<T, O>* portIn = new PortInT<T, O>(obj_, setVar_, datat);
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}

};