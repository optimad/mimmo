//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PORTS METHODS	//
//==================================================//

/*!It adds an output port to the object.
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] getVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] label Label (TAG) of the port. The plurality of the same tag is not allowed.
 * \param[in] portS ID of the port; the port will be created at portS-th slot of the output ports of the object.
 */
template<typename T, typename O>
bool
BaseManipulation::createPortOut(O* obj_, T (O::*getVar_)(), PortType label, PortID portS){
	bool check = false;
	if (m_portOut.count(portS) != 0 ) return (check);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(obj_, getVar_);
	portOut->m_label = label;
	m_mapPortOut[label] = portS;
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}


/*!It adds an input port to the object.
 * \param[in] var_ Pointer to member to fill.
 * \param[in] label Label (TAG) of the port. The plurality of the same tag is not allowed.
 * \param[in] portR ID of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] compatibilities Labels (TAG) compatible with the input port. The output ports
 * with tags in compatibilities can send data to this input port. The compatibility between ports with same tag
 * is automatic.
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(T* var_, PortType label, PortID portR, std::vector<PortType> compatibilities){
	bool check = false;
	if (m_portIn.count(portR) != 0 ) return (check);
	PortInT<T, O>* portIn = new PortInT<T, O>(var_);
	m_mapPortIn[label] = portR;
	portIn->addCompatibility(label);
	for(int i = 0; i < compatibilities.size(); i++) {
		portIn->addCompatibility(compatibilities[i]);
	}
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}


/*!It adds an input port to the object.
 * \param[in] obj_ Pointer to BaseManipulation object.
 * \param[in] setVar_ Get function of the object (copy return) linked by the sender port.
 * \param[in] label Label (TAG) of the port. The plurality of the same tag is not allowed.
 * \param[in] portR ID of the port; the port will be created at portR-th slot of the input ports of the object.
 * \param[in] compatibilities Labels (TAG) compatible with the input port. The output ports
 * with tags in compatibilities can send data to this input port. The compatibility between ports with same tag
 * is automatic.
 */
template<typename T, typename O>
bool
BaseManipulation::createPortIn(O* obj_, void (O::*setVar_)(T), PortType label, PortID portR, std::vector<PortType> compatibilities){
	bool check = false;
	if (m_portIn.count(portR) != 0 ) return (check);
	PortInT<T, O>* portIn = new PortInT<T, O>(obj_, setVar_);
	m_mapPortIn[label] = portR;
	portIn->addCompatibility(label);
	for(int i = 0; i < compatibilities.size(); i++) {
		portIn->addCompatibility(compatibilities[i]);
	}
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}


