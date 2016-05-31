//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PINS METHODS	//
//==================================================//

///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (reference return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (copy return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (pointer return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (reference return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (copy return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (pointer return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
//
///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (reference return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (copy return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (pointer return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	m_pinIn.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (reference return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (copy return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
///*!It adds an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (pointer return).
// */
//template<typename T>
//void
//BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	m_pinOut.push_back(pin);
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (reference return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (copy return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (pointer return).
// * \param[in] setVal Set function of this receiver object (copy argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (reference return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (copy return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (copy argument).
// * \param[in] getVal Get function of this sender object (pointer return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (reference return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (copy return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an input pin to the object.
// * \param[in] objIn Pointer to BaseManipulation sender object.
// * \param[in] getVal Get function of the sender object (pointer return).
// * \param[in] setVal Set function of this receiver object (pointer argument).
// */
//template<typename T>
//void
//BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setInput(objIn, getVal, setVal);
//	int idx = findPinIn(pin);
//	removePinIn(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (reference return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (copy return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//
///*!It removes an output pin to the object.
// * \param[in] objOut Pointer to BaseManipulation receiver object.
// * \param[in] setVal Set function of the receiver object (pointer argument).
// * \param[in] getVal Get function of this sender object (pointer return).
// */
//template<typename T>
//void
//BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
//	InOutT<T>* pin = new InOutT<T>();
//	pin->setOutput(objOut, setVal, getVal);
//	int idx = findPinOut(pin);
//	removePinOut(idx);
//	delete pin;
//	pin = NULL;
//};
//

template<typename T, typename O>
bool
BaseManipulation::createPortOut(O* obj, T (O::*getVar_)(), PortType label, PortID portS){
	bool check = false;
	if (getNPortsOut() < portS+1) m_portOut.resize(portS+1);
	if (m_portOut[portS] != NULL ) return (check);
	PortOutT<T, O>* portOut = new PortOutT<T, O>(obj, getVar_);
	portOut->m_label = label;
	m_mapPortOut[label] = portS+1;
	m_portOut[portS] = portOut;
	check = true;
	return(check);
}

template<typename T>
bool
BaseManipulation::createPortIn(T* var_, PortType label, PortID portR, std::vector<PortType> compatibilities){
	bool check = false;
	if (getNPortsIn() < portR+1) m_portIn.resize(portR+1);
	if (m_portIn[portR] != NULL ) return (check);
	PortInT<T>* portIn = new PortInT<T>(var_);
	m_mapPortIn[label] = portR+1;
	portIn->addCompatibility(label);
	for(int i = 0; i < compatibilities.size(); i++) {
		portIn->addCompatibility(compatibilities[i]);
	}
	m_portIn[portR] = portIn;
	check = true;
	return(check);
}


