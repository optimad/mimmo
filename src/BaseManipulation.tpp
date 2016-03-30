//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PINS METHODS	//
//==================================================//

/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (reference return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (copy return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (pointer return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (reference return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (copy return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (pointer return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};


/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (reference return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (copy return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (pointer return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (reference return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (copy return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

/*!It adds an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (pointer return).
 */
template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (reference return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (copy return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (pointer return).
 * \param[in] setVal Set function of this receiver object (copy argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (reference return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (copy return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (copy argument).
 * \param[in] getVal Get function of this sender object (pointer return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (reference return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (copy return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an input pin to the object.
 * \param[in] objIn Pointer to BaseManipulation sender object.
 * \param[in] getVal Get function of the sender object (pointer return).
 * \param[in] setVal Set function of this receiver object (pointer argument).
 */
template<typename T>
void
BaseManipulation::removePinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	int idx = findPinIn(pin);
	removePinIn(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (reference return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (copy return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};

/*!It removes an output pin to the object.
 * \param[in] objOut Pointer to BaseManipulation receiver object.
 * \param[in] setVal Set function of the receiver object (pointer argument).
 * \param[in] getVal Get function of this sender object (pointer return).
 */
template<typename T>
void
BaseManipulation::removePinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	int idx = findPinOut(pin);
	removePinOut(idx);
	delete pin;
	pin = NULL;
};



//==================================================//
// BASEMANIPULATION CLASS TEMPLATED INPUT METHODS	//
//==================================================//

/*!It sets the input member of the object.
 * \param[in] data Pointer to data to be stored in the input member.
 */
template<typename T>
void
BaseManipulation::setInput(T* data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
	m_input = std::move(dummy);
}

/*!It sets the input member of the object.
 * \param[in] data Data to be stored in the input member.
 */
template<typename T>
void
BaseManipulation::setInput(T& data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(data));
	m_input = std::move(dummy);
}

/*!It gets the input member of the object.
 * \return Pointer to data stored in the input member.
 */
template<typename T>
T*
BaseManipulation::getInput(){
	return(static_cast<IODataT<T>*>(m_input.get())->getData());
}

//==================================================//
// BASEMANIPULATION CLASS TEMPLATED RESULT METHODS	//
//==================================================//

/*!It sets the result member of the object.
 * \param[in] data Pointer to data to be stored in the result member.
 */
template<typename T>
void
BaseManipulation::setResult(T* data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
	m_result = std::move(dummy);
}

/*!It sets the result member of the object.
 * \param[in] data Data to be stored in the result member.
 */
template<typename T>
void
BaseManipulation::setResult(T& data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(data));
	m_result = std::move(dummy);
}

/*!It gets the result member of the object.
 * \return Pointer to data stored in the result member.
 */
template<typename T>
T*
BaseManipulation::getResult(){
	return(static_cast<IODataT<T>*>(m_result.get())->getData());
}

