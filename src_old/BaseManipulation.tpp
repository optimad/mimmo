//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PINS METHODS	//
//==================================================//

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};


template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};


//==================================================//
// BASEMANIPULATION CLASS TEMPLATED INPUT METHODS	//
//==================================================//

template<typename T>
void
BaseManipulation::setInput(T* data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
	m_input = std::move(dummy);
}

template<typename T>
void
BaseManipulation::setInput(T& data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(data));
	m_input = std::move(dummy);
}

template<typename T>
T*
BaseManipulation::getInput(){
	return(static_cast<IODataT<T>*>(m_input.get())->getData());
}


//==================================================//
// BASEMANIPULATION CLASS TEMPLATED RESULT METHODS	//
//==================================================//

template<typename T>
void
BaseManipulation::setResult(T* data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(*data));
	m_result = std::move(dummy);
}

template<typename T>
void
BaseManipulation::setResult(T& data){
	std::unique_ptr<IOData> dummy(new IODataT<T>(data));
	m_result = std::move(dummy);
}

template<typename T>
T*
BaseManipulation::getResult(){
	return(static_cast<IODataT<T>*>(m_result)->getData());
}

