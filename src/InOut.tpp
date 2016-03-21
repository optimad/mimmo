
//==============================================================//
// TEMPLATE DERIVED INOUT CLASS	TEMPLATE METHODS				//
//==============================================================//

/*!Default constructor of InOutT
 */
template<typename T>
InOutT<T>::InOutT(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_getValP 	= NULL;
	m_setVal 	= NULL;
	m_setValP 	= NULL;
};


/*!Default destructor of InOutT
 */
template<typename T>
InOutT<T>::~InOutT(){
	m_getVal 	= NULL;
	m_getValR 	= NULL;
	m_getValP 	= NULL;
	m_setVal 	= NULL;
	m_setValP 	= NULL;
};

/*!Copy constructor of InOutT.
 */
template<typename T>
InOutT<T>::InOutT(const InOutT<T> & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_getValP 	= other.m_getValP;
	m_setVal 	= other.m_setVal;
	m_setValP 	= other.m_setValP;
};

/*!Assignement operator of InOutT.
 */
template<typename T>
InOutT<T> & InOutT<T>::operator=(const InOutT<T> & other){
	m_getVal 	= other.m_getVal;
	m_getValR 	= other.m_getValR;
	m_getValP 	= other.m_getValP;
	m_setVal 	= other.m_setVal;
	m_setValP 	= other.m_setValP;
	return (*this);
};

/*!Compare operator of InOutT.
 */
template<typename T>
bool InOutT<T>::operator==(const InOutT<T> & other){
	bool equal = true;
	equal &= (*this InOut::== other);
	equal &= (m_getVal == other.m_getVal);
	equal &= (m_getValR == other.m_getValR);
	equal &= (m_getValP == other.m_getValP);
	equal &= (m_setVal == other.m_setVal);
	equal &= (m_setValP == other.m_setValP);
	return (equal);
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (copy return).
 * \param[in] setVal Bound set function of the pin owner object (copy argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (copy argument).
 * \param[in] setVal Bound get function of the pin owner object (copy return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setVal	= setVal;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (reference return).
 * \param[in] setVal Bound set function of the pin owner object (copy argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T&(void)> getValR, std::function<void(T)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (copy argument).
 * \param[in] setVal Bound get function of the pin owner object (reference return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setVal	= setVal;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (pointer return).
 * \param[in] setVal Bound set function of the pin owner object (copy argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T*(void)> getValP, std::function<void(T)> setVal){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValP	= getValP;
	m_setVal	= setVal;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (copy argument).
 * \param[in] setVal Bound get function of the pin owner object (pointer return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getValP){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValP	= getValP;
	m_setVal	= setVal;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (copy return).
 * \param[in] setVal Bound set function of the pin owner object (pointer argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setValP){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getVal	= getVal;
	m_setValP	= setValP;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (pointer argument).
 * \param[in] setVal Bound get function of the pin owner object (copy return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T(void)> getVal){
	m_type 		= true;
	m_objLink	= objOut;
	m_getVal	= getVal;
	m_setValP	= setValP;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (reference return).
 * \param[in] setVal Bound set function of the pin owner object (pointer argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T&(void)> getValR, std::function<void(T*)> setValP){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValR	= getValR;
	m_setValP	= setValP;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (pointer argument).
 * \param[in] setVal Bound get function of the pin owner object (reference return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T&(void)> getValR){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValR	= getValR;
	m_setValP	= setValP;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation parent object to be linked.
 * \param[in] getVal Bound get function of the parent object (pointer return).
 * \param[in] setVal Bound set function of the pin owner object (pointer argument).
 */
template<typename T>
void
InOutT<T>::setInput(BaseManipulation* objIn, std::function<T*(void)> getValP, std::function<void(T*)> setValP){
	m_type 		= false;
	m_objLink 	= objIn;
	m_getValP	= getValP;
	m_setValP	= setValP;
};

/*!It sets an output pin for the owner of this pin.
 * \param[in] objIn Pointer to BaseManipulation child object to be linked.
 * \param[in] getVal Bound set function of the child object (pointer argument).
 * \param[in] setVal Bound get function of the pin owner object (pointer return).
 */
template<typename T>
void
InOutT<T>::setOutput(BaseManipulation* objOut, std::function<void(T*)> setValP, std::function<T*(void)> getValP){
	m_type 		= true;
	m_objLink	= objOut;
	m_getValP	= getValP;
	m_setValP	= setValP;
};

/*! Execution of the pin.
 * During the execution the get/set function of the parent/child objects are called.
 * All the pins of an object are called in execute of the owner after its own exec.
 *
 */
template<typename T>
void
InOutT<T>::exec(){
	if (m_getVal != NULL){
		T val = m_getVal();
		if (m_setVal != NULL){
			m_setVal(val);
		}else if (m_setValP != NULL){
			m_setValP(&val);
		}
	}else if (m_getValR != NULL){
		T& val = m_getValR();
		if (m_setVal != NULL){
			m_setVal(val);
		}else if (m_setValP != NULL){
			m_setValP(&val);
		}
	}else if (m_getValP != NULL){
		T* val = m_getValP();
		if (m_setVal != NULL){
			m_setVal(*val);
		}else if (m_setValP != NULL){
			m_setValP(val);
		}
	}
};
