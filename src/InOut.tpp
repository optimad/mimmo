//==============================================================//
// TEMPLATE DERIVED INOUT CLASS	TEMPLATE METHODS				//
//==============================================================//

#include <iostream>

/*!Default constructor of PortOutT
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(){
	m_var_ = NULL;
};

/*!Custom constructor of PortOutT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] var_ Pointer to variable to be streamed.
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(T *var_){
	m_obj_ 		= NULL;
	m_var_ 		= var_;
	m_getVar_ 	= NULL;
};

/*!Custom constructor of PortOutT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] getVar_ Pointer to function that gets the data to be streamed.
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(O* obj_, T (O::*getVar_)()){
	m_obj_ 		= obj_;
	m_getVar_ 	= getVar_;
	m_var_ 		= NULL;
};


/*!Default destructor of PortOutT
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::~PortOutT(){
	m_obj_ 		= NULL;
	m_var_ 		= NULL;
	m_getVar_ 	= NULL;
};

/*!Copy constructor of PortOutT.
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(const PortOutT<T,O> & other){
	*this = other;
};

/*!Assignement operator of PortOutT.
 */
template<typename T, typename O>
mimmo::PortOutT<T,O> & mimmo::PortOutT<T,O>::operator=(const PortOutT<T,O> & other){
	m_obj_ 		= other.m_obj_;
	m_var_ 		= other.m_var_;
	m_getVar_ 	= other.m_getVar_;
	return (*this);
};

/*!Compare operator of PortOutT.
 */
template<typename T, typename O>
bool mimmo::PortOutT<T,O>::operator==(const PortOutT<T,O> & other){
	bool equal = true;
	equal &= (m_obj_ == other.m_obj_);
	equal &= (m_var_ == other.m_var_);
	equal &= (m_getVar_ == other.m_getVar_);
	return (equal);
};

/*!It writes the buffer of the output port with the data to be communicated.
 * It uses the linked get function if the member pointer m_getVar_ is not NULL.
 * Alternatively it uses m_var_ directly (if not NULL).
 */
template<typename T, typename O>
void
mimmo::PortOutT<T,O>::writeBuffer(){
	if (m_getVar_ != NULL){
		T temp = ((m_obj_->*m_getVar_)());
		m_obuffer << temp;
		return;
	}
	if (m_var_ != NULL){
		m_obuffer << (*m_var_);
	}
	return;
}



/*!Default constructor of PortInT
 */
template<typename T, typename O>
mimmo::PortInT<T, O>::PortInT(){
	m_obj_ 		= NULL;
	m_var_ 		= NULL;
	m_setVar_ 	= NULL;
};

/*!Custom constructor of PortInT
 * \param[in] var_ Pointer to variable to be streamed.
 */
template<typename T, typename O>
mimmo::PortInT<T, O>::PortInT(T *var_){
	m_obj_ 		= NULL;
	m_var_ 		= var_;
	m_setVar_ 	= NULL;
};

/*!Custom constructor of PortInT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] setVar_ Pointer to function that sets members with the data in buffer.
 */
template<typename T, typename O>
mimmo::PortInT<T,O>::PortInT(O* obj_, void (O::*setVar_)(T)){
	m_obj_ 		= obj_;
	m_setVar_ 	= setVar_;
	m_var_ 		= NULL;
};

/*!Default destructor of PortInT
 */
template<typename T, typename O>
mimmo::PortInT<T, O>::~PortInT(){
	m_obj_ 		= NULL;
	m_var_ 		= NULL;
	m_setVar_ 	= NULL;
};

/*!Copy constructor of PortInT.
 */
template<typename T, typename O>
mimmo::PortInT<T, O>::PortInT(const PortInT<T, O> & other){
	*this = other;
};

/*!Assignement operator of PortInT.
 */
template<typename T, typename O>
mimmo::PortInT<T, O> & mimmo::PortInT<T, O>::operator=(const PortInT<T, O> & other){
	m_obj_ 		= other.m_obj_;
	m_var_ 		= other.m_var_;
	m_setVar_	= other.m_setVar_;
	return (*this);
};

/*!Compare operator of PortInT.
 */
template<typename T, typename O>
bool mimmo::PortInT<T, O>::operator==(const PortInT<T, O> & other){
	bool equal = true;
	equal &= (m_obj_ == other.m_obj_);
	equal &= (m_var_ == other.m_var_);
	equal &= (m_setVar_ == other.m_setVar_);
	return (equal);
};

/*!It reads the buffer of the output port with the data to be communicated.
 * It stores the read values in the linked m_var_ by casting in the stream operator.
 */
template<typename T, typename O>
void
mimmo::PortInT<T, O>::readBuffer(){
	T temp;
	m_ibuffer >> temp;
	if (m_setVar_ != NULL){
		(m_obj_->*m_setVar_)(temp);
		return;
	}
	if (m_var_ != NULL){
		(*m_var_) = temp;
	}
	return;
}

