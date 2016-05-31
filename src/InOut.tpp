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

template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(T *var_){
	m_var_ = var_;
};

template<typename T, typename O>
mimmo::PortOutT<T,O>::PortOutT(O* obj_, T (O::*getVar_)()){
	m_obj_ = obj_;
	m_getVar_ = getVar_;
};


/*!Default destructor of PortOutT
 */
template<typename T, typename O>
mimmo::PortOutT<T,O>::~PortOutT(){
	m_var_ = NULL;
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
	m_var_ = other.m_var_;
	return (*this);
};

/*!Compare operator of PortOutT.
 */
template<typename T, typename O>
bool mimmo::PortOutT<T,O>::operator==(const PortOutT<T,O> & other){
	bool equal = true;
	equal &= (m_var_ == other.m_var_);
	return (equal);
};


template<typename T, typename O>
void
mimmo::PortOutT<T,O>::writeBuffer(){
	if (m_getVar_ != NULL){
		T temp = ((m_obj_->*m_getVar_)());
		m_obuffer << temp;
		return;
	}
	m_obuffer << (*m_var_);
	return;

}

/*!Default constructor of PortInT
 */
template<typename T>
mimmo::PortInT<T>::PortInT(){
	m_var_ = NULL;
};

template<typename T>
mimmo::PortInT<T>::PortInT(T *var_){
	m_var_ = var_;
};


/*!Default destructor of PortInT
 */
template<typename T>
mimmo::PortInT<T>::~PortInT(){
	m_var_ = NULL;
};

/*!Copy constructor of PortInT.
 */
template<typename T>
mimmo::PortInT<T>::PortInT(const PortInT<T> & other){
	*this = other;
};

/*!Assignement operator of PortInT.
 */
template<typename T>
mimmo::PortInT<T> & mimmo::PortInT<T>::operator=(const PortInT<T> & other){
	m_var_ = other.m_var_;
	return (*this);
};

/*!Compare operator of PortInT.
 */
template<typename T>
bool mimmo::PortInT<T>::operator==(const PortInT<T> & other){
	bool equal = true;
	equal &= (m_var_ == other.m_var_);
	return (equal);
};

template<typename T>
void
mimmo::PortInT<T>::readBuffer(){
	T temp;
	m_ibuffer >> temp;
	(*m_var_) = temp;
}

