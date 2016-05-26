//==============================================================//
// TEMPLATE DERIVED INOUT CLASS	TEMPLATE METHODS				//
//==============================================================//

#include <iostream>

/*!Default constructor of PinOutT
 */
template<typename T>
mimmo::PinOutT<T>::PinOutT(){
	m_var_ = NULL;
};

template<typename T>
mimmo::PinOutT<T>::PinOutT(T *var_){
	m_var_ = var_;
};


/*!Default destructor of PinOutT
 */
template<typename T>
mimmo::PinOutT<T>::~PinOutT(){
	m_var_ = NULL;
};

/*!Copy constructor of PinOutT.
 */
template<typename T>
mimmo::PinOutT<T>::PinOutT(const PinOutT<T> & other){
	*this = other;
};

/*!Assignement operator of PinOutT.
 */
template<typename T>
mimmo::PinOutT<T> & mimmo::PinOutT<T>::operator=(const PinOutT<T> & other){
	m_var_ = other.m_var_;
	return (*this);
};

/*!Compare operator of PinOutT.
 */
template<typename T>
bool mimmo::PinOutT<T>::operator==(const PinOutT<T> & other){
	bool equal = true;
	equal &= (m_var_ == other.m_var_);
	return (equal);
};

template<typename T>
void
mimmo::PinOutT<T>::writeBuffer(){
	m_obuffer << (*m_var_);
}

/*!Default constructor of PinInT
 */
template<typename T>
mimmo::PinInT<T>::PinInT(){
	m_var_ = NULL;
};

template<typename T>
mimmo::PinInT<T>::PinInT(T *var_){
	m_var_ = var_;
};


/*!Default destructor of PinInT
 */
template<typename T>
mimmo::PinInT<T>::~PinInT(){
	m_var_ = NULL;
};

/*!Copy constructor of PinInT.
 */
template<typename T>
mimmo::PinInT<T>::PinInT(const PinInT<T> & other){
	*this = other;
};

/*!Assignement operator of PinInT.
 */
template<typename T>
mimmo::PinInT<T> & mimmo::PinInT<T>::operator=(const PinInT<T> & other){
	m_var_ = other.m_var_;
	return (*this);
};

/*!Compare operator of PinInT.
 */
template<typename T>
bool mimmo::PinInT<T>::operator==(const PinInT<T> & other){
	bool equal = true;
	equal &= (m_var_ == other.m_var_);
	return (equal);
};

template<typename T>
void
mimmo::PinInT<T>::readBuffer(){
	m_ibuffer >> (*m_var_);
}

