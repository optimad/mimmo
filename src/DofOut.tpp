//==============================================================//
// TEMPLATE DERIVED DOFOUT CLASS TEMPLATE METHODS				//
//==============================================================//

#include <iostream>

/*!Default constructor of DofOutT
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O>::DofOutT(){
	m_getVal 	= NULL;
	m_getValP 	= NULL;
};

/*!Custom constructor of DofOutT
 * \param[in] obj Pointer to linked object.
 * \param[in] fget Get function of the object (copy return).
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O>::DofOutT(O* obj, T (C::*fget) ()){
	m_obj		=  static_cast<BaseManipulation*>(obj);
	setGetFunction(fget, obj);
}

/*!Custom constructor of DofOutT
 * \param[in] obj Pointer to linked object.
 * \param[in] fget Get function of the object (pointer return).
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O>::DofOutT(O* obj, T* (C::*fget) ()){
	m_obj		= static_cast<BaseManipulation*>(obj);
	setGetFunction(fget, obj);
}


/*!Default destructor of DofOutT
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O>::~DofOutT(){
	m_getVal 	= NULL;
	m_getValP 	= NULL;
};

/*!Copy constructor of DofOutT.
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O>::DofOutT(const DofOutT<T,C,O> & other){
	*this = other;
};

/*!Assignement operator of DofOutT.
 */
template<typename T, typename C, typename O>
mimmo::DofOutT<T,C,O> & mimmo::DofOutT<T,C,O>::operator=(const DofOutT<T,C,O> & other){
	//inherited members
	this->m_obj 	= other.m_obj;
	//its own members
	this->m_getVal 	= other.m_getVal;
	this->m_getValP = other.m_getValP;
	return (*this);
};

/*!Compare operator of DofOutT.
 */
template<typename T, typename C, typename O>
bool mimmo::DofOutT<T,C,O>::operator==(const DofOutT<T,C,O> & other){
	bool equal = true;
	equal &= (m_obj == other.m_obj);
	equal &= (m_getVal == other.m_getVal);
	equal &= (m_getValP == other.m_getValP);
	return (equal);
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] getVal Bound get function of the parent object (copy return).
 */
template<typename T, typename C, typename O>
void
mimmo::DofOutT<T,C,O>::setGetFunction(T (C::*fget) (), O* obj){
	std::function<T(void)> res = std::bind(fget, obj);
	m_getVal	= res;
	m_getValP	= NULL;
};

/*!It sets an input pin for the owner of this pin.
 * \param[in] getValP Bound get function of the parent object (pointer return).
 */
template<typename T, typename C, typename O>
void
mimmo::DofOutT<T,C,O>::setGetFunction(T* (C::*fget) (), O* obj){
	std::function<T*(void)> res = std::bind(fget, obj);
	m_getValP	= res;
	m_getVal	= NULL;
};

template<typename T, typename C, typename O>
T
mimmo::DofOutT<T,C,O>::get(){
	return (m_getVal());
};


template<typename T, typename C, typename O>
T*
mimmo::DofOutT<T,C,O>::getP(){
	return (m_getValP());
};

template<typename T, typename C, typename O>
dvector1D
mimmo::DofOutT<T,C,O>::getDof(){
	T value;
	if (m_getVal != NULL){
		value = m_getVal();
	}
	if (m_getValP != NULL){
		value = (*m_getValP());
	}
	dvector1D dofs = getValue(value);
	dvector1D res;
	for (int i=0; i<dofs.size(); i++){
		if (m_actives[i]) res.push_back(dofs[i]);
	}
	return res;
}

template<typename T, typename C, typename O>
dvector1D
mimmo::DofOutT<T,C,O>::getValue(T value){
	dvector1D res, vval(value.size());
	for (int i=0; i<value.size(); i++){
		vval = getValue(value[i]);
		res.insert(res.end(), vval.begin(), vval.end());
	}
	return res;
}


template<typename T, typename C, typename O>
dvector1D
mimmo::DofOutT<T,C,O>::getValue(double value){
	dvector1D res(1,value);
	return res;
}


template<typename T, typename C, typename O>
dvector1D
mimmo::DofOutT<T,C,O>::getValue(darray3E value){
	dvector1D res(3);
	res[0] = value[0];
	res[1] = value[1];
	res[2] = value[2];
	return res;
}

