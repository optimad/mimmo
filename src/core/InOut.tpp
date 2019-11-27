/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
 \ *---------------------------------------------------------------------------*/

//==============================================================//
// TEMPLATE DERIVED INOUT CLASS	TEMPLATE METHODS                //
//==============================================================//
#include <binary_stream.hpp>

namespace mimmo {

/*!
 * Default constructor of PortOutT
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(){
    m_var_ = NULL;
};

/*!
 * Custom constructor of PortOutT
 * \param[in] var_ Pointer to variable to be streamed.
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(T *var_){
    m_obj_ 		= NULL;
    m_var_ 		= var_;
    m_getVar_ 	= NULL;
};

/*!
 * Custom constructor of PortOutT
 * \param[in] var_ Pointer to variable to be streamed.
 * \param[in] datatype TAG of datat type communicated.
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(T *var_, DataType datatype){
    m_obj_ 		= NULL;
    m_var_ 		= var_;
    m_getVar_ 	= NULL;
    m_datatype	= datatype;
};

/*!
 * Custom constructor of PortOutT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] getVar_ Pointer to function that gets the data to be streamed.
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(O* obj_, T (O::*getVar_)()){
    m_obj_ 		= obj_;
    m_getVar_ 	= getVar_;
    m_var_ 		= NULL;
};

/*!
 * Custom constructor of PortOutT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] getVar_ Pointer to function that gets the data to be streamed.
 * \param[in] datatype TAG of datat type communicated.
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(O* obj_, T (O::*getVar_)(), DataType datatype){
    m_obj_ 		= obj_;
    m_getVar_ 	= getVar_;
    m_var_ 		= NULL;
    m_datatype	= datatype;
};


/*!
 * Default destructor of PortOutT
 */
template<typename T, typename O>
PortOutT<T,O>::~PortOutT(){
    m_obj_ 		= NULL;
    m_var_ 		= NULL;
    m_getVar_ 	= NULL;
};

/*!
 * Copy constructor of PortOutT.
 */
template<typename T, typename O>
PortOutT<T,O>::PortOutT(const PortOutT<T,O> & other):PortOut(other){
    m_obj_      = other.m_obj_;
    m_var_      = other.m_var_;
    m_getVar_   = other.m_getVar_;
};

/*!
 * Compare operator of PortOutT.
 */
template<typename T, typename O>
bool  PortOutT<T,O>::operator==(const PortOutT<T,O> & other){
    bool equal = true;
    equal &= (m_obj_ == other.m_obj_);
    equal &= (m_var_ == other.m_var_);
    equal &= (m_getVar_ == other.m_getVar_);
    return (equal);
};

/*!
* It writes the buffer of the output port with the data to be communicated.
* It uses the linked get function if the member pointer m_getVar_ is not NULL.
* Alternatively it uses m_var_ directly (if not NULL).
*/
template<typename T, typename O>
void
PortOutT<T,O>::writeBuffer(){
    if (m_getVar_ != NULL){
        T temp = ((m_obj_->*m_getVar_)());
        m_obuffer << temp;
        return;
    }
    if (m_var_ != NULL){
        m_obuffer << (*m_var_);
    }
}



/*!
 * Default constructor of PortInT
 */
template<typename T, typename O>
PortInT<T, O>::PortInT(){
    m_obj_ 		= NULL;
    m_var_ 		= NULL;
    m_setVar_ 	= NULL;
};

/*!
 * Custom constructor of PortInT
 * \param[in] var_ Pointer to variable to be streamed.
 */
template<typename T, typename O>
PortInT<T, O>::PortInT(T *var_){
    m_obj_ 		= NULL;
    m_var_ 		= var_;
    m_setVar_ 	= NULL;
    m_mandatory = false;
};

/*!
 * Custom constructor of PortInT
 * \param[in] var_ Pointer to variable to be streamed.
 * \param[in] datatype TAG of datat type communicated.
 * \param[in] mandatory mandatory port?
 * \param[in] family mandatory family tag.
 */
template<typename T, typename O>
PortInT<T, O>::PortInT(T *var_, DataType datatype, bool mandatory, int family){
    m_obj_ 		= NULL;
    m_var_ 		= var_;
    m_setVar_ 	= NULL;
    m_datatype	= datatype;
    m_mandatory = mandatory;
    m_familym   = family;
};

/*!
 * Custom constructor of PortInT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] setVar_ Pointer to function that sets members with the data in buffer.
 * \param[in] mandatory mandatory port?
 * \param[in] family mandatory family tag.
 */
template<typename T, typename O>
PortInT<T,O>::PortInT(O* obj_, void (O::*setVar_)(T), bool mandatory, int family){
    m_obj_ 		= obj_;
    m_setVar_ 	= setVar_;
    m_var_ 		= NULL;
    m_mandatory = mandatory;
    m_familym   = family;
};

/*!
 * Custom constructor of PortInT
 * \param[in] obj_ Pointer to object owner of the port.
 * \param[in] setVar_ Pointer to function that sets members with the data in buffer.
 * \param[in] datatype TAG of datat type communicated.
 * \param[in] mandatory mandatory port?
 * \param[in] family mandatory family tag.
 */
template<typename T, typename O>
PortInT<T,O>::PortInT(O* obj_, void (O::*setVar_)(T), DataType datatype, bool mandatory, int family){
    m_obj_ 		= obj_;
    m_setVar_ 	= setVar_;
    m_var_ 		= NULL;
    m_datatype	= datatype;
    m_mandatory = mandatory;
    m_familym   = family;
};

/*!
 * Default destructor of PortInT
 */
template<typename T, typename O>
PortInT<T, O>::~PortInT(){
    m_obj_ 		= NULL;
    m_var_ 		= NULL;
    m_setVar_ 	= NULL;
};

/*!
 * Copy constructor of PortInT.
 */
template<typename T, typename O>
PortInT<T, O>::PortInT(const PortInT<T, O> & other):PortIn(other){
    m_obj_      = other.m_obj_;
    m_var_      = other.m_var_;
    m_setVar_   = other.m_setVar_;
};

/*!
 * Compare operator of PortInT.
 */
template<typename T, typename O>
bool  PortInT<T, O>::operator==(const PortInT<T, O> & other){
    bool equal = true;
    equal &= (m_obj_ == other.m_obj_);
    equal &= (m_var_ == other.m_var_);
    equal &= (m_setVar_ == other.m_setVar_);
    return (equal);
};

/*!
 * It reads the buffer of the output port with the data to be communicated.
 * It stores the read values in the linked m_var_ by casting in the stream operator.
 */
template<typename T, typename O>
void
PortInT<T, O>::readBuffer(){
    T temp;
    m_ibuffer >> temp;
    if (m_setVar_ != NULL){
        (m_obj_->*m_setVar_)(temp);
        return;
    }
    if (m_var_ != NULL){
        (*m_var_) = temp;
    }
}


} //end of mimmo namespace

//BINARY STREAMS
/*!
    Output stream operator for vector<T>
    \param[in] buffer is the output stream
    \param[in] var is the element to be streamed
    \result Returns the same output stream received in input.
*/
template<typename T>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const std::vector<T> &var)
{
    std::size_t nP = var.size();
    buffer << nP;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer << var[i];
    }
    return buffer;
}


/*!
    Input stream operator for vector<T>
    \param[in] buffer is the input stream
    \param[in] var is the element to be streamed
    \result Returns the same input stream received in input.
*/
template<typename T>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector<T> &var)
{
    std::size_t nP;
    buffer >> nP;
    var.resize(nP);
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> var[i];
    }
    return buffer;
}


/*!
*	Output stream operator for std::array\<T,d\> enum
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
template <typename T, std::size_t d>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const std::array<T,d> &var)
{
    for(std::size_t i=0; i<d; ++i){
        buffer << var[i];
    }
    return buffer;
}


/*!
*	Input stream operator for std::array\<T,d\> enum
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template <typename T, std::size_t d>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::array<T,d> &var)
{
    for(std::size_t i=0; i<d; ++i){
        buffer >> var[i];
    }
    return buffer;
}

/*!
*	Input stream operator for std::pair\<T,Q\>
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template<typename T, typename Q>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::pair<T, Q>& element){

    T geo;
    Q data;
    buffer >> geo >> data ;
    element = std::make_pair(geo, data);
    return buffer;
};

/*!
*	Output stream operator for std::pair\<T, Q\>
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template<typename T, typename Q>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::pair<T, Q>& element){

    buffer<<element.first<<element.second;
    return buffer;
};


/*!
*	Input stream operator for std::unordered_map\<T,Q\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template<typename T, typename Q>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::unordered_map<T, Q >&  var){

    T key;
    Q value;
    std::size_t nP;
    buffer >> nP;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> key>> value;
        var[key] = value;
    }
    return buffer;
};

/*!
*	Output stream operator for std::unordered_map\<T, Q\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
template<typename T, typename Q>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::unordered_map<T, Q>& var){

    std::size_t nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee.first<<ee.second;
    }
    return buffer;
};

/*!
*	Input stream operator for std::map\<T,Q\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template<typename T, typename Q>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::map<T, Q >&  var){

    T key;
    Q value;
    std::size_t nP;
    buffer >> nP;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> key>> value;
        var[key] = value;
    }
    return buffer;
};

/*!
*	Output stream operator for std::map\<T, Q\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
template<typename T, typename Q>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::map<T, Q>& var){

    std::size_t nP = var.size();
    buffer << nP;
    for (std::pair<T,Q> & ee : var) {
        buffer << ee.first<<ee.second;
    }
    return buffer;
};

/*!
* Input stream operator for bitpit::MimmoPiercedVector\< T \>
* \param[in] buffer is the input stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
template< typename T>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, mimmo::MimmoPiercedVector<T>& element){

    element.clear();
    mimmo::MimmoObject* geo;
    buffer >> geo;
    element.setGeometry(geo);

    std::string name;
    // std::size_t length;
    // buffer >> length;
    // std::vector<char> vectorname(length);
    // for (char & a : vectorname){
    // 	buffer >> a;
    // }
    // std::memcpy(&name, vectorname.data(), length);
    buffer>>name;
    element.setName(name);

    int loc_;
    buffer >> loc_;
    element.setDataLocation(static_cast<mimmo::MPVLocation>(loc_));
    std::size_t nP;
    buffer >> nP;
    T val;
    long int id;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> id;
        buffer >> val;
        element.insert(id, val);
    }
    return buffer;
};

/*!
* Output stream operator for bitpit::MimmoPiercedVector\< T \>
* \param[in] buffer is the output stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
template<typename T>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const mimmo::MimmoPiercedVector<T>& element){

    buffer << element.getGeometry();

    std::string name = element.getName();
    // std::size_t length = name.size();
    // buffer << length;
    // for (std::size_t i=0; i<length; i++){
    // 	char a = name.at(i);
    // 	buffer << a;
    // }
    buffer << name;
    buffer << static_cast<int>(element.getConstDataLocation());
    buffer << (std::size_t)element.size();
    auto itE = element.cend();
    for (auto it=element.cbegin(); it!=itE; it++){
        buffer << it.getId()<<*it;
    }
    return buffer;
};
