/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
// TEMPLATE BINARY STREAMS
//==============================================================//

/*!
    Output stream operator for vector<T>
    \param[in] buffer is the output stream
    \param[in] var is the element to be streamed
    \result Returns the same output stream received in input.
*/
template<typename T>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream  &buffer, const std::vector<T> &var)
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
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, std::vector<T> &var)
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
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream  &buffer, const std::array<T,d> &var)
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
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, std::array<T,d> &var)
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
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, std::pair<T, Q>& element){

    T geo;
    Q data;
    buffer >> geo;
    buffer >> data;
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
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const std::pair<T, Q>& element){

    buffer<<element.first;
    buffer<<element.second;
    return buffer;
};


/*!
*	Input stream operator for std::unordered_map\<T,Q\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
template<typename T, typename Q>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer,std::unordered_map<T, Q >&  var){

    T key;
    Q value;
    std::size_t nP;
    buffer >> nP;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> key;
        buffer >> value;
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
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const std::unordered_map<T, Q>& var){

    std::size_t nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee.first;
        buffer <<ee.second;
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
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer,std::map<T, Q >&  var){

    T key;
    Q value;
    std::size_t nP;
    buffer >> nP;
    for (std::size_t i = 0; i < nP; ++i) {
        buffer >> key;
        buffer >> value;
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
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const std::map<T, Q>& var){

    std::size_t nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee.first;
        buffer <<ee.second;
    }
    return buffer;
};
