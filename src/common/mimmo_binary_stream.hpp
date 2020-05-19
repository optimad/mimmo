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
\*---------------------------------------------------------------------------*/
#ifndef __MIMMOBINARYSTREAM_HPP__
#define __MIMMOBINARYSTREAM_HPP__

#include <binary_stream.hpp>
#include <array>
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>

namespace mimmo{

 /*!
     @class IBinaryStream
     @ingroup binaryStream
     @brief mimmo custom derivation of bitpit IBinaryStream (see relative doc)
 */
 class IBinaryStream : public bitpit::IBinaryStream {

 public:
     IBinaryStream(void);
     IBinaryStream(std::size_t size);
     IBinaryStream(const char *buffer, std::size_t size);
     IBinaryStream(const std::vector<char> &buffer);
 };

 /*!
     @class OBinaryStream
     @ingroup binaryStream
     @brief mimmo custom derivation of bitpit OBinaryStream (see relative doc)
 */
 class OBinaryStream : public bitpit::OBinaryStream {

 public:
     OBinaryStream();
     OBinaryStream(std::size_t size);
 };


}

/*!
 * \ingroup binaryStream
 * \{
 */
//BASIC TYPES
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, std::string & element);
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::string & element);

//TEMPLATE STRUCTURES

template<typename T>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf,std::vector<T>& element);
template<typename T>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::vector<T>& element);

template<typename T, std::size_t d>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, std::array<T,d>& element);
template<typename T, std::size_t d>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::array<T,d>& element);

template<typename T, typename Q>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, std::pair<T, Q>& element);
template<typename T, typename Q>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::pair<T, Q>& element);

template<typename T, typename Q>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, std::unordered_map<T, Q>&  element);
template<typename T, typename Q>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::unordered_map<T, Q >& element);

template<typename T, typename Q>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, std::map<T,Q>& element);
template<typename T, typename Q>
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const std::map<T, Q>& element);

/*!
 *\}
 */

#include "mimmo_binary_stream.tpp"



#endif /* __MIMMOBINARYSTREAM_HPP__ */
