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
#include "mimmo_binary_stream.hpp"

namespace mimmo{
    /*!
        base constructor
    */
    IBinaryStream::IBinaryStream(void) : bitpit::IBinaryStream(){};
    /*!
        custom constructor
        @param[in] size of the buffer
    */
    IBinaryStream::IBinaryStream(std::size_t size) : bitpit::IBinaryStream(size){};
    /*!
        custom constructor
        @param[in] buffer pointer to buffer
        @param[in] size size of the buffer
    */
    IBinaryStream::IBinaryStream(const char *buffer, std::size_t size) : bitpit::IBinaryStream(buffer,size){};
    /*!
        custom constructor
        @param[in] buffer as vector of char
    */
    IBinaryStream::IBinaryStream(const std::vector<char> &buffer): bitpit::IBinaryStream(buffer){};

    /*!
        base constructor
    */
    OBinaryStream::OBinaryStream(): bitpit::OBinaryStream(){};
    /*!
        custom constructor
        @param[in] size size of the buffer
    */
    OBinaryStream::OBinaryStream(std::size_t size): bitpit::OBinaryStream(size){};

};

/*!
    Output stream operator for std::string
    \param[in] buffer is the output stream
    \param[in] element is the element to be streamed
    \result Returns the same output stream received in input.
*/
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const std::string& element){

    // Copy until +1 size term, the last one is the end character (old version)
	// std::vector<char> inputss(element.c_str(), element.c_str()+element.size()+1);

	std::vector<char> inputss(element.c_str(), element.c_str()+element.size());
    buffer << (std::size_t)inputss.size();
    for (char & pp: inputss){
        buffer<<pp;
    }

    return buffer;

}

/*!
    Input stream operator for std::string
    \param[in] buffer is the input stream
    \param[in] element is the element to be streamed
    \result Returns the same input stream received in input.
*/
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, std::string& element){

	std::size_t nids;
    buffer >> nids;
    std::vector<char> inputss(nids);
    for (std::size_t i = 0; i < nids; ++i){
        buffer >> inputss[i];
    }

    // Copy until -1 term, the last one is the end character (old version)
    //  element = std::string(inputss.begin(), inputss.end()-1);

    element = std::string(inputss.begin(), inputss.end());


    return buffer;

}
