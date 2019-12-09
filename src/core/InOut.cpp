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
#include "InOut.hpp"
#include "MimmoNamespace.hpp"
#include "BaseManipulation.hpp"

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
*	Output stream operator for mimmo::ShapeType enum
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream  &buffer, const mimmo::ShapeType &var)
{
    buffer << static_cast<int> (var);
    return buffer;
}


/*!
*	Input stream operator for mimmo::ShapeType enum
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, mimmo::ShapeType &var)
{
    int val;
    buffer >> val;
    var = static_cast<mimmo::ShapeType>	(val);
    return buffer;
}

/*!
*	Output stream operator for mimmo::CoordType enum
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream  &buffer, const mimmo::CoordType &var)
{
    buffer << static_cast<int> (var);
    return buffer;
}


/*!
*	Input stream operator for mimmo::CoordType enum
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, mimmo::CoordType &var)
{
    int val;
    buffer >> val;
    var = static_cast<mimmo::CoordType>	(val);
    return buffer;
}

/*!
*	Input stream operator for mimmo::FileDataInfo
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, mimmo::FileDataInfo&  var){

    buffer >> var.ftype;
    buffer >> var.fdir;
    buffer >> var.fname;

    return buffer;
};

/*!
*	Output stream operator for mimmo::FileDataInfo
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const mimmo::FileDataInfo& var){

    buffer<<var.ftype;
    buffer<<var.fdir;
    buffer<<var.fname;
    return buffer;


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

//==============================================================//
// DATA TYPE  CLASS	IMPLEMENTATION								//
// //==============================================================//
namespace mimmo{

/*!
 * Default constructor of DataType
 */
DataType::DataType(){
    m_conType = "MC_SCALAR";
    m_dataType = "MD_INT";
};

/*!
 * Custom constructor of DataType.
* \param[in] conType string of type of port container.
* \param[in] dataType string of type of port data.
*/
DataType::DataType(containerTAG conType, dataTAG dataType){
    m_conType 	= conType;
    m_dataType	= dataType;
    return;
};

/*!
 * Default destructor of DataType
 */
DataType::~DataType(){};

/*!
 * Copy constructor of DataType.
 */
DataType::DataType(const DataType & other){
    m_conType 	= other.m_conType;
    m_dataType	= other.m_dataType;
    return;
};

/*!
 * Compare operator of DataType.
 */
bool DataType::operator==(const DataType & other){
    bool check = true;
    check = check && (m_conType == other.m_conType);
    check = check && (m_dataType == other.m_dataType);
    return(check);
};


//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION                              //
//==============================================================//

/*!
 * Default constructor of PortOut
 */
PortOut::PortOut(){
    m_objLink.clear();
};

/*!
 * Default destructor of PortOut
 */
PortOut::~PortOut(){};

/*!
 * Copy constructor of PortOut.
 */
PortOut::PortOut(const PortOut & other){
    m_objLink 	= other.m_objLink;
    m_obuffer	= other.m_obuffer;
    m_portLink	= other.m_portLink;
    m_datatype	= other.m_datatype;
    return;
};

/*!
 * Compare operator of PortOut.
 */
bool PortOut::operator==(const PortOut & other){
    bool check = true;
    check = check && (m_portLink == other.m_portLink);
    check = check && (m_objLink == other.m_objLink);
    check = check && (m_datatype == other.m_datatype);
    return(check);
};

/*!
 * It gets the objects linked by this port.
* \return vector of pointer to linked objects.
*/
std::vector<mimmo::BaseManipulation*>
PortOut::getLink(){
    return(m_objLink);
}

/*!
 * It gets the input port Markers of the objects linked by this port.
* \return Vector of PortID.
*/
std::vector<PortID>
PortOut::getPortLink(){
    return(m_portLink);
}

/*!
 * It gets the TAG of container/data types communicated by this port.
* \return TAG of container/data types communicated.
*/
DataType
PortOut::getDataType(){
    return(m_datatype);
}

/*!
 * It empties the output buffer.
 */
void
mimmo::PortOut::cleanBuffer(){
    m_obuffer.seekg(0);
}

/*!
 * It clears the links to objects and the related ports.
 */
void
mimmo::PortOut::clear(){
    m_objLink.clear();
    m_portLink.clear();
}

/*!
 * It removes the link to an object and the related port Marker.
* \param[in] j Index of the linked object in the links vector of this port.
*/
void
mimmo::PortOut::clear(int j){
    if (j < (int)m_objLink.size() && j >= 0){
        m_objLink.erase(m_objLink.begin() + j);
        m_portLink.erase(m_portLink.begin() + j);
    }
}

/*!
 * Execution of the PIN.
 * All the pins are called in execution of the sending owner after its own execution.
 * Reading stage of pin linked receivers is automatically performed within this execution.
 */
void
mimmo::PortOut::exec(){
    if (m_objLink.size() > 0){
        writeBuffer();
        mimmo::IBinaryStream input(m_obuffer.data(), m_obuffer.getSize());
        cleanBuffer();
        for (int j=0; j<(int)m_objLink.size(); j++){
            if (m_objLink[j] != NULL){
                m_objLink[j]->setBufferIn(m_portLink[j], input);
                m_objLink[j]->readBufferIn(m_portLink[j]);
                m_objLink[j]->cleanBufferIn(m_portLink[j]);
            }
        }
    }
};

/*!
 * Default constructor of PortIn
 */
PortIn::PortIn(){
    m_mandatory =false;
    m_familym = 0;
};

/*!
 * Default destructor of PortIn
 */
PortIn::~PortIn(){};

/*!
 * Copy constructor of PortIn.
*/
PortIn::PortIn(const PortIn & other){
    m_objLink 	= other.m_objLink;
    m_ibuffer	= other.m_ibuffer;
    m_datatype  = other.m_datatype;
    m_mandatory = other.m_mandatory;
    m_familym   = other.m_familym;
    return;
};

/*!
 * Compare operator of PortIn.
 */
bool PortIn::operator==(const PortIn & other){
    bool check = true;
    check = check && (m_objLink == other.m_objLink);
    check = check && (m_datatype == other.m_datatype);
    check = check && (m_mandatory = other.m_mandatory);
    check = check && (m_familym   = other.m_familym);
    return(check);
};

/*!
 * It gets the linked object by this port.
 * \return Pointer to linked object.
 */
std::vector<mimmo::BaseManipulation*>
PortIn::getLink(){
    return(m_objLink);
}

/*!
 * It gets the TAG of container/data types communicated by this port.
 * \return TAG of container/data type communicated.
 */
DataType
PortIn::getDataType(){
    return(m_datatype);
}

/*!
 * It gets if this port has to be mandatorily linked.
 * \return mandatory?.
 */
bool
PortIn::isMandatory(){
    return(m_mandatory);
}

/*!
 * It gets the family of this port.
 * \return mandatory family.
 */
int
PortIn::getFamily(){
    return(m_familym);
}

/*!
 * It clears all linked objects into this port.
 */
void
mimmo::PortIn::clear(){
    m_objLink.clear();
}

/*!
 * It removes the link to an object and the related port ID.
 * \param[in] j Index of the linked object in the links vector of this port.
 */
void
mimmo::PortIn::clear(int j){
    if (j < (int)m_objLink.size() && j >= 0){
        m_objLink.erase(m_objLink.begin() + j);
    }
}


/*!
 * It releases the memory occupied by the input buffer.
 */
void
mimmo::PortIn::cleanBuffer(){
    m_ibuffer.seekg(0);
}

}
