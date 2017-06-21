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
#ifndef __INOUT_HPP__
#define __INOUT_HPP__

#include "mimmoTypeDef.hpp"
#include "binary_stream.hpp"
#include "MimmoNamespace.hpp"
#include "BasicShapes.hpp"
#include "MimmoObject.hpp"
#include "MimmoPiercedVector.hpp"
#include "TrackingPointer.hpp"
#include <functional>

namespace mimmo{

class BaseManipulation;

typedef short int                   PortID; /**< port ID typedef */
typedef mimmo::pin::containerTAG    containerTAG; /**< containerTAG enum typedef*/
typedef mimmo::pin::dataTAG         dataTAG; /**< dataTAG enum typedef */

typedef mimmo::MimmoPiercedVector<double>  dmpvector1D;   /**< mimmo custom typedef*/
typedef mimmo::MimmoPiercedVector<darray3E>  dmpvecarr3E;   /**< mimmo custom typedef*/

/*!
* \class DataType
* \brief Class DataType defines the container and the type of data communicated by ports.
* \ingroup core
* 
* Class retains two members, m_conType and m_dataType to indentify container and data type,
* respectively(see containerTAG and dataTAG enums). 
*/
class DataType{
public:
    containerTAG	m_conType; /**< type of container*/
    dataTAG			m_dataType;/**< type of data*/

public:
    DataType();
    DataType(containerTAG conType, dataTAG dataType);
    virtual ~DataType();

    DataType(const DataType & other);
    bool operator==(const DataType & other);

};

/*!
* \class PortOut
* \brief PortOut is the abstract PIN base class dedicated to exchange data from a target class to other ones (output).
* \ingroup core
* 
* A PIN is an object member of BaseManipulation object.
* Through a PIN two base manipulation objects are linked together. One of these two
* objects is the parent object that gives an output to the other one that takes
* this value as input. 
* 
* The class store the following data:
* 
* - a buffer to communicate output data (m_obuffer)
* - a list of pointer to BaseManipulation receivers (m_objLink)
* - a list of Ports integer identifiers, marking the input ports of receivers, where the data will be sent (m_portLink)
* - information on the container and data type exchanged (m_datatype)
*  
* In general, a set of data (still not specified in this abstract class) of type m_datatype, written in a buffer stream m_obuffer, 
* will be sent to a list of BaseManipulation objects/receivers. Input ports of receivers are responsible to decode the 
* data buffer sent, and make available the data to their respective receveir. 
* 
* The execution of the output PortT will automatically
* exchange the buffer data, pass it to the input ports connected and makes them reading and decoding the data.
*/
class PortOut{
public:
    //members
    bitpit::OBinaryStream           m_obuffer;	/**<Output buffer to communicate data.*/
    std::vector<BaseManipulation*>  m_objLink;	/**<Outputs object to which communicate the data.*/
    std::vector<PortID>             m_portLink;	/**<ID of the input ports of the linked objects.*/
    DataType                        m_datatype;	/**<TAG of type of data communicated.*/

public:
    PortOut();
    virtual ~PortOut();

    PortOut(const PortOut & other);
    bool operator==(const PortOut & other);

    std::vector<BaseManipulation*>	getLink();
    std::vector<PortID>				getPortLink();
    DataType						getDataType();

    /*!
     * Pure virtual function to write a buffer.
     */
    virtual void	writeBuffer() = 0;
    void 			cleanBuffer();

    void clear();
    void clear(int j);

    void exec();

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
* \class PortOutT
* \brief PortOutT is the PIN class to exchange output data from an object to others.
* \ingroup core
* 
* PortOutT is the template derived class of PortOut specifying the set of data 
* that need to exchanged.
* 
* PortOut stores the following members:
* 
* - pointer to the sender object, owner of the PIN (m_obj)
* - pointer to the type of sender object variable containing the data that need to be exchanged (m_var)
* - pointer to a "get" method of the sender object that recovers the data that need to be exchanged (O::*m_getVar)
* 
* The last two members are alternative to each other, depending on the sender object interface design.
* Once the data is recovered from its sender object, it is redistributed as in PortOut base class (see PortOut doc).
* 
* "Get" methods have to be function objects of the standard library (functional include) created by the bind method.
* The data value must be returned by this methods as a copy or pointer. 
*
*/
template<typename T, typename O>
class PortOutT: public PortOut {

public:

    O*  m_obj_;             /**<Object owner of the port.*/
    T   *m_var_;            /**<Linked variable to communicate.*/
    T   (O::*m_getVar_)();  /**<Pointer to function that recovers the data to communicate (alternative to linked variable).*/

public:
    PortOutT();
    PortOutT(T *var_);
    PortOutT(T *var_, DataType datatype);
    PortOutT(O* obj_, T (O::*getVar_)());
    PortOutT(O* obj_, T (O::*getVar_)(), DataType datatype);
    virtual ~PortOutT();

    PortOutT(const PortOutT & other);
    bool operator==(const PortOutT & other);

    void writeBuffer();

};



/*!
* \class PortIn
* \brief PortIn is the abstract PIN base class dedicated to carry data to a target class from other ones (input).
* \ingroup core
* 
* A PIN is an object member of BaseManipulation object.
* Through a PIN two base manipulation objects are linked together. One of these two
* objects is the parent object that gives an output to the other one that takes
* this value as input. 
* 
* The class store the following data:
* 
* - a buffer to communicate input data (m_ibuffer)
* - a list of pointer to BaseManipulation senders (m_objLink)
* - information on the container and data type exchanged (m_datatype)
*  
* In general, a set of data of type m_datatype, is sent from one or more senders 
* and read as a buffer stream m_ibuffer. This class is responsible to decode the data and handle with the problem to manage multiple data 
* coming from multiple senders and makes it available (how to read the data, handle with its multiplicity and makes it available its still
* not specified in this abstract class).
*/
class PortIn{
public:
    //members
    bitpit::IBinaryStream               m_ibuffer;          /**<input buffer to recover data.*/
    std::vector<BaseManipulation*>      m_objLink;          /**<Input objects from which recover the data. */
    DataType                            m_datatype;         /**<TAG of type of data communicated.*/
    bool                                m_mandatory;        /**<Does the port have to be mandatorily linked?.*/
    int                                 m_familym;          /**<Tag of family of mandatory alternative ports.
                                                                 At least one of the ports of the same family has to be linked.
                                                                 Family TAG = 0 -> no family, this port has to be mandatory linked.*/

public:
    PortIn();
    virtual ~PortIn();

    PortIn(const PortIn & other);
    bool operator==(const PortIn & other);

    std::vector<mimmo::BaseManipulation*>   getLink();
    DataType                                getDataType();
    bool                                    isMandatory();
    int                                     getFamily();

    void clear();
    void clear(int j);

    /*!
     * Pure virtual function to read a buffer.
     */
    virtual void    readBuffer() = 0;
    void            cleanBuffer();

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
* \class PortInT
* \brief PortInT is the PIN class to get input data arriving to an object from other objects.
* \ingroup core
* 
* PortInT is the template derived class of PortIn specifying the set of data 
* that need to exchanged. 
* 
* PortIn stores the following members:
* 
* - pointer to the receiver object, owner of the PIN (m_obj)
* - pointer to the type of receiver object variable, that will contain the data that need to be exchanged (m_var)
* - pointer to a "set" method of the receiver object that will store the data that need to be exchanged (O::*m_setVar)
* 
* The last two members are alternative to each other, depending on the receiver object interface design.
* If receiver's m_var or set method are capable to handle with multiple inputs,all inputs passing through the ports will be stored
* in the receiver class, otherwise multiple data will be progressively overwritten each other, till the last, who will remain. 
* Activation of the class will be automatically triggered by sender port PortOutT associated to it. See more in PortOut,PortOutT, 
* PortIn documentation.
* 
* "Set" methods have to be function objects of the standard library (functional include) created by the bind method.
* The data value must be passed as argument of the function by copy or pointer. 
*/
template<typename T, typename O>
class PortInT: public PortIn {

public:

    O*      m_obj_;             /**<Object owner of the port.*/
    T       *m_var_;            /**<Linked variable to fill with communicated data.*/
    void    (O::*m_setVar_)(T); /**<Pointer to function that fills members with the communicated (alternative to linked variable).*/

public:
    PortInT();
    PortInT(T *var_);
    PortInT(T *var_, DataType datatype, bool mandatory = false, int family = 0);
    PortInT(O* obj_, void (O::*setVar_)(T), bool mandatory = false, int family = 0);
    PortInT(O* obj_, void (O::*setVar_)(T), DataType datatype, bool mandatory = false, int family = 0);
    virtual ~PortInT();

    PortInT(const PortInT & other);
    bool operator==(const PortInT & other);

    void readBuffer();

};

}

/*!
 * \ingroup binaryStream
 * \{
 */
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,dvector1D& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const dvector1D& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,livector1D& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const livector1D& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,livector2D& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const livector2D& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,shivector1D& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const shivector1D& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,dvecarr3E& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const dvecarr3E& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::array<mimmo::CoordType,3>& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::array<mimmo::CoordType,3>& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, mimmo::ShapeType& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const mimmo::ShapeType& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::pair<mimmo::MimmoObject*, dvecarr3E *>& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::pair<mimmo::MimmoObject*, dvecarr3E *>& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::pair<mimmo::MimmoObject*, dvector1D *>& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::pair<mimmo::MimmoObject*, dvector1D *>& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::vector<mimmo::TrackingPointer * >& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector<mimmo::TrackingPointer * >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::pair<mimmo::BaseManipulation*, double>& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::pair<mimmo::BaseManipulation *, double>& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,dvecarr2E& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const dvecarr2E& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,ivecarr2E& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const ivecarr2E& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,std::vector< dvecarr2E >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector< dvecarr2E >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::unordered_map< std::string, std::pair<int, mimmo::MimmoObject* > >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::unordered_map< std::string, std::pair<int, mimmo::MimmoObject* > >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::unordered_map< long, std::pair<int, long> >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::unordered_map< long, std::pair<int, long > >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::vector< mimmo::MimmoObject* >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector< mimmo::MimmoObject* >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, mimmo::FileDataInfo&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const mimmo::FileDataInfo& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::vector<mimmo::FileDataInfo>&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector<mimmo::FileDataInfo>& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::unordered_map< mimmo::MimmoObject*, dvector1D* >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::unordered_map< mimmo::MimmoObject*, dvector1D* >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::unordered_map< mimmo::MimmoObject*, dvecarr3E* >&  element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::unordered_map< mimmo::MimmoObject*, dvecarr3E* >& element);


bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::vector< std::pair<mimmo::MimmoObject*, dvector1D *> >& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector< std::pair<mimmo::MimmoObject*, dvector1D *> >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::vector< std::pair<mimmo::MimmoObject*, dvecarr3E *> >& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::vector< std::pair<mimmo::MimmoObject*, dvecarr3E *> >& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,mimmo::dmpvector1D& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const mimmo::dmpvector1D& element);

bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf,mimmo::dmpvecarr3E& element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const mimmo::dmpvecarr3E& element);

/*!
 *\}
 */

#include "InOut.tpp"



#endif /* __INOUT_HPP__ */
