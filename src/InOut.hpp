/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __INOUT_HPP__
#define __INOUT_HPP__

#include "MiMMO_TypeDef.hpp"
#include "binary_stream.hpp"
#include "MimmoNamespace.hpp"
#include "BasicShapes.hpp"
#include "MimmoObject.hpp"
#include "TrackingPointer.hpp"
#include <functional>

namespace mimmo{

class BaseManipulation;

typedef short int 					PortID;
typedef mimmo::pin::containerTAG	containerTAG;
typedef mimmo::pin::dataTAG			dataTAG;


/*!
 *  \class DataType
 *	\brief Class DataType defines the container and the type of data communicated by ports.
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
	DataType & operator=(const DataType & other);
	bool operator==(const DataType & other);

};

/*!
 *  \class PortOut
 *	\brief PortOut is the input-output PIN base class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 */
class PortOut{
public:
	//members
	bitpit::OBinaryStream				m_obuffer;	/**<Output buffer to communicate data.*/
	std::vector<BaseManipulation*>		m_objLink;	/**<Outputs object to which communicate the data.*/
	std::vector<PortID>					m_portLink;	/**<ID of the input ports of the linked objects.*/
	DataType							m_datatype;	/**<TAG of type of data communicated.*/

public:
	PortOut();
	virtual ~PortOut();

	PortOut(const PortOut & other);
	PortOut & operator=(const PortOut & other);
	bool operator==(const PortOut & other);

	std::vector<BaseManipulation*>	getLink();
	std::vector<PortID>				getPortLink();
	DataType						getDataType();

	virtual void	writeBuffer() = 0;
	void 			cleanBuffer();

	void clear();
	void clear(int j);

	/*!Execution method.
	 */
	void exec();

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
 *  \class PortOutT
 *	\brief PortOutT is the input-output templated PIN class derived from base PortOut pin class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 *	The parent object gives the output value by a function set by the user
 *	during the creation of the pin (getVal functions), while the child object
 *	use this value by a set function (setVal function).
 *	The functions are given to the pin by get/set function objects with any
 *	type of value (the same for the two functions of the linked object) by using
 *	the template formulation. The functions have to be function objects of the standard library
 *	(functional include) created by the bind method.
 *	The input/output value can be returned by copy, reference or pointer by the get function
 *	of the parent and it can be passed by copy or pointer to the set value of the child object.
 *	The output pins of an object are executed after its own execution.
 *
 */
template<typename T, typename O>
class PortOutT: public PortOut {

public:

	O*		m_obj_;				/**<Object owner of the port.*/
	T		*m_var_;			/**<Linked variable to communicate.*/
	T		(O::*m_getVar_)();	/**<Pointer to function that recovers the data to communicate (alternative to linked variable).*/

public:
	PortOutT();
	PortOutT(T *var_);
	PortOutT(T *var_, DataType datatype);
	PortOutT(O* obj_, T (O::*getVar_)());
	PortOutT(O* obj_, T (O::*getVar_)(), DataType datatype);
	virtual ~PortOutT();

	PortOutT(const PortOutT & other);
	PortOutT & operator=(const PortOutT & other);
	bool operator==(const PortOutT & other);

	void writeBuffer();

};



/*!
 *  \class PortIn
 *	\brief PortOut is the input-output PIN base class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 */
class PortIn{
public:
	//members
	bitpit::IBinaryStream				m_ibuffer;	/**<input buffer to recover data.*/
	std::vector<BaseManipulation*>		m_objLink;	/**<Input objects from which	recover the data. */
	DataType							m_datatype;	/**<TAG of type of data communicated.*/

public:
	PortIn();
	virtual ~PortIn();

	PortIn(const PortIn & other);
	PortIn & operator=(const PortIn & other);
	bool operator==(const PortIn & other);

	std::vector<mimmo::BaseManipulation*>	getLink();
	DataType						getDataType();

	void clear();
	void clear(int j);


	virtual void	readBuffer() = 0;
	void 			cleanBuffer();

};


//==============================================================//
// TEMPLATE DERIVED INOUT CLASS									//
//==============================================================//

/*!
 *  \class PortInT
 *	\brief PortInT is the input-output templated PIN class derived from base PortIn pin class.
 *
 *	A pin is an object member of BaseManipulation object.
 *	Through a pin two base manipulation objects are linked together. One of these two
 *	objects is the parent object that gives an output to the other one that takes
 *	this value as input. Therefore a pin can be OR input OR output pin for an object.
 *
 *	The parent object gives the output value by a function set by the user
 *	during the creation of the pin (getVal functions), while the child object
 *	use this value by a set function (setVal function).
 *	The functions are given to the pin by get/set function objects with any
 *	type of value (the same for the two functions of the linked object) by using
 *	the template formulation. The functions have to be function objects of the standard library
 *	(functional include) created by the bind method.
 *	The input/output value can be returned by copy, reference or pointer by the get function
 *	of the parent and it can be passed by copy or pointer to the set value of the child object.
 *	The output pins of an object are executed after its own execution.
 *
 */
template<typename T, typename O>
class PortInT: public PortIn {

public:

	O*		m_obj_;				/**<Object owner of the port.*/
	T		*m_var_;			/**<Linked variable to fill with communicated data.*/
	void	(O::*m_setVar_)(T);	/**<Pointer to function that fills members with the communicated (alternative to linked variable).*/

public:
	PortInT();
	PortInT(T *var_);
	PortInT(T *var_, DataType datatype);
	PortInT(O* obj_, void (O::*setVar_)(T));
	PortInT(O* obj_, void (O::*setVar_)(T), DataType datatype);
	virtual ~PortInT();

	PortInT(const PortInT & other);
	PortInT & operator=(const PortInT & other);
	bool operator==(const PortInT & other);

	void readBuffer();

};

}

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

/*
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buf, std::string *   element);
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buf, const std::string * element);*/

#include "InOut.tpp"



#endif /* __INOUT_HPP__ */
