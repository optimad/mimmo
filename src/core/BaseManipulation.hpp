/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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
#ifndef __BASEMANIPULATION_HPP__
#define __BASEMANIPULATION_HPP__

#include "factory.hpp"
#include "MimmoNamespace.hpp"
#include "MimmoObject.hpp"
#include "InOut.hpp"


#include <string>
#include <functional>
#include <unordered_map>
#include <typeinfo>

namespace mimmo{

using namespace pin;

class IOData;

/*!
 * \class BaseManipulation 
 * \brief BaseManipulation is the base class of any object (derived class) for manipulation of the geometry.
 *
 * BaseManipulation is the base class used to build each generic or particular manipulation object.
 * This base class has some common interface methods, as the base get/set methods.
 * The only method to be called to execute the manipulation object is the method exec() that calls the pure virtual execute().
 * Each manipulation base object has a linked geometry (a target mimmo object) and one or more linked manipulation
 * objects from wich recovering/distributing relevant data (as number of degrees of freedom, initial displacements or other). 
 * The exchange of such data is realized through Ports (input/output link to same class objects). See mimmo::pin namespace
 * for further information about the linking procedure of BaseManipulation objects.
 * Please note, when you get a copy by deputed operators/constructors of the class, pins, parental linking as well as custom 
 * input/result slots are not copied and left empty. Copy of BaseManipulation members in itself or in its derivations retains 
 * only the following parameters of BaseManipulation: link to target geometry, its own name and its supported Type of Ports.
 *
 *
 * PortType specification : please refer to PortType enum documentation.
 *
 */
class BaseManipulation{

	friend bool pin::addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced);
	friend void pin::removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);
	friend void pin::removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);
	friend bool pin::checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);
	friend PortOut;

public:
	//type definitions
	typedef std::unordered_map<BaseManipulation*, int>	bmumap;			/**<Unordered map type used for parent/child storing.*/
	typedef pin::ConnectionType					        ConnectionType;	/**<Connection type specification for Manipulation object.*/
	typedef	short int									PortID;			/**<Port ID (position of slot).*/
	
protected:
	uint                        m_priority;         /**<Flag marking priority of execution of class for multichain exec >*/
	std::string					m_name;				/**<Name of the manipulation object.*/
	int 						m_counter;			/**<Counter associated to the object */
	MimmoObject*				m_geometry;			/**<Pointer to manipulated geometry. */
	bmumap						m_parent;			/**<Pointers list to manipulation objects FATHER of the current class. List retains for each
														pointer a counter. When this counter is 0, pointer is released*/
	bmumap						m_child;			/**<Pointers list to manipulation objects CHILD of the current class.List retains for each
														pointer a counter. When this counter is 0, pointer is released*/

	ConnectionType				m_portsType;		/**<Type of ports of the object: BOTH (bidirectional),
														BACKWARD (only input) or FORWARD (only output).*/
	std::map<PortID, PortIn*>	m_portIn;			/**<Input ports map. */
	std::map<PortID, PortOut*>	m_portOut;			/**<Output ports map. */
	bool						m_arePortsBuilt;	/**<True or false is the ports are already set or not.*/

	bool						m_active;			/**<True/false to activate/disable the object.*/
	bool						m_execPlot; 		/**<Activate plotting of optional result directly in execution*/
	std::string					m_outputPlot;		/**<Define path for plotting optional results in execution*/


	//static members
    static  int                 sm_baseManipulationCounter;    /**<Current global number of BaseManipulation object in the instance. */

	
public:
	BaseManipulation();
	virtual ~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	bool				arePortsBuilt();

	uint 				getPriority();
	std::string			getName();
	MimmoObject*		getGeometry();
	int					getNParent();
	BaseManipulation*	getParent(int i = 0);
	bool				isParent(BaseManipulation *, int&);
	int					getNChild();
	BaseManipulation*	getChild(int i = 0);
	bool				isChild(BaseManipulation *, int&);
	ConnectionType		getConnectionType();
	int 				getNPortsIn();
	int 				getNPortsOut();
	std::map<PortID, PortIn*> getPortsIn();
	std::map<PortID, PortOut*>getPortsOut();
	
	bool	isPlotInExecution();
	bool	isActive();
    int     getClassCounter();
    int     getId();

	void 	setPriority(uint priority);
	void	setName(std::string name);
	void 	setGeometry(MimmoObject* geometry);
	void	setPlotInExecution(bool);
	void	setOutputPlot(std::string path);
    void    setClassCounter(int );
    void    setId(int );
	
	void	activate();
	void	disable();

	void 	unsetGeometry();
	void 	removePins();
	void 	removePinsIn();
	void 	removePinsOut();
	void	clear();
	
	void 	exec();
	
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
	

protected:

	virtual void 	buildPorts() = 0;
	void			deletePorts();
	
	template<typename T, typename O>
	bool	createPortOut(O* obj, T (O::*getVar_)(), PortID portS, containerTAG conType, dataTAG dataType);

	template<typename T, typename O>
	bool	createPortOut(T* var_, PortID portS, containerTAG conType, dataTAG dataType);

	template<typename T, typename O>
	bool	createPortIn(T* var_, PortID portR, containerTAG conType, dataTAG dataType);

	template<typename T, typename O>
	bool	createPortIn(O* obj, void (O::*setVar_)(T), PortID portR, containerTAG conType, dataTAG dataType);

	void	setBufferIn(PortID port, bitpit::IBinaryStream& input);
	void	readBufferIn(PortID port);
	void	cleanBufferIn(PortID port);

	void	addParent(BaseManipulation* parent);
	void	addChild(BaseManipulation* child);
	void 	unsetParent(BaseManipulation * parent);
	void 	unsetChild(BaseManipulation * child);
	
	PortID	findPinIn(PortIn& pin);
	PortID	findPinOut(PortOut& pin);

	void	addPinIn(BaseManipulation* objIn, PortID portR);
	void	addPinOut(BaseManipulation* objOut, PortID portS, PortID portR);

	void	removePinIn(BaseManipulation* objIn, PortID portR);
	void	removePinOut(BaseManipulation* objOut, PortID portS);

	void 	removePinIn(PortID portR, int j);
	void 	removePinOut(PortID portS, int j);

	/*!Execution method.
	 * Pure virtual method, it has to be implemented in derived classes.
	 */
	virtual void 	execute() = 0;
	
	virtual void 	plotOptionalResults();
};


//==============================//
//DATA CLASS
//==============================//

//==============================//
//DATA BASE CLASS
//==============================//
/*!
 *
 * \class IOData
 * \brief IOData is the base class of generic data stored as input or result in a manipulation object.
 *
 * IOData represents a generic data used in manipulation objects. Get and set methods of
 * this class are templated methods on the type of attached data.
 *
 */
class IOData{

public:
	/*!Default constructor of IOData.
	 */
	IOData(){};

	/*!Copy constructor of IOData.
	 */
	IOData(const IOData & other){
		*this = other;
	}

	/*!Assignement operator of IOData.
	 */
	IOData & operator=(const IOData & other){
		BITPIT_UNUSED(other);
		return (*this);
	}

	
	/*!It gets the data stored in the object.
	 * Even if not declared as pure virtual its behavior is
	 * like a pure virtual method, i.e. an analogous method
	 * is implemented in the template derived class IODataT.
	 * \return Pointer to data stored.
	 */
	template<typename T>
	T* getData();

	/*!It sets the data stored in the object.
	 * Even if not declared as pure virtual its behavior is
	 * like a pure virtual method, i.e. an analogous method
	 * is implemented in the template derived class IODataT.
	 * \param[in] data Data to be stored.
	 */
	template<typename T>
	void setData(T data);

};

//==============================//
//DATA DERIVED TEMPLATE CLASS
//==============================//
/*!
 * \class IODataT
 * \brief IODataT is the templated class of generic data derived from IOData base class.
 *
 * IODataT stores a generic data used in manipulation objects.
 */
template<typename T>
class IODataT: public IOData{
public:
	T 				m_data;	/**<Data contained in the object.*/

public:
	/*!Default constructor of IODataT.
	 */
	IODataT(){};

	/*!Custom constructor of IODataT.
	 * \param[in] data Data to be stored.
	 */
	IODataT(T data){
		m_data = data;
	};

	/*!Default destructor of IODataT.
	 */
	~IODataT();

	/*!Copy constructor of IODataT.
	 */
	IODataT(const IODataT & other){
		*this = other;
	}

	/*!Assignement operator of IODataT.
	 */
	IODataT & operator=(const IODataT & other){
		this->m_data 	= other.m_data;
		return (*this);
	}

	/*!It sets the data stored in the object.
	 * \param[in] data Data to be stored.
	 */
	void setData(T data){
		m_data = data;
	}

	/*!It gets the data stored in the object.
	 * \return Pointer to data stored.
	 */
	T* getData(){
		return(&m_data);
	}

};

};

#include "BaseManipulation.tpp"

#endif /* __BASEMANIPULATION_HPP__ */
