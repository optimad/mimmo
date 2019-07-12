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
#include <iterator>
#include <utility>
#include "BaseManipulation.hpp"

using namespace std;

namespace mimmo {

int BaseManipulation::sm_baseManipulationCounter(1);

/*!
 * Default constructor of BaseManipulation.
 * It sets to zero/null each member/pointer.
 */
BaseManipulation::BaseManipulation(){
    m_geometry      = NULL;
    m_portsType     = ConnectionType::BOTH;
    m_name          = "mimmo";
    m_active        = true;
    m_arePortsBuilt = false;
    m_execPlot      = false;
    m_outputPlot    = "./";
    m_counter       = sm_baseManipulationCounter;
    m_priority      = 0;
    m_apply         = false;
    sm_baseManipulationCounter++;

#if MIMMO_ENABLE_MPI
	int initialized;
	MPI_Initialized(&initialized);
	if (!initialized)
	   MPI_Init(NULL, NULL);

    //Fixed MPI comm world
	MPI_Comm_dup(MPI_COMM_WORLD, &m_communicator);

	MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &m_nprocs);

#endif

    bool logexists  = bitpit::log::manager().exists(MIMMO_LOG_FILE);
    if (m_counter == 1)
        initializeLogger(logexists);
    else
    	m_log = &bitpit::log::cout(MIMMO_LOG_FILE);

};

/*!
 * Default destructor of BaseManipulation.
 */
BaseManipulation::~BaseManipulation(){
    clear();
    deletePorts();
};

/*!
 * Copy constructor of BaseManipulation. Preexistent port connections will not be copied.
 * \param[in] other class of type BaseManipulation
 */
BaseManipulation::BaseManipulation(const BaseManipulation & other){
    m_geometry      = other.m_geometry;
    m_portsType     = other.m_portsType;
    m_name          = other.m_name;
    m_active        = other.m_active;
    m_arePortsBuilt = false;
    m_execPlot      = other.m_execPlot;
    m_outputPlot    = other.m_outputPlot;

    //only for copy construction
    m_counter       = sm_baseManipulationCounter;
    sm_baseManipulationCounter++;

    m_priority      = other.m_priority;
    m_apply         = other.m_apply;

    //logger is ready, since another BaseManipulation other, is instantiated.
    m_log           = &bitpit::log::cout(MIMMO_LOG_FILE);

#if MIMMO_ENABLE_MPI
	m_rank = other.m_rank;
	m_nprocs = other.m_nprocs;
	MPI_Comm_dup(other.m_communicator, &m_communicator);
#endif
};

/*!
 * Assignement operator of BaseManipulation.
 * \param[in] other class of type BaseManipulation
 */
BaseManipulation & BaseManipulation::operator=(const BaseManipulation & other){

    //NOT CAS type, since the class is abstract.
    //Please note static members or link to static members are not copied.
    //Port connection and port building members are not copied.
    m_geometry      = other.m_geometry;
    m_portsType     = other.m_portsType;
    m_name          = other.m_name;
    m_active        = other.m_active;
    m_execPlot      = other.m_execPlot;
    m_outputPlot    = other.m_outputPlot;
    m_priority      = other.m_priority;
    m_apply         = other.m_apply;
#if MIMMO_ENABLE_MPI
	MPI_Comm_dup(other.m_communicator, &m_communicator);
	m_rank			= other.m_rank;
	m_nprocs		= other.m_nprocs;
#endif
    return  *this;
};

/*!
 * Exchange data between the current class and an external BaseManipulation object.
 * Data will be swapped (data of A assigned to B and vice versa). Port connections/building info
 * will not be exchanged, as well as class unique-id.
 * \param[in] x BaseManipulation object to be swapped with current class.
 */
void BaseManipulation::swap(BaseManipulation & x) noexcept
{
    std::swap(m_priority, x.m_priority);
    std::swap(m_name, x.m_name);
    std::swap(m_geometry, x.m_geometry);
    std::swap(m_portsType, x.m_portsType);
    std::swap(m_active, x.m_active);
    std::swap(m_execPlot, x.m_execPlot);
    std::swap(m_apply, x.m_apply);
    std::swap(m_outputPlot, x.m_outputPlot);
#if MIMMO_ENABLE_MPI
    std::swap(m_communicator, x.m_communicator);
    std::swap(m_rank, x.m_rank);
    std::swap(m_nprocs, x.m_nprocs);
#endif

}

/*! Initialize the logger.
 *  \param[in] logexist does a logger already exist?
 * NOTE: console verbosity set to NORMAL as default (file verbosity set to DEBUG). If a logger already exists
 * the verbosities of the external logger are used.
 */
void
BaseManipulation::initializeLogger(bool logexist){
	if (!logexist){
		bitpit::log::manager().setMode(bitpit::log::Mode::COMBINED);
#if MIMMO_ENABLE_MPI
		bitpit::log::manager().create(MIMMO_LOG_FILE, false, m_nprocs, m_rank);
#else
		bitpit::log::manager().create(MIMMO_LOG_FILE, false);
#endif
	}
	m_log = &bitpit::log::cout(MIMMO_LOG_FILE);
	bitpit::log::setVisibility((*m_log), bitpit::log::Visibility::MASTER);
	bitpit::log::setConsoleVerbosity((*m_log), bitpit::log::NORMAL);
	bitpit::log::setFileVerbosity((*m_log), bitpit::log::NORMAL);
	if (m_counter == 1){
#if MIMMO_ENABLE_MPI
		if (m_rank == 0){
#endif
		(*m_log) << bitpit::log::priority(bitpit::log::NORMAL);
		(*m_log) << bitpit::log::context("mimmo");
		(*m_log) << "--------------------------------------------------" << std::endl;
		(*m_log) << "- mimmo - Surface manipulation and mesh morphing -" << std::endl;
		(*m_log) << "--------------------------------------------------" << std::endl;
		(*m_log) << "    " << std::endl;
		(*m_log) << bitpit::log::priority(bitpit::log::DEBUG);
#if MIMMO_ENABLE_MPI
		}
#endif
	}
}

/*!Get the logger.
 * \return Reference to logger object.
 */
bitpit::Logger&
BaseManipulation::getLog(){
    return (*m_log);
}

/*!
 * It gets if the ports of this object are already built.
 * \return True/false if ports are set.
 */
bool
BaseManipulation::arePortsBuilt(){
    return(m_arePortsBuilt);
}

/*!
 * Get current unsigned int marking priority of execution of class in a multi-chain
 * frame.
 * \return priority mark
 */
uint
BaseManipulation::getPriority(){
    return m_priority;
}

/*!
 * It gets the name of the manipulator object.
 * \return Name of the manipulator object.
 */
string
BaseManipulation::getName(){
    return m_name;
};

/*!
 * It gets the geometry linked by the manipulator object.
 * \return Pointer to geometry to be deformed by the manipulator object.
 */
MimmoObject*
BaseManipulation::getGeometry(){
    return m_geometry;
};

/*!
 * It gets the number of parents linked to the manipulator object.
 * \return Number of parents.
 */
int
BaseManipulation::getNParent(){
    return m_parent.size();
};

/*!
 * It gets the manipulator object linked by this object.
 * \param[in] i Index of target parent.
 * \return Pointer to i-th parent manipulator object.
 */
BaseManipulation*
BaseManipulation::getParent(int i){
    if (i>(int)m_parent.size()-1) return NULL;
    return next(m_parent.begin(), i)->first;
};

/*!
 * Return true if the target is contained in the parent list.
 * \param[in] target BaseManipulation target object
 * \param[out] index Actual position in the list (may vary if the parent list is modified).
 * \return false if not found.
 */
bool
BaseManipulation::isParent(BaseManipulation * target, int &index){
    unordered_map<BaseManipulation *, int>::iterator it;
    it = m_parent.find(target);
    index = -1;
    if(it == m_parent.end()) return false;

    index = distance(m_parent.begin(), it);
    return true;
};

/*!
 * It gets the number of children linked to the manipulator object.
 * \return Number of children.
 */
int
BaseManipulation::getNChild(){
    return m_child.size();
};

/*!
 * It gets one child object linked by this object.
 * \param[in] i Index of target child.
 * \return Pointer to i-th child manipulator object.
 */
BaseManipulation*
BaseManipulation::getChild(int i){
    if (i>(int)m_child.size()-1) return NULL;
    return next(m_child.begin(), i)->first;
};

/*!
 * Return true if the target is contained in the child list.
 * \param[in] target BaseManipulation target object
 * \param[out] index Actual position in the list.
 * \return False if not found.
 */
bool
BaseManipulation::isChild(BaseManipulation * target, int &index){
    unordered_map<BaseManipulation *, int>::iterator it;
    it = m_child.find(target);
    index = -1;
    if(it == m_child.end()) return false;

    index = distance(m_child.begin(), it);
    return true;
};

/*!
 * It gets the ports-type of the manipulation object.
 * \return ConnectionType of manipulation object (bi-directional, only backward, only forward).
 */
BaseManipulation::ConnectionType
BaseManipulation::getConnectionType(){
    return (m_portsType);
}

/*!
 * It gets the number of input ports of the object.
 * \return Number of input ports of the object.
 */
int
BaseManipulation::getNPortsIn(){
    return (m_portIn.size());
}

/*!
 * It gets the number of output ports of the object.
 * \return Number of output ports of the object.
 */
int
BaseManipulation::getNPortsOut(){
    return (m_portOut.size());
}

/*!
 * It gets all the input ports of the object
 * \return Map with PortID as key and pointer to input ports as value.
 */
std::unordered_map<PortID, PortIn*>
BaseManipulation::getPortsIn(){
    return (m_portIn);
}

/*!
 * It gets all the output ports of the object
 * \return Map with PortID as key and pointer to output ports as value.
 */
std::unordered_map<PortID, PortOut*>
BaseManipulation::getPortsOut(){
    return (m_portOut);
}

/*!
 * \return true if the feature to plot optional results of the block is active.
 */
bool
BaseManipulation::isPlotInExecution(){
    return (m_execPlot);
}

/*!
 * \return true if the feature to apply the results of the block is active.
 */
bool
BaseManipulation::isApply(){
    return (m_apply);
}

/*!
 * It gets if the object is activates or disable during the execution.
 * \return True/false if the object is activates or disable during the execution.
 */
bool
BaseManipulation::isActive(){
    return (m_active);
}

/*!
 * \return integer identifier of the object (deprecated : use getId())
 */
int
BaseManipulation::getClassCounter(){
    return getId();
}

/*!
 * \return integer identifier of the object
 */
int
BaseManipulation::getId(){
    return m_counter;
}

/*!Set the logger.
 * \param[in] Pointer to logger object.
 */
void
BaseManipulation::setLog(bitpit::Logger& log){
	m_log = &log;
}

/*!
 * Set unsigned int priority of execution of class in a multi-chain
 * frame. value must be > 0;
 * \param[in] priority priority mark
 */
void
BaseManipulation::setPriority(uint priority){
    m_priority = priority;
};

/*!
 * It sets the name of the manipulator object.
 * \param[in] name Name of the manipulator object.
 */
void
BaseManipulation::setName(std::string name){
    m_name = name;
};

/*!It sets the geometry linked by the manipulator object.
 * \param[in] geometry Pointer to geometry to be deformed by the manipulator object.
 */
void
BaseManipulation::setGeometry(MimmoObject* geometry){
    m_geometry = geometry;
};

/*!
 * Activates the feature to plot optional results of the block.
 * This features is meant for debug purpose only.
 * \param[in] flag true/false to activate/deactivate the feature
 */
void
BaseManipulation::setPlotInExecution( bool flag){
    m_execPlot = flag;
}

/*!
 * Set path to directory where the optional results will be stored,
 * if setPlotInExecution feature is set active.
 * \param[in] path absolute path to specified directory, if empty use the default value "."
 */
void
BaseManipulation::setOutputPlot(std::string path){
    if(path.empty())	path = ".";
    m_outputPlot = path;
}

/*!
 * Set (force) integer identifier of the object (deprecated : use setId())
 * \param[in] id integer identifier
 */
void
BaseManipulation::setClassCounter(int id){
    setId(id);
}

/*!
 * Activates the feature to apply directly the results of the block if possible.
 * \param[in] flag true/false to activate/deactivate the feature
 */
void
BaseManipulation::setApply( bool flag){
    m_apply = flag;
}

/*!
 * Set (force) integer identifier of the object
 * \param[in] id integer identifier
 */
void
BaseManipulation::setId(int id){
    m_counter = id;
}

/*!
 * It activates the object during the execution.
 */
void
BaseManipulation::activate(){
    m_active = true;
};

/*!
 * It disables the object during the execution.
 */
void
BaseManipulation::disable(){
    m_active = false;
};

/*!
 * It clears the pointer to the geometry linked by the object.
 */
void
BaseManipulation::unsetGeometry(){
    m_geometry = NULL;
};

/*!
 * It removes all the pins (connections) of the object and the
 * related pins of the linked objects.
 */
void
BaseManipulation::removePins(){
    removePinsIn();
    removePinsOut();
}

/*!
 * It removes all the input pins (connections) of the object and the related
 * output pins (connections) of the linked objects.
 */
void
BaseManipulation::removePinsIn(){
    unordered_map<BaseManipulation*, int>::iterator it;
    while(m_parent.size()){
        it = m_parent.begin();
        mimmo::pin::removeAllPins(it->first, this);
    }
}

/*!
 * It removes all the output (connections) pins of the object and the related
 * input pins (connections) of the linked objects.
 */
void
BaseManipulation::removePinsOut(){
    unordered_map<BaseManipulation*, int>::iterator it;
    while(m_child.size()){
        it = m_child.begin();
        mimmo::pin::removeAllPins(this, it->first);
    }
}

/*!
 * It clears the object, by setting to zero/NULL each member/pointer in the object.
 */
void
BaseManipulation::clear(){
    unsetGeometry();
    removePins();
    m_execPlot = false;
    m_outputPlot = ".";
};

/*!
 * Execution command. exec() runs the execution of output pins (connections) at the end of the execution.
 * execute is pure virtual and it has to be implemented in a derived class.
 */
void
BaseManipulation::exec(){

    if (!MIMMO_EXPERT){
        std::map<int, vector<PortIn*> > families;
        std::map<int, vector<PortID> > familiesID;
        for (std::unordered_map<PortID, PortIn*>::iterator i=m_portIn.begin(); i!=m_portIn.end(); i++){
            if (i->second->isMandatory()){
                families[i->second->getFamily()].push_back(i->second);
                familiesID[i->second->getFamily()].push_back(i->first);
            }
        }
        std::map<int, vector<PortID> >::iterator itID = familiesID.begin();
        for (std::map<int, vector<PortIn*> >::iterator i=families.begin(); i!=families.end(); i++){
            if (i->first == 0){
                int count = 0;
                for (auto ip : i->second){
                    std::vector<BaseManipulation*>  linked = ip->getLink();
                    if (linked.size() == 0){
                        (*m_log) << "error: " << m_name << " mandatory port " << itID->second[count] << " not linked -> exit! " << std::endl;
                        throw std::runtime_error (m_name + " mandatory port " + itID->second[count] + " not linked");
                    }
                    count++;
                }
            }
            else{
                bool check = false;
                for (auto ip : i->second){
                    std::vector<BaseManipulation*>  linked = ip->getLink();
                    check = check || linked.size();
                }
                if (!check){
                    (*m_log) << "error: " << m_name << " none of mandatory ports " << itID->second << " linked -> exit! " << std::endl;
                    std::string ports;
                    for (auto p : itID->second)
                        ports += " " + p + " ";
                    throw std::runtime_error (m_name + " none of mandatory ports " + ports + " linked");
                }
            }
            itID++;
        }
    }

    if (m_active) execute();

    for (std::unordered_map<PortID, PortOut*>::iterator i=m_portOut.begin(); i!=m_portOut.end(); i++){
        std::vector<BaseManipulation*>	linked = i->second->getLink();
        if (linked.size() > 0){
            i->second->exec();
        }
    }

    if(isPlotInExecution())	plotOptionalResults();
    if(isApply()) apply();
}

/*!
 * Base method to absorb parameter infos from an XML parser class of bitpit.This method
 * absorbs only the base attributes specified in the BaseManipulation class general documentation.
 * Derived classes using their own attributes, need to write a personal absorbSectionXML.
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name      string name associated to the slot. OPTIONAL
 */
void
BaseManipulation::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){
    BITPIT_UNUSED(name);

    std::string input;
    if(slotXML.hasOption("Priority")){
        input = slotXML.get("Priority");
        int value =0;
        if(!input.empty()){
            std::stringstream ss(bitpit::utils::string::trim(input));
            ss>>value;
        }
        setPriority(value);
    };


    if(slotXML.hasOption("Apply")){
        std::string input = slotXML.get("Apply");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setApply(value);
    }

    if(slotXML.hasOption("PlotInExecution")){
        std::string input = slotXML.get("PlotInExecution");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setPlotInExecution(value);
    }

    if(slotXML.hasOption("OutputPlot")){
        std::string input = slotXML.get("OutputPlot");
        input = bitpit::utils::string::trim(input);
        std::string temp = ".";
        if(!input.empty())  setOutputPlot(input);
        else                setOutputPlot(temp);
    }

}

/*!
 * Base method to flush parameter infos to an XML parser class of bitpit. This method
 * flushes only the base attributes specified in the BaseManipulation class general documentation.
 * Derived classes using their own attributes, need to write a personal flushSectionXML.
 * \param[in] slotXML	reference to a Section slot of bitpit::Config class.
 * \param[in] name      string name associated to the slot.OPTIONAL
 */
void
BaseManipulation::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);

    slotXML.set("ClassName", m_name);
    slotXML.set("Priority", std::to_string(getPriority()));

    if(isApply()){
        slotXML.set("Apply", std::to_string(1));
    }
    if(isPlotInExecution()){
        slotXML.set("PlotInExecution", std::to_string(1));
        slotXML.set("OutputPlot", m_outputPlot);
    }
}

/*!
 * Protected utility to delete all port of a class BaseManipulation
 */
void
BaseManipulation::deletePorts(){
    for (std::unordered_map<PortID, PortOut*>::iterator i = m_portOut.begin(); i != m_portOut.end(); i++){
        delete i->second;
        i->second = NULL;
    }
    for (std::unordered_map<PortID, PortIn*>::iterator i = m_portIn.begin(); i != m_portIn.end(); i++){
        delete i->second;
        i->second = NULL;
    }
}

/*!
 * It sets the buffer stored in an input port of the object.
 * \param[in] port ID of the input port.
 * \param[in] input Reference to bitpit::IBinaryStream to store in m_ibuffer member of the port.
 */
void
BaseManipulation::setBufferIn(PortID port, bitpit::IBinaryStream& input){
    m_portIn[port]->m_ibuffer = input;
}

/*!
 * It reads the buffer stored in an input port of the object.
 * \param[in] port ID of the port that reads the buffer and stores the
 * value in the related variable.
 */
void
BaseManipulation::readBufferIn(PortID port){
    m_portIn[port]->readBuffer();
}

/*!
 * It cleans the buffer stored in an input port of the object.
 * \param[in] port ID of the port.
 */
void
BaseManipulation::cleanBufferIn(PortID port){
    m_portIn[port]->cleanBuffer();
}

/*!
 * It adds a manipulator object linked by this object.
 * \param[in] parent Pointer to parent manipulator object.
 */
void
BaseManipulation::addParent(BaseManipulation* parent){

    if(!m_parent.count(parent)){
        m_parent.insert(pair<BaseManipulation*,int>(parent,1)); //add new parent with counter 1;
    }else{
        m_parent[parent]++; //just incrementing pre-existent parent counter;
    }
};

/*!
 * It adds a child manipulator object to the children linked by this object.
 * \param[in] child Pointer to child manipulator object.
 */
void
BaseManipulation::addChild(BaseManipulation* child){
    if(!m_child.count(child)){
        m_child.insert(pair<BaseManipulation*,int>(child,1)); //add new child with counter 1;
    }else{
        m_child[child]++; //just incrementing pre-existent child counter;
    }
};

/*!
 * Decrement target parent multiplicity, contained in member m_parent.
 * If multiplicity is zero, erase target from list.
 * The method is meant to be used together to manual cut off of object pins.
 * \param[in] parent Pointer to BaseManipulation object
 */
void
BaseManipulation::unsetParent(BaseManipulation * parent){

    unordered_map<BaseManipulation*, int>::iterator got = m_parent.find(parent);
    if(got != m_parent.end()){
        m_parent[parent]--;
        if(m_parent[parent] <1)	m_parent.erase(parent);
    }
};

/*!
 * Decrement target child multiplicity, contained in member m_child.
 * If multiplicity is zero, erase target from list.
 * The method is meant to be used in together to manual cut off of object pins.
 * \param[in] child Pointer to BaseManipulation object
 */
void
BaseManipulation::unsetChild(BaseManipulation * child){
    unordered_map<BaseManipulation*, int>::iterator got = m_child.find(child);
    if(got != m_child.end()){
        m_child[child]--;
        if(m_child[child] <1) m_child.erase(child);
    }
};

/*!
 * It finds an input pin (connection) of the object
 * \param[in] pin Target pin (connection).
 * \return Index of target pin in the input pins structure.
 * Return "NONE" if pin (connection) not found.
 */
PortID
BaseManipulation::findPinIn(PortIn& pin){
    for (std::unordered_map<PortID, PortIn*>::iterator i=m_portIn.begin(); i!=m_portIn.end(); i++){
        if (pin == *(i->second)) return(i->first);
    }
    return("NONE");
}

/*!
 * It finds an output pin (connection)  of the object
 * \param[in] pin Target pin (connection).
 * \return Index of target pin in the output pins structure.
 * Return "NONE" if pin (connection) not found.
 */
PortID
BaseManipulation::findPinOut(PortOut& pin){
    for (std::unordered_map<PortID, PortOut*>::iterator i=m_portOut.begin(); i!=m_portOut.end(); i++){
        if (pin == *(i->second)) return(i->first);
    }
    return("NONE");
}

/*!It adds an input pin (connection) of the object.
 * \param[in] objIn Pointer to sender BaseManipulation object.
 * \param[in] portR ID of target input port.
 */
void
BaseManipulation::addPinIn(BaseManipulation* objIn, PortID portR){
    if (objIn != NULL && m_portIn.count(portR) !=0 ){
        m_portIn[portR]->m_objLink.push_back(objIn);
    }
};

/*!It adds an output pin (connection) of the object.
 * \param[in] objOut Pointer to receiver BaseManipulation object.
 * \param[in] portS ID of target output port of sender.
 * \param[in] portR ID of target input port of receiver.
 */
void
BaseManipulation::addPinOut(BaseManipulation* objOut, PortID portS, PortID portR){
    if (objOut != NULL && m_portOut.count(portS) != 0 ){
        if (objOut->m_portIn.count(portR) != 0){
            m_portOut[portS]->m_objLink.push_back(objOut);
            m_portOut[portS]->m_portLink.push_back(portR);
        }
    }
};

/*!It removes an input pin (connection) of the object and the related output pin (connection) of the linked object.
 * \param[in] objIn Pointer to sender BaseManipulation object.
 * \param[in] portR ID of target input port of receiver.
 */
void
BaseManipulation::removePinIn(BaseManipulation* objIn, PortID portR){
    if (objIn != NULL && m_portIn.count(portR) != 0){
        std::vector<BaseManipulation*>	linked = m_portIn[portR]->getLink();
        for (int i=0; i<(int)linked.size(); i++){
            if (linked[i] == objIn){
                m_portIn[portR]->clear(i);
            }
        }
    }
};

/*!It removes an output pin (connection) of the object and the related input pin (connection) of the linked object.
 * \param[in] objOut Pointer to receiver BaseManipulation object.
 * \param[in] portS ID of target output port of sender.
 */
void
BaseManipulation::removePinOut(BaseManipulation* objOut, PortID portS){
    if (objOut != NULL && m_portOut.count(portS) != 0){
        std::vector<BaseManipulation*>	linked = m_portOut[portS]->getLink();
        for (int i=0; i<(int)linked.size(); i++){
            if (linked[i] == objOut){
                m_portOut[portS]->clear(i);
            }
        }
    };
}

/*!
 * It removes an input pin (connection) of the object and the
 * related output pin (connection) of the linked object.
 * \param[in] portR Port ID of target input pin (connection).
 * \param[in] j Index of target output pin (connection) of the portR input port.
 */
void
BaseManipulation::removePinIn(PortID portR, int j){
    if ( m_portIn.count(portR) != 0 ){
        if (j<(int)m_portIn[portR]->getLink().size() && j >= 0){
            m_portIn[portR]->clear(j);
        }
    }
}

/*!It removes an output pin (connection) of the object and the
 * related input pin (connection) of the linked object.
 * \param[in] portS Port ID of target output pin (connection).
 * \param[in] j Index of target output pin (connection) of the portS output port.
 */
void
BaseManipulation::removePinOut(PortID portS, int j){
    if ( m_portOut.count(portS) != 0 ){
        if (j<(int)m_portOut[portS]->getLink().size() && j >= 0){
            m_portOut[portS]->clear(j);
        }
    }
}

/*!
 * Plot optional data related to your blocks. The current method for the BaseManipulation class
 * do nothing. For Programmers only: to customize your optional results plotting create your own reimplementation
 * of this method in each BaseManipulation derived block.
 */
void 	BaseManipulation::plotOptionalResults(){};

/*!
 * Apply directly the result of your blocks without the need of a follower applier block. The current method for the BaseManipulation class
 * do nothing. For Programmers only: to customize your optional direct apply create your own implementation
 * of this method in each BaseManipulation derived block.
 */
void    BaseManipulation::apply(){};


/*!
 * Method to return all possible executable bloks of type BaseManipulation
 * embedded in the current class. For example create a BaseManipulation executable object with its own ports
 * which wraps a certain number of other BaseManipulation executable objects with their own i/o ports,
 * independent from their father.
 * By default this method return nothing. It needs to be customized according to functionality of your executable
 * derived class.
 * \return list of children BaseManipulation object pointers
 */
std::vector<BaseManipulation*>    BaseManipulation::getSubBlocksEmbedded(){
    return std::vector<BaseManipulation *>();
};

};
