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
#ifndef __BASEMANIPULATION_HPP__
#define __BASEMANIPULATION_HPP__

#include "MimmoNamespace.hpp"
#include "MimmoObject.hpp"
#include "InOut.hpp"
#include "MimmoPiercedVector.hpp"

#include <factory.hpp>
#include <portManager.hpp>

#include <bitpit_common.hpp>

#include <string>
#include <functional>
#include <unordered_map>
#include <typeinfo>
#include <type_traits>

#if MIMMO_ENABLE_MPI
    #include <mpi.h>
#endif

#if defined(_WIN32)
    #define uint unsigned int
#endif

namespace mimmo{

/*!
 * \class BaseManipulation
 * \ingroup core
 * \brief BaseManipulation is the base class of any manipulation object of the library.
 *
 * BaseManipulation is the base class used to derive each manipulation object of mimmo API. \n
 * This base class provides some common interface methods, as the base get/set methods. \n
 * To execute an object derived from BaseManipulation call the
 * method exec() of the base class.
 * In the exec() function the pure virtual execute() method is called.
 * Pure virtual execute() is the real working function of a manipulation object and has to be
 * implemented in each derived class. \n
 * A manipulation base object has a linked geometry, a MimmoObject, that is
 * the target geometry to manipulate. \n
 * A set of manipulation objects can be rearranged in an execution chain and linked
 * in order to communicate each other useful data and geometries.
 * The exchange of such data is realized through ports of input/output;
 * two objects (sender/receiver) can be linked by two ports (output/input)
 * that communicate the same type of data. \n
 * The ports of a manipulation object are built in buildPorts() method. This function is
 * a pure virtual method for the base class, i.e. buildPorts() has to be implemented for each
 * derived. In buildPorts() of the derived class function will be defined all the input/output allowed connections
 * of a manipulation object. \n
 * See mimmo::pin namespace and examples for further
 * information about the linking procedure of BaseManipulation derived objects. \n
 * Please note, when copying a derived object by operators/constructors of the base class,
 * ports and linking are not
 * copied and left empty. Copy of BaseManipulation members in itself
 * or in its derivations retains only the following parameters of
 * BaseManipulation:
 * - link to target geometry;
 * - name;
 * - supported type of Ports.
 *
 *
 * For further information about PortType specification
 * please refer to PortType enum documentation. \n
 *
 * BaseManipulation controls a initial set of xml attributes which can be read from a xml file interface or written to it,
 * through absorbSectionXML/flushSectionXML methods. Such parameters are:

 * - <B>ClassName</B>: specific name identifying the class as "mimmo.XXXXX"
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply result directly in execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class, for debugging purpose.
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * All BaseManipulation derived classes inherite these attributes.
 */
class BaseManipulation{

    /*!
     * see pin::addPin
     */
    friend bool pin::addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced);
    /*!
     * see pin::removePin
     */
    friend void pin::removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);
    /*!
     * see pin::removeAllPins
     */
    friend void pin::removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);
    /*!
     * see pin::checkCompatibility
     */
    friend bool pin::checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);
    /*!
     * see PortOut::exec
     */
    friend void PortOut::exec();
    /*!
     * see mimmo::setLogger
     */
    friend void mimmo::setLogger(std::string log);
    /*!
     * see mimmo::setLoggerDirectory
     */
    friend void mimmo::setLoggerDirectory(std::string dir);

public:
    //type definitions
    /*! \ingroup typedefs{ */
    typedef std::unordered_map<BaseManipulation*, int>    bmumap;            /**<Unordered map type used for parent/child storing.*/
    typedef pin::ConnectionType                           ConnectionType;    /**<Connection type specification for Manipulation object.*/
    typedef std::string                                   PortID;            /**<Port ID (identifier of the port).*/
    /*! } */
protected:
    uint                            m_priority;         /**<Flag marking priority of execution of the object (0 - highest priority) >*/
    std::string                     m_name;             /**<Name of the manipulation object.*/
    int                             m_counter;          /**<Counter ID associated to the object */
    MimmoSharedPointer<MimmoObject> m_geometry;         /**<Mimmo shared pointer to manipulated geometry. */
    bmumap                          m_parent;           /**<Pointers list to manipulation objects FATHER of the current class. List retains for each
                                                        pointer a counter. When this counter is 0, pointer is released*/
    bmumap                          m_child;            /**<Pointers list to manipulation objects CHILD of the current class.List retains for each
                                                        pointer a counter. When this counter is 0, pointer is released*/

    ConnectionType                  m_portsType;        /**<Type of ports of the object: BOTH (bidirectional),
                                                            BACKWARD (only input) or FORWARD (only output).*/
    std::unordered_map<PortID, PortIn*>   m_portIn;         /**<Input ports map. */
    std::unordered_map<PortID, PortOut*>  m_portOut;        /**<Output ports map. */
    bool                                  m_arePortsBuilt;  /**<True or false is the ports are already set or not.*/

    bool                        m_active;        /**<True/false to activate/disable the object during the execution.*/
    bool                        m_execPlot;      /**<Activate plotting of optional result directly in execution.*/
    bool                        m_apply;         /**<Activate apply result directly in execution.*/
    std::string                 m_outputPlot;    /**<Define path for plotting optional results in execution.*/

    bitpit::Logger*             m_log;           /**<Pointer to logger.*/

    //static members
    static  int                 sm_baseManipulationCounter;     /**<Current global number of BaseManipulation object in the instance. */

#if MIMMO_ENABLE_MPI
    int							m_nprocs;			/**<Total number of processors.*/
    int							m_rank;				/**<Current rank number.*/
	MPI_Comm 					m_communicator; 	/**<MPI communicator.*/
#endif


public:
    BaseManipulation();
    virtual ~BaseManipulation();

    BaseManipulation(const BaseManipulation & other);
    BaseManipulation & operator=(const BaseManipulation & other);

    bitpit::Logger& getLog();

    bool                 arePortsBuilt();

    uint                 getPriority();
    std::string          getName();
    MimmoSharedPointer<MimmoObject>         getGeometry();
    MimmoSharedPointer<MimmoObject>&        getGeometryReference();
    int                  getNParent();
    BaseManipulation*    getParent(int i = 0);
    bool                 isParent(BaseManipulation *, int&);
    int                  getNChild();
    BaseManipulation*    getChild(int i = 0);
    bool                 isChild(BaseManipulation *, int&);
    ConnectionType       getConnectionType();
    int                  getNPortsIn();
    int                  getNPortsOut();

#if MIMMO_ENABLE_MPI
    int                 getProcessorCount();
    int                 getRank();
    MPI_Comm &          getCommunicator();
#endif


    std::unordered_map<PortID, PortIn*> getPortsIn();
    std::unordered_map<PortID, PortOut*>getPortsOut();

    bool    isPlotInExecution();
    bool    isActive();
    bool    isApply();
    int     getId();

    void	setLog(bitpit::Logger& log);
    void    setPriority(uint priority);
    void    setName(std::string name);
    void    setGeometry(MimmoSharedPointer<MimmoObject> geometry);
    void    setPlotInExecution(bool);
    void    setOutputPlot(std::string path);
    void    setId(int );
    void    setApply(bool flag = true);

    void    activate();
    void    disable();

    void    unsetGeometry();
    void    removePins();
    void    removePinsIn();
    void    removePinsOut();
    void    clear();

    void    exec();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

    virtual std::vector<BaseManipulation*> getSubBlocksEmbedded();

protected:

    void swap(BaseManipulation & x) noexcept;
    void initializeLogger(bool logexists);

    /*!
     * Build ports of the class.
     * Pure virtual method.
     */
    virtual void     buildPorts() = 0;
    void             deletePorts();

    template<typename T, typename O>
    bool    createPortOut(O* obj, T (O::*getVar_)(), PortID portS);

    template<typename T, typename O>
    bool    createPortOut(T* var_, PortID portS);

    template<typename T, typename O>
    bool    createPortIn(T* var_, PortID portR, bool mandatory = false, int family = 0);

    template<typename T, typename O>
    bool    createPortIn(O* obj, void (O::*setVar_)(T), PortID portR, bool mandatory = false, int family = 0);

    void    setBufferIn(PortID port, mimmo::IBinaryStream& input);
    void    readBufferIn(PortID port);
    void    cleanBufferIn(PortID port);

    void    addParent(BaseManipulation* parent);
    void    addChild(BaseManipulation* child);
    void    unsetParent(BaseManipulation * parent);
    void    unsetChild(BaseManipulation * child);

    PortID   findPinIn(PortIn& pin);
    PortID   findPinOut(PortOut& pin);

    void    addPinIn(BaseManipulation* objIn, PortID portR);
    void    addPinOut(BaseManipulation* objOut, PortID portS, PortID portR);

    void    removePinIn(BaseManipulation* objIn, PortID portR);
    void    removePinOut(BaseManipulation* objOut, PortID portS);

    void    removePinIn(PortID portR, int j);
    void    removePinOut(PortID portS, int j);

    /*!Execution method.
     * Pure virtual method, it has to be implemented in derived classes.
     */
    virtual void    execute() = 0;
    virtual void    plotOptionalResults();
    virtual void	apply();
    void            _apply(MimmoPiercedVector<darray3E> & displacements);


    void		    write(MimmoSharedPointer<MimmoObject> geometry);

    template<typename mpv_t, typename... Args>
    void		    write(MimmoSharedPointer<MimmoObject> geometry, MimmoPiercedVector<mpv_t> & data, Args ... args);

    template<typename mpv_t>
    void		    write(MimmoSharedPointer<MimmoObject> geometry, MimmoPiercedVector<mpv_t> & data);

    template<typename mpv_t, typename... Args>
    void		    write(MimmoSharedPointer<MimmoObject> geometry, std::vector<MimmoPiercedVector<mpv_t>> & data, Args ... args);

    template<typename mpv_t>
    void		    write(MimmoSharedPointer<MimmoObject> geometry, std::vector<MimmoPiercedVector<mpv_t>> & data);

    template<typename mpv_t, typename... Args>
    void            write(MimmoSharedPointer<MimmoObject> geometry, std::vector<MimmoPiercedVector<mpv_t>*> & data, Args ... args);

    template<typename mpv_t>
    void            write(MimmoSharedPointer<MimmoObject> geometry, std::vector<MimmoPiercedVector<mpv_t>*> & data);

};


};

#include "BaseManipulation.tpp"


#endif /* __BASEMANIPULATION_HPP__ */
