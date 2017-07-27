#include "MimmoNamespace.hpp"
#include "BaseManipulation.hpp"

std::string mimmo::MIMMO_LOG_FILE = "mimmo"; /**<Default name of logger file.*/
bool        mimmo::MIMMO_EXPERT = false;    /**<Flag that defines expert mode (true) or safe mode (false).
                                                In case of expert mode active the mandatory ports are not checked. */

namespace mimmo{

/*!
 * Default constructor od FileDataInfo
 */
FileDataInfo::FileDataInfo(){
    ftype = 0;
};

/*!
 * Default destructor of FileDataInfo.
 */
FileDataInfo::~FileDataInfo(){};

/*!
 * Copy constructor of FileDataInfo;
 */
FileDataInfo::FileDataInfo(const FileDataInfo & other){
    *this = other;
};

/*!
 * Assignement operator of FileDataInfo.
 */
FileDataInfo & FileDataInfo::operator=(const FileDataInfo & other){
    ftype = other.ftype;
    fname = other.fname;
    fdir = other.fdir;
    return *this;
};

namespace pin{

/*!
 * It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] portS Port ID of the output port of sender object.
 * \param[in] portR Port ID of the input port of receiver object.
 * \param[in] forced If true it forces to build the connection without checking the compatibility between ports.
 * \return true if the pin is added.
 */
bool
addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced){
    bitpit::Logger* log = &bitpit::log::cout(MIMMO_LOG_FILE);
    bool done = false;
    if (!objSend->arePortsBuilt()){
        objSend->buildPorts();
        if (!objSend->arePortsBuilt()){
            (*log) << "error: " << objSend->m_name << " cannot build ports -> exit! " << std::endl;
            throw std::runtime_error (objSend->m_name + " cannot build ports");
        }
    }
    if (!objRec->arePortsBuilt()){
        objRec->buildPorts();
        if (!objRec->arePortsBuilt()){
            (*log) << "error: " << objRec->m_name << " cannot build ports -> exit! " << std::endl;
            throw std::runtime_error (objRec->m_name + " cannot build ports");
        }
    }
    if (!(objSend->getConnectionType() == ConnectionType::BACKWARD) && !(objRec->getConnectionType() == ConnectionType::FORWARD) ){
        if (objSend->m_portOut.count(portS) != 0 && objRec->m_portIn.count(portR) != 0){
            if (forced || checkCompatibility(objSend, objRec, portS, portR)){
                objSend->addPinOut(objRec, portS, portR);
                objRec->addPinIn(objSend, portR);
                objSend->addChild(objRec);
                objRec->addParent(objSend);
                done = true;
            }
        }
    }
    log->setPriority(bitpit::log::NORMAL);
    if (!done) (*log)<<"warning: pin connection " << objSend->getName() << "[" << portS << "] --> " << objRec->getName() << "[" << portR << "] NOT linked. "<< std::endl;
    log->setPriority(bitpit::log::DEBUG);
    return done;
}

/*!
 * It removes all pins between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 */
void
removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec){

    std::map<PortID, PortOut*> pinsOut = objSend->getPortsOut();
    for (std::map<PortID, PortOut*>::iterator i = pinsOut.begin(); i != pinsOut.end(); i++){
        if (i->second != NULL){
            std::vector<BaseManipulation*> linked = i->second->getLink();
            for (int j=0; j<(int)linked.size(); j++){
                if (linked[j] == objRec){
                    objSend->removePinOut(i->first,j);
                    objSend->unsetChild(objRec);
                }
            }
        }
    }

    std::map<PortID, PortIn*> pinsIn = objRec->getPortsIn();
    for (std::map<PortID, PortIn*>::iterator i = pinsIn.begin(); i != pinsIn.end(); i++){
        std::vector<BaseManipulation*> linked = i->second->getLink();
        if (i->second != NULL){
            for (int j=0; j<(int)linked.size(); j++){
                if (linked[j] == objSend){
                    objRec->removePinIn(i->first,j);
                    objRec->unsetParent(objSend);
                }
            }
        }
    }

}

/*!
 * It removes a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] portS Port ID of the output port of sender object.
 * \param[in] portR Port ID of the input port of receiver object.
 */
void
removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR){

    objSend->removePinOut(objRec, portS);
    objRec->removePinIn(objSend, portR);
    objSend->unsetChild(objRec);
    objRec->unsetParent(objSend);

}

/*!
 * It checks the compatibility between input and output ports of two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] portS Port ID of the output port of sender object.
 * \param[in] portR Port ID of the input port of receiver object.
 * \return True if the pin is added.
 */
bool
checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR){
    bool check = false;
    PortIn*		pinin 	= objRec->getPortsIn()[portR];
    PortOut*	pinout 	= objSend->getPortsOut()[portS];
    DataType typeR = pinin->getDataType();
    DataType typeS = pinout->getDataType();
    check = (typeS == typeR);
    return(check);
}

}//end pin namespace


/*!Change the name of the logger.
 * \param[in] log name of the logger file.
 */
void    setLogger(std::string log){
    if (BaseManipulation::sm_baseManipulationCounter > 1){
        bitpit::Logger* log = &bitpit::log::cout(MIMMO_LOG_FILE);
        (*log) << "warning: logger already set -> cannot change logger name " << std::endl;
        return;
    }
    MIMMO_LOG_FILE = log;
}

/*!Base warning for no data found in xml dictionary.
 * \param[in] log pointer to logger file.
 * \param[in] name name of the mimmo block.
 */
void    warningXML(bitpit::Logger* log, std::string name){
    (*log)<<"warning in custom xml " << name << " constructor. No valid xml data found"<<std::endl;
}

/*!
 * Active/inactive the expert mode for mimmo.
 */
void setExpertMode(bool flag){
    MIMMO_EXPERT = flag;
}


/*!Maximum value function for MimmoPiercedVector<double>.
 * \param[in] field MimmoPiercedVector<double>
 * \return Maximum value
 */
double  maxvalmp(const MimmoPiercedVector<double> & field){
    double val = 1.0e-18;
    for (const auto v : field){
        val = std::max(val,v);
    }
    return val;
}

}//end mimmo namespace
