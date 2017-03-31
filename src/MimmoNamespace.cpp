#include "MimmoNamespace.hpp"
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * Default constructor od FileDataInfo
 */
FileDataInfo::FileDataInfo(){};
	
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

/**< @namespace mimmo::pin
 * pin namespace
 */

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] portS Port ID of the output port of sender object.
 * \param[in] portR Port ID of the input port of receiver object.
 * \param[in] forced If true it forces to build the connection without checking the compatibility between ports.
 * \return True if the pin is added.
 */
bool
addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced){
	bool done = false;
	if (!objSend->arePortsBuilt()){
		objSend->buildPorts();
		if (!objSend->arePortsBuilt()){
			std::cout << "MiMMO : error " << objSend->m_name << " cannot build ports -> exit! " << std::endl;
			exit(11);
		}
	}
	if (!objRec->arePortsBuilt()){
		objRec->buildPorts();
		if (!objRec->arePortsBuilt()){
			std::cout << "MiMMO : error " << objRec->m_name << " cannot build ports -> exit! " << std::endl;
			exit(11);
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
	return done;
}

/*!It remove all pins between two objects.
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

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
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

/*!It checks the compatibility between input and output ports of two objects.
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

}//end mimmo namespace
