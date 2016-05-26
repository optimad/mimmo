#include "MimmoNamespace.hpp"
#include "BaseManipulation.hpp"

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (copy return).
 * \param[in] fset Set function of the receiver object (copy argument).
 * \return True if the pin is added.
 */
bool
mimmo::pin::addPin(mimmo::BaseManipulation* objSend, BaseManipulation* objRec, int portS, int portR){
	bool done = false;
	if (!objSend->isPinSet()) objSend->setPins();
	if (!objRec->isPinSet()) objRec->setPins();
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){
		objSend->addPinOut(objRec, portS, portR);
		objRec->addPinIn(objSend, portR);
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;
}
