#include "InOut.hpp"

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (copy return).
 * \param[in] fset Set function of the receiver object (copy argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGet(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSet(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;
}

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (reference return).
 * \param[in] fset Set function of the receiver object (copy argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;

}

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (pointer return).
 * \param[in] fset Set function of the receiver object (copy argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSet(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;

}

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (copy return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGet(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;

}

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (reference return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;

}

/*!It adds a pin between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (pointer return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 * \return True if the pin is added.
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
bool mimmo::pin::addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){
	bool done = false;
	if (!(objSend->getPinType() == PinsType::BACKWARD) && !(objRec->getPinType() == PinsType::FORWARD) ){

		objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
		objRec->addPinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
		objSend->addChild(objRec);
		objRec->addParent(objSend);
		done = true;
	}
	return done;

}

/*!It binds the get function of an object.
 * \param[in] obj Pointer to BaseManipulation object.
 * \param[in] fget Get function of the object (copy return).
 */
template<typename T, typename U, typename VAL>
std::function<VAL(void)> mimmo::pin::pinGet(VAL (T::*fget) (), U* obj){
	std::function<VAL(void)> res = std::bind(fget, obj);
	return res;
}

/*!It binds the get function of an object.
 * \param[in] obj Pointer to BaseManipulation object.
 * \param[in] fget Get function of the object (reference return).
 */
template<typename T, typename U, typename VAL>
std::function<VAL&(void)> mimmo::pin::pinGetR(VAL& (T::*fget) (), U* obj){
	std::function<VAL&(void)> res = std::bind(fget, obj);
	return res;
}

/*!It binds the get function of an object.
 * \param[in] obj Pointer to BaseManipulation object.
 * \param[in] fget Get function of the object (pointer return).
 */
template<typename T, typename U, typename VAL>
std::function<VAL*(void)> mimmo::pin::pinGetP(VAL* (T::*fget) (), U* obj){
	std::function<VAL*(void)> res = std::bind(fget, obj);
	return res;
}

/*!It binds the set function of an object.
 * \param[in] obj Pointer to BaseManipulation object.
 * \param[in] fset Set function of the object (copy argument).
 */
template<typename T, typename U, typename VAL>
std::function<void(VAL)> mimmo::pin::pinSet(void (T::*fset) (VAL), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

/*!It binds the set function of an object.
 * \param[in] obj Pointer to BaseManipulation object.
 * \param[in] fset Set function of the object (pointer argument).
 */
template<typename T, typename U, typename VAL>
std::function<void(VAL*)> mimmo::pin::pinSetP(void (T::*fset) (VAL*), U* obj){
	std::function<void(VAL*)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

/*!It remove all pins between two objects.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 */
template<typename OO, typename OI>
void mimmo::pin::removeAllPins(OO* objSend, OI* objRec){

	std::vector<InOut*> pinsOut = objSend->getPinsOut();
	int removed = 0;
	for (int i=0; i<objSend->getNPinsOut(); i++){
		if (pinsOut[i]->getLink() == objRec){
			objSend->removePinOut(i-removed);
			removed++;
			objSend->unsetChild(objRec);
		}
		pinsOut[i] = NULL;
	}

	std::vector<InOut*> pinsIn = objRec->getPinsIn();
	removed = 0;
	for (int i=0; i<objRec->getNPinsIn(); i++){
		if (pinsIn[i]->getLink() == objSend){
			objRec->removePinIn(i-removed);
			removed++;
			objRec->unsetParent(objSend);
		}
		pinsIn[i] = NULL;
	}

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (copy return).
 * \param[in] fset Set function of the receiver object (copy argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (reference return).
 * \param[in] fset Set function of the receiver object (copy argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (pointer return).
 * \param[in] fset Set function of the receiver object (copy argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (copy return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (reference return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

/*!It remove a pin between two objects. If the pin is found it is removed, otherwise
 * nothing is done.
 * \param[in] objSend Pointer to BaseManipulation sender object.
 * \param[in] objRec Pointer to BaseManipulation receiver object.
 * \param[in] fget Get function of the sender object (pointer return).
 * \param[in] fset Set function of the receiver object (pointer argument).
 */
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}


