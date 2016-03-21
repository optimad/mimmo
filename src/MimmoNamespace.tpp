
template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}


template<typename T, typename U, typename VAL>
std::function<VAL(void)> mimmo::pin::pinGet(VAL (T::*fget) (), U* obj){
	std::function<VAL(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<VAL&(void)> mimmo::pin::pinGetR(VAL& (T::*fget) (), U* obj){
	std::function<VAL&(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<VAL*(void)> mimmo::pin::pinGetP(VAL* (T::*fget) (), U* obj){
	std::function<VAL*(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL)> mimmo::pin::pinSet(void (T::*fset) (VAL), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL*)> mimmo::pin::pinSetP(void (T::*fset) (VAL*), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}


template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removeAllPins(OO* objSend, OI* objRec){

	vector<Inout*> pinsOut = objSend->getPinsOut();
	int removed = 0;
	for (int i=0; i<objSend->getNPinsOut(); i++){
		if (pinsOut[i]->getLink() == objRec){
			objSend->removePinOut(i-removed);
			removed++;
			objSend->unsetChild(objRec);
		}
		pinsOut[i] = NULL;
	}

	vector<Inout*> pinsIn = objRec->getPinsIn();
	int removed = 0;
	for (int i=0; i<objRec->getNPinsIn(); i++){
		if (pinsIn[i]->getLink() == objSend){
			objRec->removePinIn(i-removed);
			removed++;
			objRec->unsetParent(objSend);
		}
		pinsIn[i] = NULL;
	}

}


template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){

	objSend->removePinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::removePin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->removePinOut(objRec, mimmo::pin::pinSetP(fset, objRec), mimmo::pin::pinGetP(fget, objSend));
	objRec->removePinIn(objSend, mimmo::pin::pinGetP(fget, objSend), mimmo::pin::pinSetP(fset, objRec));
	objSend->unsetChild(objRec);
	objRec->unsetParent(objSend);

}


