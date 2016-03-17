
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

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGet(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGet(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void mimmo::pin::addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, mimmo::pin::pinSet(fset, objRec), mimmo::pin::pinGetR(fget, objSend));
	objRec->addPinIn(objSend, mimmo::pin::pinGetR(fget, objSend), mimmo::pin::pinSet(fset, objRec));
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
std::function<VAL*(void)> mimmo::pin::pinGetR(VAL* (T::*fget) (), U* obj){
	std::function<VAL*(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL)> mimmo::pin::pinSet(void (T::*fset) (VAL), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL)> mimmo::pin::pinSet(void (T::*fset) (VAL*), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

