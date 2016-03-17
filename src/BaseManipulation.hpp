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
#ifndef __BASEMANIPULATION_HPP__
#define __BASEMANIPULATION_HPP__

#include "MimmoObject.hpp"
#include "Info.hpp"
#include "InOut.hpp"
#include <string>
#include <functional>

class Input;
class Result;

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief BaseManipulation is the base class of any object (derived class) for manipulation of the geometry.
 *
 *	BaseManipulation is the base class used to build each generic or particular manipulation object.
 *	This base class has some common interface methods, as the base get/set methods, and one virtual method.
 *	The only method to be called to execute the manipulation object is the method exec() that calls the pure virtual execute().
 *	Each manipulation base object has a linked geometry (a target MiMMO object) and one or more linked manipulation
 *	objects from wich recovering some info (as number of degrees of freedom, initial displacements or other).
 *
 */
class BaseManipulation{
public:
	//members
	int								m_ndeg;			/**<Number of degrees of freedom used as input. */
	dvecarr3E						m_displ;		/**<Displacements of degrees of freedom used as input. */
protected:
	std::vector<BaseManipulation*>	m_parent;		/**<Pointers to manipulation objects parent giving info to current class. */
	std::vector<BaseManipulation*>	m_child;		/**<Pointers to manipulation objects child giving/receiving info (degrees of freedom and its displacements) to current class. */
	MimmoObject*					m_geometry;		/**<Pointer to manipulated geometry. */

	bool							m_relInfo;		/**<Is this object a "release Info" object?.*/
	Info*							m_info;			/**<Pointer to related object of class Info.*/
	std::vector<InOut*>				m_pinIn;		/**<Input pins vector. */
	std::vector<InOut*>				m_pinOut;		/**<Output pins vector. */

	Input*							m_input;		/**<Pointer to a bas class object Input (derived class is template).*/
	Result*							m_result;		/**<Pointer to a bas class object Result (derived class is template).*/


public:
	BaseManipulation();
	BaseManipulation(MimmoObject* geometry, BaseManipulation* child = NULL);
	BaseManipulation(BaseManipulation* child);
	~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	//get/set methods
	int					getNDeg();
	dvecarr3E&			getDisplacements();
//	int					getNDegOut(int i = 0);
	int					getNParent();
	BaseManipulation*	getParent(int i = 0);
	int					getNChild();
	BaseManipulation*	getChild(int i = 0);
	MimmoObject*		getGeometry();
	bool			 	getReleaseInfo();
	Info*			 	getInfo();

	void	setNDeg(int ndeg);
	void	setDisplacements(dvecarr3E displacements);
//	void	setNDegOut(int i, int ndeg);
//	void	setDisplacementsOut(int i, dvecarr3E & displacements);
	void 	setGeometry(MimmoObject* geometry);
	void 	setReleaseInfo(bool flag = true);

	//template get/set methods
	template<typename T>
	void 	setInput(T* data);
	template<typename T>
	void 	setInput(T& data);
	template<typename T>
	void 	setResult(T* data);
	template<typename T>
	void 	setResult(T& data);

	template<typename T>
	T* 	getInput();
	template<typename T>
	T* 	getResult();


	void 	unsetGeometry();
	void	clearDisplacements();
//	void	clearDisplacementsOut();
//	void	clearDisplacementsOut(int i = 0);
	void	clearInput();
	void	clearResult();
	void	clear();

	//TODO MAKE CHILD/PARENT ADD/UNSET METHODS PRIVATE
	void 	addParent(BaseManipulation* parent);
	void 	addChild(BaseManipulation* child);
	void 	unsetParent();
	void 	unsetChild();

	void 	removePins();
	void 	removePinsIn();
	void 	removePinsOut();
	void 	removePinIn(int i);
	void 	removePinOut(int i);

	//templated pins methods
	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal);
	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal);
	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal);

	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal);
	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal);
	template<typename T>
	void	addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void	addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal);

public:
	//relationship methods
	void 	exec();
	void	releaseInfo();
	Info* 	recoverInfo();

protected:

public:
	virtual void	setInfo();
	virtual void	useInfo();
	virtual void 	execute() = 0;				//called in exec

};

//==============================//
//TEMPLATED AND EXTERNAL METHODS
//==============================//

//==============================//
//EXTERNAL TEMPLATE METHODS
//==============================//

//==============================//
// TEMPLATED PINS METHODS		//
//==============================//

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGet(fget, objSend));
	objRec->addPinIn(objSend, pinGet(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGetR(fget, objSend));
	objRec->addPinIn(objSend, pinGetR(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGetR(fget, objSend));
	objRec->addPinIn(objSend, pinGetR(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGet(fget, objSend));
	objRec->addPinIn(objSend, pinGet(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGetR(fget, objSend));
	objRec->addPinIn(objSend, pinGetR(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGetR(fget, objSend));
	objRec->addPinIn(objSend, pinGetR(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}


template<typename T, typename U, typename VAL>
std::function<VAL(void)> pinGet(VAL (T::*fget) (), U* obj){
	std::function<VAL(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<VAL&(void)> pinGetR(VAL& (T::*fget) (), U* obj){
	std::function<VAL&(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<VAL*(void)> pinGetR(VAL* (T::*fget) (), U* obj){
	std::function<VAL*(void)> res = std::bind(fget, obj);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL)> pinSet(void (T::*fset) (VAL), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}

template<typename T, typename U, typename VAL>
std::function<void(VAL)> pinSet(void (T::*fset) (VAL*), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}


//==================================================//
// BASEMANIPULATION CLASS TEMPLATED PINS METHODS	//
//==================================================//


template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};



template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setInput(objIn, getVal, setVal);
	m_pinIn.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};

template<typename T>
void
BaseManipulation::addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal){
	InOutT<T>* pin = new InOutT<T>();
	pin->setOutput(objOut, setVal, getVal);
	m_pinOut.push_back(pin);
};



//==============================//
//RESULT BASE CLASS
//==============================//
class Result{
	template<typename T>
	T* getResult(){};
};

//==============================//
//RESULT DERIVED TEMPLATE CLASS
//==============================//
template<typename T>
class ResultT: public Result{
	T m_data;

public:
	ResultT();
	ResultT(T data){
		m_data = data;
	};
	~ResultT();

	ResultT(const ResultT & other){
		m_data 	= other.m_data;
	}

	ResultT & operator=(const ResultT & other){
		m_data 	= other.m_data;
		return (*this);
	}

	void setResult(T data){
		m_data = data;
	}

	T* getResult(){
		return(&m_data);
	}

};


//==============================//
//INPUT BASE CLASS
//==============================//
class Input{

	template<typename T>
	T* getInput(){};

};

//==============================//
//INPUT DERIVED TEMPLATE CLASS
//==============================//
template<typename T>
class InputT: public Input{
	T m_data;

public:
	InputT();
	InputT(T data){
		m_data = data;
	};
	~InputT();

	InputT(const InputT & other){
		m_data 	= other.m_data;
	}

	InputT & operator=(const InputT & other){
		m_data 	= other.m_data;
		return (*this);
	}

	void setInput(T data){
		m_data = data;
	}

	T* getInput(){
		return(&m_data);
	}

};



//==================================================//
// BASEMANIPULATION CLASS TEMPLATED INPUT METHODS	//
//==================================================//

template<typename T>
void
BaseManipulation::setInput(T* data){
	m_input = new InputT<T>(*data);
}

template<typename T>
void
BaseManipulation::setInput(T& data){
	m_input = new InputT<T>(data);
}

template<typename T>
T*
BaseManipulation::getInput(){
	return(static_cast<InputT<T>*>(m_input)->getInput());
}


//==================================================//
// BASEMANIPULATION CLASS TEMPLATED RESULT METHODS	//
//==================================================//

template<typename T>
void
BaseManipulation::setResult(T* data){
	m_result = new ResultT<T>(*data);
}

template<typename T>
void
BaseManipulation::setResult(T& data){
	m_result = new ResultT<T>(data);
}

template<typename T>
T*
BaseManipulation::getResult(){
	return(static_cast<ResultT<T>*>(m_result)->getResult());
}




#endif /* __BASEMANIPULATION_HPP__ */
