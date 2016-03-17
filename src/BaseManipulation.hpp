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

	//friendship declaration
	template<typename OO, typename G, typename OI, typename S, typename VAL>
	friend void addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL)) ;
	
	
	
protected:
	MimmoObject*					m_geometry;		/**<Pointer to manipulated geometry. */
	std::vector<BaseManipulation*>	m_parent;		/**<Pointers to manipulation objects parent giving info to current class. */
	std::vector<BaseManipulation*>	m_child;		/**<Pointers to manipulation objects child giving/receiving info (degrees of freedom and its displacements) to current class. */

	std::vector<InOut*>				m_pinIn;		/**<Input pins vector. */
	std::vector<InOut*>				m_pinOut;		/**<Output pins vector. */

	//TODO 1)intended for providing a slot to store a generic data input/output. Maybe is not necessary. if not, please clean it
	Input*							m_input;		/**<Pointer to a base class object Input (derived class is template).*/
	Result*							m_result;		/**<Pointer to a base class object Result (derived class is template).*/

	//TODO 2)maybe can be TEMPORARY STRUCTURE
	int								m_ndeg;			/**<Number of degrees of freedom used as input. */
	dvecarr3E						m_displ;		/**<Displacements of degrees of freedom used as input. */

	//TODO 3)Info exchange can be useless later. Check it out and clean eventually
	bool							m_relInfo;		/**<Is this object a "release Info" object?.*/
	Info*							m_info;			/**<Pointer to related object of class Info.*/
	

public:
	BaseManipulation();
	~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	//get methods
	MimmoObject*		getGeometry();
	int					getNParent();
	BaseManipulation*	getParent(int i = 0);
	int					getNChild();
	BaseManipulation*	getChild(int i = 0);
	int 				getNPinIn();
	int 				getNPinOut();

	//TODO see 1)
	template<typename T>	
	T*					getInput();
	template<typename T>
	T* 					getResult();
	
	//TODO see 2)	
	int					getNDeg();
	dvecarr3E&			getDisplacements();
	
	//TODO see 3)
	bool				getReleaseInfo();
	Info*				getInfo();

	//set methods
	void 				setGeometry(MimmoObject* geometry);

	
	//TODO see 1)
	template<typename T>
	void 				setInput(T* data);
	template<typename T>
	void 				setInput(T& data);
	template<typename T>
	void 				setResult(T* data);
	template<typename T>
	void 				setResult(T& data);
	
	
	//TODO see 2)
	void				setNDeg(int ndeg);
	void				setDisplacements(dvecarr3E displacements);
	
	//TODO see 3)
	void 				setReleaseInfo(bool flag = true);
	virtual void		setInfo();
	

	//cleaning/unset/remove
	//TODO 5) Not complete list of cleaning methods. To be reviewed and reorganized
	void 	unsetGeometry();
	void	clearDisplacements();
	void	clearInput();
	void	clearResult();
	void	clear();
	void 	unsetParent();
	void 	unsetChild();
	void 	removePins();
	void 	removePinsIn();
	void 	removePinsOut();
	void 	removePinIn(int i);
	void 	removePinOut(int i);

	//execution utils
	void 	exec();

	//TODO see 3)
	void	releaseInfo();
	Info* 	recoverInfo();
	virtual void	useInfo();
	
	
protected:
	//TODO 4) public but not intended to user interface
	void				addParent(BaseManipulation* parent); 
	void				addChild(BaseManipulation* child);
	
	//TODO 4)
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T(void)> getVal);
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T&(void)> getVal);
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T)> setVal, std::function<T*(void)> getVal);
	
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T(void)> getVal);
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T&(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T&(void)> getVal);
	template<typename T>
	void				addPinIn(BaseManipulation* objIn, std::function<T*(void)> getVal, std::function<void(T*)> setVal);
	template<typename T>
	void				addPinOut(BaseManipulation* objOut, std::function<void(T*)> setVal, std::function<T*(void)> getVal);
	
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
