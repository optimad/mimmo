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

	std::vector<InOut*>				m_pinIn;			/**<Input pins vector. */
	std::vector<InOut*>				m_pinOut;			/**<Output pins vector. */

public:
	BaseManipulation();
	BaseManipulation(MimmoObject* geometry, BaseManipulation* child = NULL);
	BaseManipulation(BaseManipulation* child);
	~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	//internal methods
	int					getNDeg();
	dvecarr3E&			getDisplacements();
	int					getNDegOut(int i = 0);
	dvecarr3E*			getDisplacementsOut(int i = 0);
	int					getNParent();
	BaseManipulation*	getParent(int i = 0);
	int					getNChild();
	BaseManipulation*	getChild(int i = 0);
	MimmoObject*		getGeometry();
	bool			 	getReleaseInfo();
	Info*			 	getInfo();

	void	setNDeg(int ndeg);
	void	setDisplacements(dvecarr3E displacements);
	void	setNDegOut(int i, int ndeg);
	void	setDisplacementsOut(int i, dvecarr3E & displacements);
	void 	addParent(BaseManipulation* parent);
	void 	addChild(BaseManipulation* child);
	void 	setGeometry(MimmoObject* geometry);
	void 	setReleaseInfo(bool flag = true);

	void 	unsetParent();
	void 	unsetChild();
	void 	unsetGeometry();
	void	clearDisplacements();
	void	clearDisplacementsOut();
	void	clearDisplacementsOut(int i = 0);
	void	clear();

	//overloaded pins methods
	void	addPinIn(BaseManipulation* objIn, std::function<double(void)> getVal, std::function<void(double)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(double)> setVal, std::function<double(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<double&(void)> getVal, std::function<void(double)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(double)> setVal, std::function<double&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<int(void)> getVal, std::function<void(int)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(int)> setVal, std::function<int(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<int&(void)> getVal, std::function<void(int)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(int)> setVal, std::function<int&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<dvector1D(void)> getVal, std::function<void(dvector1D)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(dvector1D)> setVal, std::function<dvector1D(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<dvector1D&(void)> getVal, std::function<void(dvector1D)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(dvector1D)> setVal, std::function<dvector1D&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<ivector1D(void)> getVal, std::function<void(ivector1D)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(ivector1D)> setVal, std::function<ivector1D(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<ivector1D&(void)> getVal, std::function<void(ivector1D)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(ivector1D)> setVal, std::function<ivector1D&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<darray3E(void)> getVal, std::function<void(darray3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(darray3E)> setVal, std::function<darray3E(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<darray3E&(void)> getVal, std::function<void(darray3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(darray3E)> setVal, std::function<darray3E&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<iarray3E(void)> getVal, std::function<void(iarray3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(iarray3E)> setVal, std::function<iarray3E(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<iarray3E&(void)> getVal, std::function<void(iarray3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(iarray3E)> setVal, std::function<iarray3E&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<dvecarr3E(void)> getVal, std::function<void(dvecarr3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(dvecarr3E)> setVal, std::function<dvecarr3E(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<dvecarr3E&(void)> getVal, std::function<void(dvecarr3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(dvecarr3E)> setVal, std::function<dvecarr3E&(void)> getVal);

	void	addPinIn(BaseManipulation* objIn, std::function<ivecarr3E(void)> getVal, std::function<void(ivecarr3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(ivecarr3E)> setVal, std::function<ivecarr3E(void)> getVal);
	void	addPinIn(BaseManipulation* objIn, std::function<ivecarr3E&(void)> getVal, std::function<void(ivecarr3E)> setVal);
	void	addPinOut(BaseManipulation* objOut, std::function<void(ivecarr3E)> setVal, std::function<ivecarr3E&(void)> getVal);

	//relationship methods
	void 	exec();
	void	releaseInfo();
	Info* 	recoverInfo();


protected:
//	virtual void	recoverDisplacementsIn();	//TODO Useful?

public:
	virtual void	setInfo();
	virtual void	useInfo();
	virtual void 	execute() = 0;				//called in exec

};

//EXTERNAL METHODS

//EXTERNAL TEMPLATE METHODS

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL)){

	objSend->addPinOut(objRec, pinSet(fset, objRec), pinGet(fget, objSend));
	objRec->addPinIn(objSend, pinGet(fget, objSend), pinSet(fset, objRec));
	objSend->addChild(objRec);
	objRec->addParent(objSend);

}

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
std::function<void(VAL)> pinSet(void (T::*fset) (VAL), U* obj){
	std::function<void(VAL)> res = std::bind(fset, obj, std::placeholders::_1);
	return res;
}


#endif /* __BASEMANIPULATION_HPP__ */
