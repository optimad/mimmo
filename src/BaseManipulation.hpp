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
#include <string>


/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief BaseManipulation is the base class of any object (derived class) for manipulation of the geometry.
 *
 *	BaseManipulation is the base class used to build each generic or particular manipulation object.
 *	This base class has some common interface methods, as the base get/set methods, and two virtual methods.
 *	The only methods to be called to execute the manipulation object is the pure virtual method exec().
 *	Each manipulation base object has a linked geometry (a target MiMMO object) and a linked manipulation
 *	object from wich recover some info (as number of degrees of freedom, initial displacements or other).
 *
 */
class BaseManipulation{
public:
	//members
	uint32_t						m_ndeg;			/**<Number of degrees of freedom used as input. */
	dvecarr3E						m_displ;		/**<Displacements of degrees of freedom used as input. */
	uint32_t						m_ndegout;		/**<Number of degrees of freedom given as output. */
	dvecarr3E						m_displout;		/**<Displacements of degrees of freedom given as output. */
	//TODO Verify its usefulness
	BaseManipulation*				m_parent;		/**<Pointer to manipulation object parent giving info to current class. */
	std::vector<BaseManipulation*>	m_child;		/**<Pointers to manipulation objects child giving/receiving info (degrees of freedom and its displacements) to current class. */
	MimmoObject*					m_geometry;		/**<Pointer to manipulated geometry. */

public:
	BaseManipulation();
	BaseManipulation(MimmoObject* geometry, BaseManipulation* child = NULL);
	BaseManipulation(BaseManipulation* child);
	~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	//internal methods
	uint32_t			getNDeg();
	dvecarr3E*			getDisplacements();
	uint32_t			getNDegOut();
	dvecarr3E*			getDisplacementsOut();
	BaseManipulation*	getParent();
	int					getNChild();
	BaseManipulation*	getChild(int i);
	MimmoObject*		getGeometry();

	void	setNDeg(uint32_t ndeg);
	void	setDisplacements(dvecarr3E & displacements);
	void	setNDegOut(uint32_t ndeg);
	void	setDisplacementsOut(dvecarr3E & displacements);
	void 	setParent(BaseManipulation* parent);
	void 	addChild(BaseManipulation* child);
	void 	setGeometry(MimmoObject* geometry);

	void 	unsetParent();
	void 	unsetChild();
	void 	unsetGeometry();
	void	clearDisplacements();
	void	clearDisplacementsOut();
	void	clear();

	//relationship methods

	void 	exec();

protected:
	virtual void	recoverDisplacementsIn();	//TODO Useful?
	virtual void	recoverDisplacementsOut();	//called in exec
	virtual void	initChild();				//called in exec
	virtual void	updateChild();				//called in exec
public:
	virtual void 	execute() = 0;				//called in exec

};

#endif /* __BASEMANIPULATION_HPP__ */
