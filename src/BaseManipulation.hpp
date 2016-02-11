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
 *	object from wich recover aome info (as number of degrees of freedom, initial displacements or other).
 *
 */
class BaseManipulation{
protected:
	//members
	dvecarr3E			m_displ			/**<Displacements of degrees of freedom. */
	uint32_t			m_ndeg			/**<Number of degrees of freedom. */
	BaseManipulation*	m_manipulator	/**<Pointer to manipulation object manipulator giving info to actual class. */
	MimmoObject*		m_geometry		/**<Pointer to manipulated geometry. */
	dvecarr3E			m_gdispl;		/**<Displacements of vertices of geometry. */

public:
	BaseManipulation();
	BaseManipulation(MimmoObject* geometry, BaseManipulation* parent = NULL);
	BaseManipulation(BaseManipulation* parent);
	~BaseManipulation();

	BaseManipulation(const BaseManipulation & other);
	BaseManipulation & operator=(const BaseManipulation & other);

	//internal methods
	uint32_t					getNDeg();
	const dvecarr3E*			getDisplacements();
	const BaseManipulation*		getManipulator();
	const MimmoObject*			getGeometry();
	const dvecarr3E*			getGeometryDisplacements();

	void	setNDeg(uint32_t ndeg);
	void	setDisplacements(dvecarr3E & displacements);
	void 	setManipulator(BaseManipulation* manipulator);
	void 	setGeometry(MimmoObject* geometry);
	void	setGeometryDisplacements(dvecarr3E & gdisplacements);

	void 	unsetManipulator();
	void	clearDisplacements();
	void	clear();

	//relationship methods
protected:
	virtual void	recoverDisplacements();   //called in exec
public:
	virtual void 	exec() = 0;

};

#endif /* __BASEMANIPULATION_HPP__ */
