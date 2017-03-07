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
#ifndef __APPLYDEFORMATION_HPP__
#define __APPLYDEFORMATION_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{
/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Apply is the class that applies the deformation resulting from a manipulation object to the geometry.
 *
 *	Apply is derived from BaseManipulation class. It uses the base member m_geometry to apply
 *	the result of the parent manipulator to the target MiMMO object.
 *	The deformation displacements have to be passed to the input member of base class through
 *	a pin linking or set by the user, i.e. one has to use setInput method to set
 *	the displacements to be applied to the geometry.
 *	After the execution of an object Apply, the original geometry will be modified.
 *
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------|
 *	|                  Port Input                             |
 *	|-------|-----------|-------------------|-----------------|
 *	|PortID | PortType  | variable/function | dataType        |
 *	|-------|-----------|-------------------|-----------------|
 *	| 11    | M_GDISPLS | setInput          | (VECARR3,FLOAT) |
 *	| 99    | M_GEOM    | setGeometry       | (SCALAR,MIMMO_) |
 *	|-------|-----------|-------------------|-----------------|
 *
 *
 *	|--------------------------------------------------------|
 *	|                  Port Output                           |
 *	|-------|----------|-------------------|-----------------|
 *	|PortID | PortType | variable/function | dataType        |
 *	|-------|----------|-------------------|-----------------|
 *  | 99    | M_GEOM   | getGeoemtry       | (SCALAR,MIMMO_) |
 *	|-------|----------|-------------------|-----------------|
 * ~~~
 *	=========================================================
 *
 */
class Apply: public BaseManipulation{
public:

	dvecarr3E	m_input;
	bool m_force;

	Apply();
	Apply(const bitpit::Config::Section & rootXML);
	
	~Apply();

	Apply(const Apply & other);
	Apply & operator=(const Apply & other);

	void buildPorts();

	bool getRefreshGeometryTrees();
	void setRefreshGeometryTrees(bool force);

	void setInput(dvecarr3E input);

	void execute();

	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");	
};

REGISTER(BaseManipulation, Apply, "MiMMO.Apply")

}

#endif /* __APPLYDEFORMATION_HPP__ */
