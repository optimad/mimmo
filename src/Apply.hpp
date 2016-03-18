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

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Apply is the class that applies the deformation resulting from a manipulation object to the geometry.
 *
 *	Apply is derived from BaseManipulation class. It uses the base member m_geometry to apply
 *	the result of the parent manipulator to the target MiMMO object.
 *	After the execution of an object Apply, the original geometry will be modified.
 *
 */
class Apply: public BaseManipulation{

public:
	Apply();
	~Apply();

	Apply(const Apply & other);
	Apply & operator=(const Apply & other);

	void 	execute();

};

#endif /* __APPLYDEFORMATION_HPP__ */
