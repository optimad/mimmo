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
#ifndef __DEFORMBOX_HPP__
#define __DEFORMBOX_HPP__

#include "BaseManipulation.hpp"
#include "FFDLattice.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief DeformBox is the class that applies the a deformation to the box of a latticeBox object.
 *
 *	The degrees of freedom are the components of the transformation matrix (ndeg = 3).
 *	The displacements are the row of the transformation matrix (matrix 3x3).
 *
 *
 */
class DeformBox: public BaseManipulation{
private:
	//members
	dvecarr3E			m_coords;	/**<Coordinates of degrees of freedom of manipulator.*/

public:
	DeformBox();
	~DeformBox();

	DeformBox(const DeformBox & other);
	DeformBox & operator=(const DeformBox & other);


	//relationship methods
protected:

public:
	void 	useInfo();
	void 	execute();

};

#endif /* __DEFORMBOX_HPP__ */
