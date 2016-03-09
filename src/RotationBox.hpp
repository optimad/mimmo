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
#ifndef __ROTATIONBOX_HPP__
#define __ROTATIONBOX_HPP__

#include "BaseManipulation.hpp"
#include "FFDLattice.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief RotationBox is the class that applies the a rotation to the box of a latticeBox object.
 *
 *	The degrees of freedom is the rotation value (ndeg = 1) and the direction of the rotation
 *	is a parameter of the object.
 *	The displacements has to be one term. As the displacements are array3E only the first term is used,
 *	i.e if alpha is the value of the rotation in the chosen direction the only array of displacement
 *	is set m_displ[0]={aplha, 0, 0};
 *
 */
class RotationBox: public BaseManipulation{
private:
	//members
	darray3E			m_origin;		/**<origin of the rotation axis.*/
	darray3E			m_direction;	/**<Components of the rotation axis.*/
	dvecarr3E			m_axes;			/**<Axes of box to be deformed (recovered in recoverInfo and used in useInfo).*/


public:
	RotationBox(darray3E origin = { {0, 0, 0} }, darray3E direction = { {0, 0, 0} });
	~RotationBox();

	RotationBox(const RotationBox & other);
	RotationBox & operator=(const RotationBox & other);

	void setAxis(darray3E origin, darray3E direction);
	void setOrigin(darray3E origin);
	void setDirection(darray3E direction);
	void setRotation(double alpha);
private:
	dmatrix44E matMul(dmatrix44E &, dmatrix44E &);

	//relationship methods
protected:

public:
	void 	useInfo();
	void 	execute();

};

#endif /* __ROTATIONBOX_HPP__ */
