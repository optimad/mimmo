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
#ifndef __MASK_HPP__
#define __MASK_HPP__

#include "BaseManipulation.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Mask is the class that applies the masking filterto a manipulation object.
 *
 *	Mask is derived from BaseManipulation class. It uses the base members m_parent and m_geometry to recover
 *	the result of the parent manipulator and apply the masking.
 *	After the execution of an object Mask (called by the linked manipulator), the original displacements will be modified.
 *
 */
class Mask: public BaseManipulation{
private:
	//members
	dvecarr3E			m_coords;	/**<Coordinates of degrees of freedom of manipulator.*/
	darray3E			m_thres;	/**<Limit of coordinates to apply the masking.*/
	std::array<bool,3>	m_forward;	/**<Condition to apply the mask (true/false to set to zero the displacements >/< the thershold).*/

public:
	Mask();
	~Mask();

	Mask(const Mask & other);
	Mask & operator=(const Mask & other);

	void	setCoords(dvecarr3E & coords);
	void	setThresholds(darray3E & thres);
	void	setForward(int i, bool forward);


	//relationship methods
protected:

public:
	void 	execute();

};

#endif /* __MASK_HPP__ */
