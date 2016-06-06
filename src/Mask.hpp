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

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Mask is the class that applies a geometrical masking filter to a set of data.
 *
 *	Mask is derived from BaseManipulation class.
 *	It uses the members m_coords as values to be compared with the thresholds
 *	fixed by the user for each coordinate.
 *	The flags into m_inside define, for each coordinate,
 *	if the object applies the masking inside or outside the thersholds.
 *
 *	=========================================================
 *	PORT TABLE
 *
 *	- Input  Ports-
 *
 *	PortID	PortType	variable/function	compatibilities
 *
 *	0 		DISPLS		m_displ					GDISPLS
 *	1 		COORDS		m_coords				DISPLS GDISPLS
 *	2 		RANGE		m_thres					-
 *	3 		BOOLS3		m_inside				-
 *
 *
 *	- Output Ports -
 *
 *	PortID	PortType	variable/function
 *
 *	0		DISPLS		getDisplacements
 *
 *	=========================================================
 *
 */
class Mask: public BaseManipulation{
private:
	dvecarr3E			m_coords;	/**<Coordinates of degrees of freedom of manipulator.*/
	dmatrix32E			m_thres;	/**<Limit of coordinates (min,max for each coordinate) to apply the masking.*/
	std::array<bool,3>	m_inside;	/**<Condition to apply the mask (true/false to set to zero the displacements inside/outside the thresholds).*/
	dvecarr3E			m_displ;	/**<Dispalcements of degree of freedom of manipulator.*/

public:
	Mask();
	~Mask();

	Mask(const Mask & other);
	Mask & operator=(const Mask & other);

	void	buildPorts();

	dvecarr3E	getCoords();
	dvecarr3E	getDisplacements();

	void	setCoords(dvecarr3E coords);
	void	setThresholds(dmatrix32E thres);
	void	setThresholds(darray2E thres, int dir);
	void	setMinThresholds(darray3E thres);
	void	setMaxThresholds(darray3E thres);
	void	setThresholdx(darray2E thres);
	void	setThresholdy(darray2E thres);
	void	setThresholdz(darray2E thres);
	void	setInside(bool inside);
	void	setInside(std::array<bool,3> inside);

	void	setInside(int i, bool inside);

	void 	execute();

};

}

#endif /* __MASK_HPP__ */
