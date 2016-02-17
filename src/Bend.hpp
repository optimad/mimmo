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
#ifndef __BEND_HPP__
#define __BEND_HPP__

#include "BaseManipulation.hpp"

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Bend is the class that applies the a polynomial bending to the displacements of a manipulation object.
 *
 */
class Bend: public BaseManipulation{
private:
	//members
	dvecarr3E			m_coords;	/**<Coordinates of degrees of freedom of manipulator.*/
	dvecarr3E			m_degree;	/**<Degree of polynom for each coordinate (each componentns of displacement is f(x,y,z) with no mixed terms)*/
	dvector3D			m_coeffs;	/**<Coeffs of polynom for each coordinate.*/

public:
	Bend();
	~Bend();

	Bend(const Bend & other);
	Bend & operator=(const Bend & other);

	void	setCoords(dvecarr3E & coords);
	void	setDegree(dvecarr3E & degree);
	void	setCoeffs(dvector3D & coeffs);


	//relationship methods
protected:

public:
	void 	execute();

};

#endif /* __BEND_HPP__ */
