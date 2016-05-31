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

namespace mimmo{

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Bend is the class that applies the a polynomial bending function of coordinates
 *	 to the displacements of a manipulation object.
 *
 *	The displacements to be bend have to be stored in the input of base class.
 *	The bend result is stored in result member of base class.
 *	For each component i the bending function of the displacement is Si = sum_jk( aijk * xj^k );
 *	where aijk is the polynomial coefficient of term of degree k related to coordinate j in the function
 *	applied to the i-th displacements.
 *
 * TODO implementation of user interface seems not ROBUST. Please recheck set*** methods.
 */
class Bend: public BaseManipulation{
private:
	dvecarr3E			m_coords;	/**<Coordinates of degrees of freedom of manipulator.*/
	umatrix33E			m_degree;	/**<Degree of polynomial law for each coordinate
										(each componentns of displacement is
										f(x,y,z) with no mixed terms)*/
	dmat33Evec			m_coeffs;	/**<Coeffs of polynomial law for each coordinate.*/
	dvecarr3E			m_displ;

public:
	Bend();
	~Bend();

	Bend(const Bend & other);
	Bend & operator=(const Bend & other);

	dvecarr3E*	getCoords();
	umatrix33E	getDegree();
	dmat33Evec	getCoeffs();

	void	setCoords(dvecarr3E coords);
	void	setDegree(umatrix33E degree);
	void	setDegree(int i, int j, uint32_t degree);
	void	setCoeffs(dmat33Evec coeffs);
	void	setCoeffs(int i, int j, dvector1D coeffs);

	void 	execute();

};

}

#endif /* __BEND_HPP__ */
