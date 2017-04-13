/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __BENDGEOMETRY_HPP__
#define __BENDGEOMETRY_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class BendGeometry
 * \brief BendGeometry is the class that applies the a polynomial bending function of coordinates
 *  to the displacements of the nodes of a MimmoObject.
 *
 * The displacements to be bend have to be stored in the input of base class.
 * The bend result is stored in result member of base class.
 * For each component i the bending function of the displacement is Si = sum_jk( aijk * xj^k );
 * where aijk is the polynomial coefficient of term of degree k related to coordinate j in the function
 * applied to the i-th displacements.
 * 
 *	=========================================================
 * ~~~
 *	|---------------------------------------------------------------|
 *	|                  Port Input                                   |
 *	|-------|-----------|-------------------|-----------------------|
 *	|PortID | PortType  | variable/function | DataType              |
 *	|-------|-----------|-------------------|-----------------------|
 *  | 12    | M_FILTER  | setFilter         | (VECTOR, FLOAT)       |
 *  | 20    | M_POINT   | setOrigin         | (ARRAY3, FLOAT)       |
 *  | 22    | M_AXES    | setRefSystem      | (ARR3ARR3, FLOAT)     |
 *	| 31    | M_BMATRIX | setDegree         | (ARR3ARR3, INT)       |
 *	| 32    | M_BCOEFFS | setCoeffs         | (ARR3ARR3VEC, FLOAT)  |
 *  | 99    | M_GEOM    | setGeometry       | (SCALAR, MIMMO_)      |
 *	|-------|-----------|-------------------|-----------------------|
 *
 *
 *	|------------------------------------------------------------|
 *	|                  Port Output                               |
 *	|-------|----------|-------------------|---------------------|
 *	|PortID | PortType | variable/function | DataType            |
 *	|-------|----------|-------------------|---------------------|
 *  | 11    | M_GDISPLS| getDisplacements  | (VECARR3, FLOAT)    |
 *	|-------|----------|-------------------|---------------------|
 * ~~~
 *	=========================================================
 *
 * TODO implementation of user interface seems not ROBUST. Please recheck set*** methods.
 */
class BendGeometry: public BaseManipulation{
private:
    darray3E            m_origin;       /**<Origin of the reference system.*/
    dmatrix33E          m_system;       /**<Local reference system w.r.t absolute one.*/
    bool                m_local;        /**<True if the reference system is set by user and not the default one.*/
	umatrix33E			m_degree;	    /**<Degree of polynomial law for each coordinate
										    (each componentns of displacement is
										    f(x,y,z) with no mixed terms)*/
	dmat33Evec			m_coeffs;	    /**<Coeffs of polynomial law for each coordinate.*/
    dvector1D           m_filter;       /**<Filter field for displacements modulation. */
    dvecarr3E           m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
	BendGeometry();
	BendGeometry(const bitpit::Config::Section & rootXML);
	~BendGeometry();

	BendGeometry(const BendGeometry & other);
	BendGeometry & operator=(const BendGeometry & other);

	void buildPorts();

	umatrix33E	getDegree();
	dmat33Evec	getCoeffs();
	dvecarr3E	getDisplacements();

    void    setFilter(dvector1D filter);
    void    setOrigin(darray3E origin);
    void    setRefSystem(dmatrix33E axes);
	void	setDegree(umatrix33E degree);
	void	setDegree(int i, int j, uint32_t degree);
	void	setCoeffs(dmat33Evec coeffs);
	void	setCoeffs(int i, int j, dvector1D coeffs);

	void 	execute();

	//XML utilities from reading writing settings to file
	virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
	virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");	
	

private:
	darray3E    toLocalCoord(darray3E point);

};

REGISTER(BaseManipulation, BendGeometry, "mimmo.BendGeometry");

}

#endif /* __BENDGEOMETRY_HPP__ */
