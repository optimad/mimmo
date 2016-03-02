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
 \ *---------------------------------------------------------------------------*/

#ifndef MIMMOTYPEDEF_HH
#define MIMMOTYPEDEF_HH

#include <vector>
#include <array>

typedef std::vector<double>			dvector1D;
typedef std::vector<float>			fvector1D;
typedef std::vector<int>			ivector1D;
typedef std::vector<uint32_t>		u32vector1D;
typedef std::vector<short int>		shivector1D;
typedef std::vector<bool>			bvector1D;
typedef std::vector<char>			cvector1D;
typedef std::vector<std::string>	svector1D;

typedef std::array<double,2>    	darray2E;
typedef std::array<double,3>		darray3E;
typedef std::array<double,4>		darray4E;
typedef std::array<float,3>			farray3E;

typedef std::vector<darray2E>		dvecarr2E;
typedef std::vector<darray3E>		dvecarr3E;
typedef std::vector<darray4E>		dvecarr4E;
typedef std::vector<farray3E>		fvecarr3E;

typedef std::vector<dvector1D>		dvector2D;
typedef std::vector<fvector1D>		fvector2D;


typedef std::vector<ivector1D>		ivector2D;
typedef std::vector<shivector1D>	shivector2D;
typedef std::vector<bvector1D>		bvector2D;
typedef std::vector<cvector1D>		cvector2D;
typedef std::vector<svector1D>		svector2D;

typedef std::vector< bvector2D >	bvector3D;
typedef std::vector< bvector3D >	bvector4D;

typedef std::vector< cvector2D >	cvector3D;
typedef std::vector< cvector3D >	cvector4D;

typedef std::vector< ivector2D >	ivector3D;
typedef std::vector< ivector3D >	ivector4D;

typedef std::vector< dvector2D >	dvector3D;
typedef std::vector< dvector3D >	dvector4D;

typedef std::vector< svector2D >	svector3D;
typedef std::vector< svector3D >	svector4D;

typedef std::array<darray2E,3>		dmatrix32E;
typedef std::array<darray3E,3>		dmatrix33E;

#endif //MIMMOTYPEDEF_HH
