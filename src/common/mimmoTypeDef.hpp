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
#include <map>

typedef std::vector<double>			dvector1D; /**< MiMMO custom typedef*/
typedef std::vector<float>			fvector1D; /**< MiMMO custom typedef*/
typedef std::vector<int>			ivector1D; /**< MiMMO custom typedef*/
typedef std::vector<long>			livector1D;/**< MiMMO custom typedef*/ 
typedef std::vector<uint32_t>		u32vector1D; /**< MiMMO custom typedef*/
typedef std::vector<short int>		shivector1D; /**< MiMMO custom typedef*/
typedef std::vector<bool>			bvector1D; /**< MiMMO custom typedef*/
typedef std::vector<char>			cvector1D; /**< MiMMO custom typedef*/
typedef std::vector<std::string>	svector1D; /**< MiMMO custom typedef*/

typedef std::array<double,2>    	darray2E; /**< MiMMO custom typedef*/
typedef std::array<double,3>		darray3E; /**< MiMMO custom typedef*/
typedef std::array<double,4>		darray4E; /**< MiMMO custom typedef*/
typedef std::array<float,3>			farray3E; /**< MiMMO custom typedef*/
typedef std::array<int,3>			iarray3E; /**< MiMMO custom typedef*/
typedef std::array<int,2>			iarray2E; /**< MiMMO custom typedef*/
typedef std::array<uint32_t,3>		uarray3E; /**< MiMMO custom typedef*/

typedef std::vector<darray2E>		dvecarr2E; /**< MiMMO custom typedef*/
typedef std::vector<darray3E>		dvecarr3E; /**< MiMMO custom typedef*/
typedef std::vector<darray4E>		dvecarr4E; /**< MiMMO custom typedef*/
typedef std::vector<farray3E>		fvecarr3E; /**< MiMMO custom typedef*/
typedef std::vector<iarray3E>		ivecarr3E; /**< MiMMO custom typedef*/
typedef std::vector<iarray2E>		ivecarr2E; /**< MiMMO custom typedef*/

typedef std::vector<dvector1D>		dvector2D; /**< MiMMO custom typedef*/
typedef std::vector<fvector1D>		fvector2D; /**< MiMMO custom typedef*/


typedef std::vector<ivector1D>		ivector2D; /**< MiMMO custom typedef*/
typedef std::vector<livector1D>		livector2D; /**< MiMMO custom typedef*/
typedef std::vector<shivector1D>	shivector2D; /**< MiMMO custom typedef*/
typedef std::vector<bvector1D>		bvector2D; /**< MiMMO custom typedef*/
typedef std::vector<cvector1D>		cvector2D; /**< MiMMO custom typedef*/
typedef std::vector<svector1D>		svector2D; /**< MiMMO custom typedef*/

typedef std::vector< bvector2D >	bvector3D; /**< MiMMO custom typedef*/
typedef std::vector< bvector3D >	bvector4D; /**< MiMMO custom typedef*/

typedef std::vector< cvector2D >	cvector3D; /**< MiMMO custom typedef*/
typedef std::vector< cvector3D >	cvector4D; /**< MiMMO custom typedef*/

typedef std::vector< ivector2D >	ivector3D; /**< MiMMO custom typedef*/
typedef std::vector< ivector3D >	ivector4D; /**< MiMMO custom typedef*/

typedef std::vector< dvector2D >	dvector3D;/**< MiMMO custom typedef*/
typedef std::vector< dvector3D >	dvector4D;/**< MiMMO custom typedef*/

typedef std::vector< svector2D >	svector3D;/**< MiMMO custom typedef*/
typedef std::vector< svector3D >	svector4D;/**< MiMMO custom typedef*/

typedef std::array<darray2E,3>		dmatrix32E;/**< MiMMO custom typedef*/
typedef std::array<darray3E,3>		dmatrix33E;/**< MiMMO custom typedef*/
typedef std::array<darray4E,4>		dmatrix44E;/**< MiMMO custom typedef*/
typedef std::array<uarray3E,3>		umatrix33E;/**< MiMMO custom typedef*/

typedef std::array<dvector1D,3>		darr3Evec;/**< MiMMO custom typedef*/
typedef std::array<darr3Evec,3>		dmat33Evec;/**< MiMMO custom typedef*/

typedef std::map<long int, int>		liimap;/**< MiMMO custom typedef*/

#endif //MIMMOTYPEDEF_HH
