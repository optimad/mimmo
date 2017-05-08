/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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
 \ *---------------------------------------------------------------------------*/

#ifndef MIMMOTYPEDEF_HH
#define MIMMOTYPEDEF_HH

#include <vector>
#include <array>
#include <map>

/*!
 * \ingroup typedefs
 * \{
 */

typedef std::vector<double>         dvector1D;   /**< mimmo custom typedef*/
typedef std::vector<float>          fvector1D;   /**< mimmo custom typedef*/
typedef std::vector<int>            ivector1D;   /**< mimmo custom typedef*/
typedef std::vector<long>           livector1D;  /**< mimmo custom typedef*/
typedef std::vector<uint32_t>       u32vector1D; /**< mimmo custom typedef*/
typedef std::vector<short int>      shivector1D; /**< mimmo custom typedef*/
typedef std::vector<bool>           bvector1D;   /**< mimmo custom typedef*/
typedef std::vector<char>           cvector1D;   /**< mimmo custom typedef*/
typedef std::vector<std::string>    svector1D;   /**< mimmo custom typedef*/

typedef std::array<double,2>        darray2E;   /**< mimmo custom typedef*/ 
typedef std::array<double,3>        darray3E;   /**< mimmo custom typedef*/
typedef std::array<double,4>        darray4E;   /**< mimmo custom typedef*/
typedef std::array<float,3>         farray3E;   /**< mimmo custom typedef*/
typedef std::array<int,3>           iarray3E;   /**< mimmo custom typedef*/
typedef std::array<int,2>           iarray2E;   /**< mimmo custom typedef*/
typedef std::array<uint32_t,3>      uarray3E;   /**< mimmo custom typedef*/

typedef std::vector<darray2E>       dvecarr2E;  /**< mimmo custom typedef*/
typedef std::vector<darray3E>       dvecarr3E;  /**< mimmo custom typedef*/
typedef std::vector<darray4E>       dvecarr4E;  /**< mimmo custom typedef*/
typedef std::vector<farray3E>       fvecarr3E;  /**< mimmo custom typedef*/
typedef std::vector<iarray3E>       ivecarr3E;  /**< mimmo custom typedef*/
typedef std::vector<iarray2E>       ivecarr2E;  /**< mimmo custom typedef*/

typedef std::vector<dvector1D>      dvector2D;  /**< mimmo custom typedef*/
typedef std::vector<fvector1D>      fvector2D;  /**< mimmo custom typedef*/


typedef std::vector<ivector1D>      ivector2D;  /**< mimmo custom typedef*/
typedef std::vector<livector1D>     livector2D; /**< mimmo custom typedef*/
typedef std::vector<shivector1D>    shivector2D;/**< mimmo custom typedef*/
typedef std::vector<bvector1D>      bvector2D;  /**< mimmo custom typedef*/
typedef std::vector<cvector1D>      cvector2D;  /**< mimmo custom typedef*/
typedef std::vector<svector1D>      svector2D;  /**< mimmo custom typedef*/

typedef std::vector< bvector2D >    bvector3D;  /**< mimmo custom typedef*/
typedef std::vector< bvector3D >    bvector4D;  /**< mimmo custom typedef*/
typedef std::vector< cvector2D >    cvector3D;  /**< mimmo custom typedef*/
typedef std::vector< cvector3D >    cvector4D;  /**< mimmo custom typedef*/
typedef std::vector< ivector2D >    ivector3D;  /**< mimmo custom typedef*/
typedef std::vector< ivector3D >    ivector4D;  /**< mimmo custom typedef*/
typedef std::vector< dvector2D >    dvector3D;  /**< mimmo custom typedef*/
typedef std::vector< dvector3D >    dvector4D;  /**< mimmo custom typedef*/
typedef std::vector< svector2D >    svector3D;  /**< mimmo custom typedef*/
typedef std::vector< svector3D >    svector4D;  /**< mimmo custom typedef*/

typedef std::array<darray2E,3>      dmatrix32E; /**< mimmo custom typedef*/
typedef std::array<darray3E,3>      dmatrix33E; /**< mimmo custom typedef*/
typedef std::array<darray4E,4>      dmatrix44E; /**< mimmo custom typedef*/
typedef std::array<uarray3E,3>      umatrix33E; /**< mimmo custom typedef*/

typedef std::array<dvector1D,3>     darr3Evec;  /**< mimmo custom typedef*/
typedef std::array<darr3Evec,3>     dmat33Evec; /**< mimmo custom typedef*/

typedef std::map<long int, int>     liimap;     /**< mimmo custom typedef*/


/*!
 * \}
 */
#endif //MIMMOTYPEDEF_HH
