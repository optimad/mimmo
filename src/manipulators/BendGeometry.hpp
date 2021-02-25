/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
 * \ingroup manipulators
 * \brief BendGeometry applies custom bending deformations along axis-directions
   of a target geometry.
 *
 * Bending deformation is defined as a vector field of displacements along
   a certain axis (0-x,1-y or 2-z) and modelized by a set of custom polynomial functions
   of N degree. Each displacement coordinate, for each axis can be influenced
   by its own bending (9 possible bending functions). In particular, for each
   displacement component \f$i\in[0,2]\f$ the bending is expressed as:

   \f$ S_i = \sum_{j,k}( a_{ijk} x_{j}^k ) \f$

   where \f$a_{ijk}\f$ is the polynomial coefficient of term of degree
   \f$k\in[0,N]\f$, related to axis \f$j\in[0,2]\f$ and applied to the \f$i\f$-th
   displacement coordinate.
 *
 * \n
 * Ports available in BendGeometry Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_FILTER  | setFilter         | (MC_SCALAR, MD_MPVECFLOAT_)       |
     | M_POINT   | setOrigin         | (MC_ARRAY3, MD_FLOAT)       |
     | M_AXES    | setRefSystem      | (MC_ARR3ARR3, MD_FLOAT)     |
     | M_BMATRIX | setDegree         | (MC_ARR3ARR3, MD_INT)       |
     | M_BCOEFFS | setCoeffs         | (MC_SCALAR, MD_MATRIXCOEFF_)  |
     | M_GEOM    | setGeometry       | (MC_SCALAR, MD_MIMMO_)      |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | M_GDISPLS| getDisplacements  | (MC_SCALAR, MD_MPVECARR3FLOAT_)    |
     | M_GEOM   | getGeometry       | (MC_SCALAR,MD_MIMMO_) |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited form BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.BendGeometry</tt>;
 * - <B>Priority</B>  : uint marking priority of class execution in multichain frame;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry
    directly in execution;
 *
 * Proper of the class:
 * - <B>DegreesMatrix(3x3)</B>: degrees \f$d_{ij}\f$ of each polynomial function
     referred to a displacement component \f$i\in[0,2]\f$ and modulating
     it along axis \f$j\in[0,2]\f$. \f$d_{ij} = 0\f$ marks a constant
     function. Written in XML as: \n\n
      <tt><B>\<DegreesMatrix\></B> \n
      &nbsp; &nbsp; &nbsp;<B>\<xDispl\></B> 1 0 0 <B>\</xDispl\></B> (linear x-displacement distribution in x bending direction) \n
      &nbsp; &nbsp; &nbsp;<B>\<yDispl\></B> 2 0 0 <B>\</yDispl\></B> (quadratic y-displacement distribution in x bending direction) \n
      &nbsp; &nbsp; &nbsp;<B>\<zDispl\></B> 0 3 0 <B>\</zDispl\></B> (cubic z-displacement in y bending direction) \n
      <B>\</DegreesMatrix\></B> </tt> \n\n

   - <B>PolyCoefficients</B>: coefficients of each 9 bending polynomial functions.
     Writing following the enumeration \f$n = 3i + j\f$, where i is the displacement
     component and j the axis of bending. For example, n=7 corresponds to i=2, j=1,
     reflecting in a z-component displacement distribution in y-axis bending direction.
     Please note the number of coefficients for a ij-bending function is equal to the function degree
     \f$d_{ij}\f$ + 1 (DegreesMatrix), ordered as \f$c_0 + c_1x +c_2x^2+...+c_nx^n\f$.
     If a polynomial function is not specified, it will be treated as 0 constant function.
     Written in XML as : \n\n
     <tt> <B>\<PolyCoefficients\></B> \n
      &nbsp; &nbsp; &nbsp;<B>\<Poly0\></B>  1.0 1.5           <B>\</Poly0\></B>  \n
      &nbsp; &nbsp; &nbsp;<B>\<Poly3\></B> -0.1 0.2 -0.01     <B>\</Poly3\></B>  \n
      &nbsp; &nbsp; &nbsp;<B>\<Poly7\></B>  1.5 0.0   0.1 0.2 <B>\</Poly7\></B>  \n
      <B>...</B> \n
     <B>\</PolyCoefficients\></B> </tt> \n\n

   - <B>Origin</B>: 3D point marking the origin of reference system (default {0,0,0});

   - <B>RefSystem</B>: define a bending custom reference system of axes.
      Written in XML as:\n\n\
      <tt><B>\<RefSystem\></B>\n
      &nbsp; &nbsp; &nbsp;<B>\<axis0\></B> 1.0 0.0 0.0 <B>\</axis0\></B> \n
      &nbsp; &nbsp; &nbsp;<B>\<axis1\></B> 0.0 1.0 0.0 <B>\</axis1\></B> \n
      &nbsp; &nbsp; &nbsp;<B>\<axis2\></B> 0.0 0.0 1.0 <B>\</axis2\></B> \n
      <B>\</RefSystem\></B> </tt> \n\n

  Geometry has to be mandatorily passed through port.

 */
class BendGeometry: public BaseManipulation{
private:
    darray3E            m_origin;       /**<Origin of the reference system.*/
    dmatrix33E          m_system;       /**<Local reference system w.r.t absolute one.*/
    bool                m_local;        /**<True if the reference system is set by user and not the default one.*/
    umatrix33E          m_degree;        /**<Degree of polynomial law for each coordinate
                                            (each componentns of displacement is
                                            f(x,y,z) with no mixed terms)*/
    dmat33Evec          m_coeffs;        /**<Coeffs of polynomial law for each coordinate.*/
    dmpvector1D           m_filter;       /**<Filter field for displacements modulation. */
    dmpvecarr3E           m_displ;        /**<Resulting displacements of geometry vertex.*/

public:
    BendGeometry();
    BendGeometry(const bitpit::Config::Section & rootXML);
    ~BendGeometry();

    BendGeometry(const BendGeometry & other);
    BendGeometry & operator =(BendGeometry other);
    void buildPorts();

    darray3E    getOrigin();
    dmatrix33E  getRefSystem();
    umatrix33E     getDegree();
    dmat33Evec*    getCoeffs();
    dmpvecarr3E*   getDisplacements();

    void    setFilter(dmpvector1D *filter);
    void    setOrigin(darray3E origin);
    void    setRefSystem(dmatrix33E axes);
    void    setDegree(umatrix33E degree);
    void    setDegree(int i, int j, uint32_t degree);
    void    setCoeffs(dmat33Evec *coeffs);
    void    setCoeffs(int i, int j, dvector1D coeffs);

    void     execute();
    void     apply();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(BendGeometry & x) noexcept;
    void    checkFilter();
    static darray3E   matmul(const darray3E & vec, const dmatrix33E & mat);
    static darray3E   matmul(const dmatrix33E & mat,const darray3E & vec);

private:
    darray3E    toLocalCoord(darray3E point);
    darray3E    toGlobalCoord(darray3E point);

};

REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_BMATRIX, MC_ARR3ARR3, MD_INT,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_BCOEFFS, MC_SCALAR, MD_MATRIXCOEFF_,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__BENDGEOMETRY_HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_,__BENDGEOMETRY_HPP__)


REGISTER(BaseManipulation, BendGeometry, "mimmo.BendGeometry")

}

#endif /* __BENDGEOMETRY_HPP__ */
