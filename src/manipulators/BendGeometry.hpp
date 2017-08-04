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
\*---------------------------------------------------------------------------*/
#ifndef __BENDGEOMETRY_HPP__
#define __BENDGEOMETRY_HPP__

#include "BaseManipulation.hpp"
#include "Apply.hpp"

namespace mimmo{

/*!
 * \class BendGeometry
 * \ingroup manipulators
 * \brief BendGeometry is the class that applies the a polynomial bending function of coordinates
 *  to the displacements of the nodes of a MimmoObject.
 *
 * The displacements to be bend have to be stored in the input of base class.
 * The bend result is stored in result member of base class.
 * For each component i the bending function of the displacement is Si = sum_jk( aijk * xj^k );
 * where aijk is the polynomial coefficient of term of degree k related to coordinate j in the function
 * applied to the i-th displacements.
 * 
 * \n
 * Ports available in BendGeometry Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_FILTER  | setFilter         | (MC_VECTOR, MD_FLOAT)       |
     | M_POINT   | setOrigin         | (MC_ARRAY3, MD_FLOAT)       |
     | M_AXES    | setRefSystem      | (MC_ARR3ARR3, MD_FLOAT)     |
     | M_BMATRIX | setDegree         | (MC_ARR3ARR3, MD_INT)       |
     | M_BCOEFFS | setCoeffs         | (MC_ARR3ARR3VEC, MD_FLOAT)  |
     | M_GEOM    | setGeometry       | (MC_SCALAR, MD_MIMMO_)      |
 
     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
     | M_GDISPLS| getDisplacements  | (MC_VECARR3, MD_FLOAT)    |
     | M_GEOM   | getGeometry       | (MC_SCALAR,MD_MIMMO_) |
 
 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited form BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.BendGeometry</tt>;
 * - <B>Priority</B>  : uint marking priority of class execution in multichain frame;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 *
 * Proper of the class:
 * - <B>DegreesMatrix(3x3)</B>: degrees of each polynomial function referred to a displacement
 *                          in direction i (x,y,z) and modulating displacement in direction j (x,y,z). Degree 0
 *                          marks a constant function.
 *                          Written in XML as: \n
 *                          <tt> \<DegreesMatrix\> \n
 *                              \<xDispl\> 1 0 0 \</xDispl\> (linear x-displacement distribution in x bending direction) \n
 *                              \<yDispl\> 2 0 0 \</yDispl\> (quadratic y-displacement distribution in x bending direction) \n
 *                              \<zDispl\> 0 3 0 \</zDispl\> (cubic z-displacement in y bending direction) \n
 *                          \</DegreesMatrix\> </tt> \n
 * - <B>PolyCoefficients</B>: coefficients of each 9 bending polynomial functions. Writing following
 *                        the enumeration n = i*3 + j, where i is the displacement direction and j the bending direction.
 *                        For example n=7 corresponds to i=2, j=1, reflecting in a z displacement distribution in y-bending direction.
 *                        Please note the number of coefficients for a ij-bending function is equal to the degree
 *                        of freedom DegreesMatrix(i,j), ordered as c0 + c1*x +c2*x^2+...
 *                        Written in xml as : \n
 *                        <tt> \<PolyCoefficients\> \n
 *                          \<Poly0\> 1.0 1.5 \</Poly0\> \n
 *                          \<Poly3\> -0.1 0.2 -0.01 \</Poly3\> \n
 *                          \<Poly7\> 1.5 0.0 0.1 0.2 \</Poly7\> \n
 *                        \</PolyCoefficients\> </tt> \n
 * - <B>Origin</B>: 3D point marking the origin of reference system (if not set O = {0,0,0});
 * - <B>RefSystem</B>: axes of local reference system. written in XML as:
 *                  <tt> \<RefSystem\> \n
 *                      \<axis0\> 1.0 0.0 0.0 \</axis0\> \n
 *                      \<axis1\> 0.0 1.0 0.0 \</axis1\> \n
 *                      \<axis2\> 0.0 0.0 1.0 \</axis2\> \n
 *                  \</RefSystem\> </tt> \n
 * 
 *
 * Geometry has to be mandatorily passed through port.
 *
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
    umatrix33E    getDegree();
    dmat33Evec    getCoeffs();
    dmpvecarr3E    getDisplacements();

    void    setFilter(dmpvector1D filter);
    void    setOrigin(darray3E origin);
    void    setRefSystem(dmatrix33E axes);
    void    setDegree(umatrix33E degree);
    void    setDegree(int i, int j, uint32_t degree);
    void    setCoeffs(dmat33Evec coeffs);
    void    setCoeffs(int i, int j, dvector1D coeffs);

    void     execute();
    void     apply();

    void    checkFilter();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
protected:
    void swap(BendGeometry & x) noexcept;

private:
    darray3E    toLocalCoord(darray3E point);
    darray3E    toGlobalCoord(darray3E point);

};


//Ports
REGISTER_PORT(M_FILTER, MC_VECTOR, MD_FLOAT)
REGISTER_PORT(M_POINT, MC_ARRAY3, MD_FLOAT)
REGISTER_PORT(M_AXES, MC_ARR3ARR3, MD_FLOAT)
REGISTER_PORT(M_BMATRIX, MC_ARR3ARR3, MD_INT)
REGISTER_PORT(M_BCOEFFS, MC_ARR3ARR3VEC, MD_FLOAT)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_)
REGISTER_PORT(M_GDISPLS, MC_VECARR3, MD_FLOAT)
REGISTER_PORT(M_PAIRVECFIELD, MC_PAIR, MD_MIMMO_VECARR3FLOAT_)

//ManipBlocks
REGISTER(BaseManipulation, BendGeometry, "mimmo.BendGeometry")

}

#endif /* __BENDGEOMETRY_HPP__ */
