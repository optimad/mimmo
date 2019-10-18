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
#ifndef __MRBF_HPP__
#define __MRBF_HPP__

#include "BaseManipulation.hpp"
#include <bitpit_RBF.hpp>

namespace mimmo{

/*!
 * @ingroup manipulators
 * @brief Enum class defining types of mimmo RBF kernel
   functions that could be used also in MRBF class, together with
   those already available in bitpit::RBF mother class
 */
enum class MRBFBasisFunction {
    HEAVISIDE10    = 101,  /**< Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 10 */
    HEAVISIDE50    = 102,  /**< Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 50 */
	HEAVISIDE100    = 103,  /**< Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 100 */
	HEAVISIDE1000    = 104,  /**< Non compact sharp heaviside 0.5*(1.+tanh(k*x)) with k = 1000 */
};

/*!
 * \ingroup manipulators
 * \brief Solver enum for MRBF data fields interpolation/ direct parameterization
 */
enum class MRBFSol{
    NONE = 0,     /**< activate class as pure parameterizator. Set freely your RBF coefficients/weights */
    WHOLE = 1,    /**< activate class as pure interpolator, with RBF coefficients evaluated solving a full linear system for all active nodes.*/
    GREEDY= 2   /**< activate class as pure interpolator, with RBF coefficients evaluated using a greedy algorithm on active nodes.*/
};

/*!
 * \class MRBF
 * \ingroup manipulators
 * \brief Radial Basis Function evaluation from clouds of control points.
 *
 * This class is derived from BaseManipulation class of mimmo and from RBF class
 * of bitpit library.
 * It evaluates the result of RBF functions built over a set of control point
   given externally or stored in a PointCloud MimmoObject container.
   Displacements (DOFs) for RBF nodes are provided externally.
   Default solver in execution is MRBFSol::NONE for direct parameterization.
   Use MRBFSol::GREEDY or MRBFSol::WHOLE to activate interpolation features.
 * See bitpit::RBF docs for further information.
 *
 * \n
 * Ports available in MRBF Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_COORDS  | setNode               | (MC_VECARR3, MD_FLOAT)      |
     | M_DISPLS  | setDisplacements      | (MC_VECARR3, MD_FLOAT)      |
     | M_FILTER  | setFilter             | (MC_SCALAR, MD_MPVECFLOAT_)       |
     | M_VALUED  | setSupportRadius      | (MC_SCALAR, MD_FLOAT)       |
     | M_VALUED2 | setSupportRadiusValue | (MC_SCALAR, MD_FLOAT)       |
     | M_GEOM    | setGeeometry            | (MC_SCALAR, MD_MIMMO_)      |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GDISPLS      | getDisplacements  | (MC_SCALAR, MD_MPVECARR3FLOAT_)             |
     | M_GEOM   | getGeometry       | (MC_SCALAR, MD_MIMMO_) |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.MRBF</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 *
 * Proper of the class:
 * - <B>Mode</B>: mode of usage of the class see MRBFSol enum;
 * - <B>SupportRadius</B>: local radius of RBF function for each nodes, expressed as ratio of local geometry bounding box;
 * - <B>SupportRadiusReal</B>: local effective radius of RBF function common to each RBF node;
 * - <B>RBFShape</B>: shape of RBF function see MRBFBasisFunction and bitpit::RBFBasisFunction enums;
 * - <B>Tolerance</B>: greedy engine tolerance (meaningful for Mode 2 only);
 *
 * Geometry, filter field, RBF nodes and displacements have to be mandatorily passed through port.
 *
 */
//TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF
class MRBF: public BaseManipulation, public bitpit::RBF {

private:
    double       m_tol;          /**< Tolerance for greedy algorithm.*/
    MRBFSol      m_solver;       /**<Type of solver specified for the class as default in execution*/
    dmpvector1D  m_filter;       /**<Filter field for displacements modulation */
    bool         m_bfilter;      /**<boolean to recognize if a filter field is applied */
    double       m_SRRatio;      /**<support Radius ratio */
    dmpvecarr3E  m_displ;        /**<Resulting displacements of geometry vertex.*/
    bool         m_supRIsValue;  /**<True if support radius is defined as absolute value, false if is ratio of bounding box diagonal.*/
    int          m_functype;     /**< Function type handler. If -1 refer to RBF getFunctionType method */

public:
    MRBF();
    MRBF(const bitpit::Config::Section & rootXML);

    virtual ~MRBF();

    //copy operators/constructors
    MRBF(const MRBF & other);
    MRBF& operator=(MRBF other);

    void buildPorts();

    dvecarr3E*      getNodes();

    MRBFSol         getMode();
    void            setMode(MRBFSol);
    void            setMode(int);
    dmpvector1D*    getFilter();
    double          getSupportRadius();
    double          getSupportRadiusValue();
    bool            getIsSupportRadiusValue();

    int             getFunctionType();

    dmpvecarr3E*     getDisplacements();

    int             addNode(darray3E);
    ivector1D       addNode(dvecarr3E);
    ivector1D       addNode(MimmoObject* geometry);

    void            setNode(darray3E);
    void            setNode(dvecarr3E);
    void            setNode(MimmoObject* geometry);
    void            setFilter(dmpvector1D * );

    ivector1D       checkDuplicatedNodes(double tol=1.0E-12);
    bool            removeDuplicatedNodes(ivector1D * list=NULL);

    void            setSupportRadius(double suppR_);
    void            setSupportRadiusValue(double suppR_);
    void            setTol(double tol);
    void            setDisplacements(dvecarr3E displ);

    void            setFunction(const MRBFBasisFunction & funct);
    void            setFunction(const bitpit::RBFBasisFunction & funct);

    void            clear();
    void            clearFilter();

    void            execute();
    void            apply();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void            setWeight(dvector2D value);
    void            plotCloud(std::string directory, std::string filename, int counterFile, bool binary, bool deformed);
    virtual void    plotOptionalResults();
    void            swap(MRBF & x) noexcept;
    void            checkFilter();

};

double	heaviside10( double dist );
double	heaviside50( double dist );
double	heaviside100( double dist );
double	heaviside1000( double dist );

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_ ,__MRBF_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_VALUED2, MC_SCALAR, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__MRBF_HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_ ,__MRBF_HPP__)

REGISTER(BaseManipulation, MRBF, "mimmo.MRBF")

};

#endif /* __MRBF_HPP__ */
