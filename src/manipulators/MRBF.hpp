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
    DSIGMOID    = 105,  /**< Non compact sharp sigmoid derivative */
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
   given externally or stored in a MimmoObject container.
   Displacements (DOFs) for RBF nodes are provided externally.
   Default solver in execution is MRBFSol::NONE for direct parameterization.
   Use MRBFSol::GREEDY or MRBFSol::WHOLE to activate interpolation features.
 * See bitpit::RBF docs for further information.

    Support radii of RBF Nodes can be set in 3 different ways:
    - setting x as Local support radius : the effective support radius will be
      calculated as x * bbox_diag, where the last is the diagonal of the RBF set
      of nodes bounding box. This value will be equal for all RBF nodes.(setSupportRadiusLocal)
 * - setting x as Real support radius: the value of x will be used as effective support
     radius, equal for all RBF nodes (setSupportRadiusReal)
   - setting variable Support Radii: provide a list of support radius values, one for
     each RBF node in the set. This option is available only in MRBFSol::NONE mode
    (setVariableSupportRadii)

     In case all three methods are connected through ports, priority list will be:
     setVariableSupportRadii, setSupportRadiusReal, setSupportRadiusLocal
 * \n
 * Ports available in MRBF Class :
 *
 *    =========================================================

     |Port Input | | |
     |-|-|-|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_COORDS   | setNode                | (MC_VECARR3, MD_FLOAT)      |
     | M_DISPLS   | setDisplacements       | (MC_VECARR3, MD_FLOAT)      |
     | M_FILTER   | setFilter              | (MC_SCALAR, MD_MPVECFLOAT_) |
     | M_DATAFIELD| setVariableSupportRadii| (MC_VECTOR, MD_FLOAT)       |
     | M_GEOM     | setGeometry            | (MC_SCALAR, MD_MIMMO_)      |

     |Port Output | | |
     |-|-|-|
     | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>|
     | M_GDISPLS       | getDisplacements  | (MC_SCALAR, MD_MPVECARR3FLOAT_)             |
     | M_GEOM          | getGeometry       | (MC_SCALAR, MD_MIMMO_) |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.MRBF</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>Apply</B>: boolean 0/1 activate apply deformation result on target geometry directly in execution;
 * - <B>PlotInExecution</B> : plot optional results in execution;
 * - <B>OutputPlot</B> : path to store optional results.

 * Proper of the class:
 * - <B>Mode</B>: 0/1/2 mode of usage of the class see MRBFSol enum;
 * - <B>SupportRadiusLocal</B>: local homogeneous radius of RBF function for each nodes,
                                expressed as ratio of local geometry bounding box;
                                see setSupportRadiusLocal method documentation.
 * - <B>SupportRadiusReal</B>: homogeneous real radius of RBF function common to each RBF node;
                                see setSupportRadiusReal method documentation.
 * - <B>RBFShape</B>: shape of RBF function see MRBFBasisFunction and bitpit::RBFBasisFunction enums;
 * - <B>Tolerance</B>: greedy engine tolerance (meaningful for Mode 2 only);
 *
    if set, SupportRadiusReal parameter bypass SupportRadiusLocal one.

 * Geometry, filter field, RBF nodes and displacements have to be mandatorily passed through port.
 * Please not if class is in MRBFSol::NONE mode and a port to setVariableSupportRadii is linked,
   SupportRadiusLocal/Real parameters are bypassed.
 *
 */
//TODO study how to manipulate supportRadius of RBF to define a local/global smoothing of RBF
class MRBF: public BaseManipulation, public bitpit::RBF {

protected:
    double       m_tol;          /**< Tolerance for greedy algorithm.*/
    MRBFSol      m_solver;       /**<Type of solver specified for the class as default in execution*/
    dmpvector1D  m_filter;       /**<Filter field for displacements modulation */
    bool         m_bfilter;      /**<boolean to recognize if a filter field is applied */
    double       m_supportRadiusValue; /**<local bounding box binded homogeneous support radius */
    bool         m_srIsReal;  /**<True if homogeneous support radius is defined as absolute value, false if is ratio of bounding box diagonal.*/
    dvector1D    m_supportRadii; /**< list of variable supportRadii for each RBF node.*/
    int          m_functype;     /**< Function type handler. If -1 refer to RBF getFunctionType method */
    dmpvecarr3E  m_displ;        /**<Resulting displacements of geometry vertex.*/

    dvector1D    m_effectiveSR; /**< INTERNAL USE list of support radii effectively used for each RBF */

public:
    MRBF(MRBFSol mode = MRBFSol::NONE);
    MRBF(const bitpit::Config::Section & rootXML);

    virtual ~MRBF();

    //copy operators/constructors
    MRBF(const MRBF & other);
    MRBF& operator=(MRBF other);

    void buildPorts();

    dvecarr3E*      getNodes();

    MRBFSol         getMode();
    dmpvector1D*    getFilter();

    bool            isSupportRadiusReal();
    bool            isVariableSupportRadiusSet();
    dvector1D &     getEffectivelyUsedSupportRadii();

    int             getFunctionType();
    dmpvecarr3E*    getDisplacements();

    int             addNode(darray3E);
    ivector1D       addNode(dvecarr3E);

    void            setNode(darray3E);
    void            setNode(dvecarr3E);
    void            setNode(MimmoObject* geometry);
    void            setFilter(dmpvector1D * );

    ivector1D       checkDuplicatedNodes(double tol=1.0E-12);
    bool            removeDuplicatedNodes(ivector1D * list=NULL);

    void            setSupportRadiusLocal(double suppR_);
    void            setSupportRadiusReal(double suppR_);
    void            setVariableSupportRadii(dvector1D sradii);

BITPIT_DEPRECATED(
    void            setSupportRadiusValue(double suppR_));

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
    void            swap(MRBF & x) noexcept;
    void            checkFilter();
    void            setWeight(dvector2D value);
    void            plotCloud(std::string directory, std::string filename, int counterFile, bool binary, bool deformed);
    virtual void    plotOptionalResults();

    void            computeEffectiveSupportRadiusList();

    livector1D searchSphereMatches(bitpit::PatchSkdTree & tree, const std::array<double,3>& sphere_center, double sphere_radius);
    bool       intersectSphereAABBox(const darray3E & center, double radius, const darray3E &bMin, const darray3E &bMax);

    //reimplemented from RBFKernel
    void            setMode(MRBFSol solver);
    std::array<double,3> evalRBF(const std::array<double,3> & val);

private:

    // redeclared to be private in the current class:
    // they behave exactly as in the base class but they are not accessible in MRBF.

    void  setSupportRadius(double val){bitpit::RBFKernel::setSupportRadius(val);};
    double getSupportRadius(){return bitpit::RBFKernel::getSupportRadius();};
    void  setMode(bitpit::RBFMode mode){bitpit::RBFKernel::setMode(mode);};

    void  setDataToNode (int n, const std::vector<double> &data){bitpit::RBFKernel::setDataToNode(n,data);};
    void  setDataToAllNodes(int n , const std::vector<double> & data){bitpit::RBFKernel::setDataToAllNodes(n,data);};

    int   addData(){return bitpit::RBFKernel::addData();};
    int   addData(const std::vector<double> & data){return bitpit::RBFKernel::addData(data);};
    bool  removeData(int n){return bitpit::RBFKernel::removeData(n);};
    bool  removeData(std::vector<int> &nn){return bitpit::RBFKernel::removeData(nn);};
    void  removeAllData(){bitpit::RBFKernel::removeAllData();};

    std::vector<double>     evalRBF(int jnode){return bitpit::RBFKernel::evalRBF(jnode);};
    
    int                     solve(){return bitpit::RBFKernel::solve();};
    int                     greedy(double tol){return bitpit::RBFKernel::greedy(tol);};

};

double	heaviside10( double dist );
double	heaviside50( double dist );
double	heaviside100( double dist );
double  heaviside1000( double dist );
double  dsigmoid( double dist );

REGISTER_PORT(M_COORDS, MC_VECARR3, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_DISPLS, MC_VECARR3, MD_FLOAT ,__MRBF_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_ ,__MRBF_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__MRBF_HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_ ,__MRBF_HPP__)
REGISTER_PORT(M_DATAFIELD, MC_VECTOR, MD_FLOAT, __MRBF_HPP_)

REGISTER(BaseManipulation, MRBF, "mimmo.MRBF")

};

#endif /* __MRBF_HPP__ */
