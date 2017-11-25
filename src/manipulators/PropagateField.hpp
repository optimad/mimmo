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
#ifndef __PROPAGATEFIELD_HPP__
#define __PROPAGATEFIELD_HPP__

#include "BaseManipulation.hpp"
#include "system.hpp"

namespace mimmo{
/*!
 * \class PropagateField
 * \ingroup manipulators
 * \brief Executable block that provides the computation of a field
 * over a 3D volume mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Class/BaseManipulation Object managing field defined on the boundaries of a 3D volume mesh.
 * It uses MimmoObject informations as input geometry.
 * The key to handle with constraints is an explicit calculation of the solution of a
 * Laplacian problem. 
 * In Laplacian solver the weights used to compute the stencil are of the form:
 * w(x) = (1/d(x))^n, where d(x) is the distance of the neighbour point to the
 * center point of the stencil and n is the weight factor, that can be tuned by the User.
 * Quality of the mesh during deformation is controlled by introducing an artificial diffusivity D
 * in the laplacian calculation: div( D*grad(Phi)) = 0. Such diffusivity varies from cell to cell, according to 
 * their typical size and aspect ratio. The smaller and high-AR the cell, the lesser the boundary deformation is 
 * propagated in that cell, according to the formula:
 * 
 *          D_cell = 1.0/(1 + T_cell*AR_cell) ^ eta.
 * 
 * Where D_cell is the value of diffusivity in the cell, T_cell is the cell volume factor defined as:
 *         
 *         T_cell = (V_max - V_min) / V_cell,
 * 
 * where V_max and V_min are the maximum and minimum volume cell in the whole mesh, and V_cell 
 * is the current volume of the target cell. AR_cell is the Aspect Ratio of each cell, while eta is 
 * a dumping factor that can be imposed to control the effects of the artificial diffusivity D in 
 * the bulk mesh. In particular, eta is a function of the distance of a point x from mesh 
 * boundaries selected by the User (typically those boundaries which are deformed). 
 * The function is of the form: 
 * 
 *          eta(x) = r/p              d(x)<= p,
 *          eta(x) = (r/d(x))         p < d(x) < r ,
 *          eta(x) = 0.0              d(x) => r
 * 
 * where p and r are limit distances choosen by the User.
 * In general the effects of the artificial variable diffusivity are maximum in a the region with d(x) < p, 
 * and completely shut off at d(x)> r.
 * Result field is stored in m_field member and returned as data field through ports.
 *
 * Ports available in PropagateField Class :
 *
 *    =========================================================
 *
    |Port Input|||
    ||||
    | <B>PortType</B>|<B>variable/function</B>      |<B>DataType</B>|
    | M_GEOM         | setGeometry                  | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setDirichletBoundarySurface  | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM3        | setDumpingBoundarySurface    | (MC_SCALAR, MD_MIMMO_) |
 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.PropagateField</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : plot optional results in execution;
 * - <B>OutputPlot</B> : path to store optional results.
 *
 * Proper of the class:
 * - <B>Solver</B>        : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
 * - <B>WeightConstant</B>  : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it. 
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>SmoothingSteps</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>Convergence</B> : convergence flag for smoothing solver;
 * - <B>Tolerance</B> : convergence tolerance for laplacian smoothing and direct solver;
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateField: public mimmo::BaseManipulation {

protected:

    btpvector1D   m_isbp;            /**< list of std::pair< bool, uint8_t> to mark points. bool true identifies boundary points, false internal points.
                                          uint8_t marks as 1 Dirichlet condition on boundary points, 2 other boundary conditions, 0 no condition.*/
    int           m_np;              /**< Number of points of the cloud. */
    int           m_nbp;             /**< Number of boundary points of the cloud. */
    limpvector2D  m_conn;            /**< Neighbor points connectivity. */
    double        m_gamma;           /**< Coefficient used to weight the points with distance in stencil computing. */
    dmpvector2D   m_weights;         /**< Weights used in stencil computing. */
    bool          m_laplace;         /**<Set true for laplace solver, false for smoothing solver (Laplacian solver not implemented in this version).*/
    int           m_sstep;           /**<Number of smoothing steps [default =10].*/
    bool          m_convergence;     /**<Convergence flag. If true the laplacian smoothing is solved until convergences is reached [default false].  */
    double        m_tol;             /**<Convergence tolerance. [default tol = 1.0e-05 on maximum differences (relative to maximum value on boundary conditions) of solution between two iterations].*/
    MimmoObject*  m_bsurface;        /**<Pointer to MimmoObject boundary patch identifying Dirichlet boundaries.*/
    MimmoObject*  m_dsurface;        /**<Pointer to MimmoObject with selected boundaries for dumping calculation.*/
    double        m_decayFactor;     /**<Dumping exponential factor for hybrid diffusivity region.*/
    dmpvector1D   m_dumping;         /**<Dumping field used for weights computing.*/
    double        m_radius;          /**<Outer limit distance of dumping function. At distance >= m_radius from boundary with bc != 0
                                         the stencil during the laplacian computing is the original one.*/
    double        m_plateau;        /**<Inner limit distance of dumping function. At distance <= m_plateau from boundary with bc != 0
                                         the stencil during the laplacian computing account of the maximum artificial diffusivity.*/
    bool          m_dumpingActive;  /**< true the dumping control is active, false otherwise.*/

    std::unique_ptr<mimmo::SystemSolver> m_solver; /**! linear system solver for laplace */


public:

    PropagateField();
    virtual ~PropagateField();
    PropagateField(const PropagateField & other);
    void swap(PropagateField & x) noexcept;

    void buildPorts();

    //get-set methods
    int         getNPoints();
    int         getNBoundaryPoints();

    void    setGeometry(MimmoObject * geometry_);
    void    setDirichletBoundarySurface(MimmoObject*);
    void    setDumpingBoundarySurface(MimmoObject*);
    
    void    setWeightConstant(double gamma);
    void    setSmoothingSteps(int ns);

    void    setSolver(bool solveLaplacian);
    void    setDumping(bool flag);
    void    setDumpingOuterDistance(double radius);
    void    setDumpingInnerDistance(double plateau);
    void    setDecayFactor(double decay);
    void    setConvergence(bool convergence);
    void    setTolerance(double tol);

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
    
protected:
    //cleaners and setters
    void clear();
    void setDefaults();

    //execute
    void        execute();
    void        computeConnectivity();
    void        computeWeights();

    virtual void computeDumpingFunction();
    /*!Propagate field with a "Laplacian smoothing" iterative solver*/
    virtual void solveSmoothing(int nstep) = 0;
    /*!Propagate field solving directly the Laplace equation */ 
    virtual void solveLaplace() = 0;

private:
    virtual bool checkBoundariesCoherence() = 0;
};

/*!
 * \class PropagateScalarField
 * \ingroup manipulators
 * \brief Executable block that provides the computation of a scalar field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet/Neumann boundary conditions.
 * 
 * Neumann and Dirichlet boundary conditions are explicitly provided by the User, identifying boundaries
 * through MimmoObject patches and associating to them the value of the field as Dirichlet condition, 
 * or the value of boundary normal flux as Neumann condition.
 * Dirichlet info are mandatory, Neumann ones are optional.
 *
 * Class/BaseManipulation Object specialization of class PropagateField
 * for the propagation in a volume mesh of a scalar field.
 *
 * Ports available in PropagateScalarField Class :
 *
 *    =========================================================
 *
    |Port Input|||
    ||||
    | <B>PortType</B>| <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER       | setDirichletConditions    | (MC_MPVECTOR, MD_FLOAT) |
    | M_GEOM4        | setNeumannBoundarySurface | (MC_SCALAR, MD_MIMMO_)  |

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER         | getPropagatedField    | (MC_MPVECTOR, MD_FLOAT)              |

  Inherited from PropagateField Class :

    |Port Input|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setDirichletBoundarySurface           | (MC_SCALAR, MD_MIMMO_) |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.PropagateScalarField</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B> : plot optional results in execution;
 * - <B>OutputPlot</B> : path to store optional results.
 *
 * Inherited from PropagateField:
 * - <B>Solver</B>        : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
 * - <B>WeightConstant</B>  : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it. 
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>SmoothingSteps</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>Convergence</B> : convergence flag for smoothing solver;
 * - <B>Tolerance</B> : convergence tolerance for laplacian smoothing and direct solver;
 
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateScalarField: public mimmo::PropagateField {

private:

    dmpvector1D   m_bc_dir;              /**< Dirichlet Boundary conditions. */
    dmpvector1D   m_field;               /**< Resulting field.*/
    MimmoObject * m_neu_surface;         /**< MimmoObject boundary patch identifying Neumann conditions */
public:

    PropagateScalarField();
    PropagateScalarField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateScalarField();
    PropagateScalarField(const PropagateScalarField & other);
    PropagateScalarField & operator=(PropagateScalarField other);
    void swap(PropagateScalarField & x) noexcept;

    void buildPorts();

    dmpvector1D getPropagatedField();
    
    void    setNeumannBoundarySurface(MimmoObject * );
    
    void    setDirichletConditions(dmpvector1D bc);

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

private:
    //cleaners and setters
    void clear();
    void setDefaults();

    bool checkBoundariesCoherence();
    void correctStencilForNeumann();
    //execute
    void solveSmoothing(int nstep);
    void solveLaplace();

    virtual void plotOptionalResults();

};


/*!
 * \class PropagateVectorField
 * \ingroup manipulators
 * \brief Executable block that provides the computation of a 3D array field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given boundary conditions.
 * 
 * Neumann type or Prescribed (Dirichlet) field boundary conditions are explicitly provided by the User, 
 * identifying boundaries through MimmoObject patches and associating to them the value of each component 
 * of the field as Dirichlet condition, or Neumann condition for each component of the field at the boundary.
 * 
 * Dirichlet info are mandatory, Neumann ones are optional.
 * 
 * Class/BaseManipulation Object specialization of class PropagateField
 * for the propagation in a volume mesh of a 3D array field.
 *
 * Ports available in PropagateVectorField Class :
 *
 *    =========================================================
 *
    | Port Input|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS      | setDirichletConditions    | (MC_MPVECARR3, MD_FLOAT)|
    | M_GEOM4        | setNeumannBoundarySurface | (MC_SCALAR, MD_MIMMO_)  |

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS         | getPropagatedField   | (MC_MPVECARR3, MD_FLOAT)              |

  Inherited from PropagateField Class :

    |Port Input|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setBoundarySurface                    | (MC_SCALAR, MD_MIMMO_) |


 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B> : name of the class as <tt>mimmo.PropagateVectorField</tt>;
 * - <B>Priority</B>  : uint marking priority in multi-chain execution;
 * - <B>Apply</B> : if set to 1, apply propagated field to target geometry as a deformation field;
 * - <B>PlotInExecution</B> : plot optional results in execution;
 * - <B>OutputPlot</B> : path to store optional results.
 *
 * Inherited from PropagateField:
 * - <B>Solver</B>        : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
 * - <B>WeightConstant</B>  : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it. 
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>SmoothingSteps</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>Convergence</B> : convergence flag for smoothing solver;
 * - <B>Tolerance</B> : convergence tolerance for laplacian smoothing and direct solver;
 
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateVectorField: public mimmo::PropagateField {

private:

    dmpvecarr3E   m_bc_dir;              /**< Dirichlet Boundary conditions. */
    dmpvecarr3E   m_field;           /**< Resulting field.*/
    MimmoObject * m_neu_surface;         /**< MimmoObject boundary patch identifying Neumann conditions */
public:

public:

    PropagateVectorField();
    PropagateVectorField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateVectorField();
    PropagateVectorField(const PropagateVectorField & other);
    PropagateVectorField & operator=(PropagateVectorField other);
    void swap(PropagateVectorField & x) noexcept;

    void buildPorts();

    dmpvecarr3E getPropagatedField();
    
    void    setNeumannBoundarySurface(MimmoObject *);
    
    void    setDirichletConditions(dmpvecarr3E bc);

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

private:
    //cleaners and setters
    void clear();
    void setDefaults();

    bool checkBoundariesCoherence();
    void correctStencilForNeumann();
    
    //execute
    void solveSmoothing(int nstep);
    void solveLaplace();

    void apply();

    virtual void plotOptionalResults();

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM3, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM4, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_FILTER, MC_MPVECTOR, MD_FLOAT,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GDISPLS, MC_MPVECARR3, MD_FLOAT,__PROPAGATEFIELD_HPP__)

REGISTER(BaseManipulation, PropagateScalarField, "mimmo.PropagateScalarField")
REGISTER(BaseManipulation, PropagateVectorField, "mimmo.PropagateVectorField")

};

#endif /* __PROPAGATEFIELD_HPP__ */
