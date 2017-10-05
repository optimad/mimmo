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
 * \ingroup core
 * \brief Executable block that provides the computation of a field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Class/BaseManipulation Object managing field defined on the boundaries of a 3D volume mesh.
 * It uses MimmoObject informations as input geometry.
 * The key to handle with constraints is an explicit calculation of the solution of a
 * Laplacian problem. Only Dirichlet boundary condition are allowed.
 * Boundary conditions are explicitly provided by the User, identifying boundaries
 * through MimmoObject patches and associating to them a field used as Dirichlet conditions.
 * In Laplacian solver the weights used to compute the stencil are of the form:
 * w(x) = (1/d(x))^n, where d(x) is the distance of the neighbour point to the
 * center point of the stencil and n is the weight factor, that can be tuned by the user.
 * A dumping function can be imposed during the Laplacian solution. The dumping function
 * is an exponential law of the distance from the boundaries with Dirichlet condition
 * different from zero. The function is of the form: f(x) = (r/d(x))^p where
 * d(x) is the distance of x from the boundaries with condition !=0,
 * r is the dumping support radius (it can be set by the user) and
 * p is the dumping factor (tunable by the user). A dumping factor = 0.0 turns off
 * the dumping function.
 * Result field is stored in m_field member and returned as data field through ports.
 *
 * Ports available in PropagateField Class :
 *
 *    =========================================================
 *
    |                   Port Input       |||                               |
    |-------|----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setBoundarySurface                    | (MC_SCALAR, MD_MIMMO_) |

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
 * - <B>WeightFactor</B> : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>DumpingFactor</B> : dumping exponential factor for weights computing ;
 * - <B>DumpingRadius</B> : support radius of dumping function;
 * - <B>Solver</B>  : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
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

    bmpvector1D   m_isbp;            /**< Boolean flags to mark boundary points. */
    int           m_np;              /**< Number of points of the cloud. */
    int           m_nbp;             /**< Number of boundary points of the cloud. */
    limpvector2D  m_conn;            /**< Neighbor points connectivity. */
    double        m_gamma;           /**< Coefficient used to weight the points with distance in stencil computing. */
    dmpvector2D   m_weights;         /**< Weights used in stencil computing. */
    bool          m_laplace;         /**<Set true for laplace solver, false for smoothing solver (Laplacian solver not implemented in this version).*/
    int           m_sstep;           /**<Number of smoothing steps [default =10].*/
    bool          m_convergence;     /**<Convergence flag. If true the laplacian smoothing is solved until convergences is reached [default false].  */
    double        m_tol;             /**<Convergence tolerance. [default tol = 1.0e-05 on maximum differences (relative to maximum value on boundary conditions) of solution between two iterations].*/
    MimmoObject*  m_bsurface;        /**<Pointer to MimmoObject with boundaries vertices.*/
    bool          m_bPointsToCompute;/**<Auxiliary variable to compute bPoints at right time.*/
    double        m_dumpingFactor;   /**<Dumping exponential factor for weights computing (0.0 = dumping inactive).*/
    dmpvector1D   m_dumping;         /**<Dumping field used for weights computing.*/
    double        m_radius;          /**<Support radius of dumping function. At distance = m_radius from boundary with bc != 0
                                         the stencil during the laplacian computing is the original one.*/

    std::unique_ptr<mimmo::SystemSolver> m_solver;


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
    void    setBoundarySurface(MimmoObject*);
    void    setWeightConstant(double gamma);
    void    setSmoothingSteps(int ns);

    void    setSolver(bool solveLaplacian);
    void    setDumpingFactor(double dump);
    void    setDumpingRadius(double radius);
    void    setConvergence(bool convergence);
    void    setTolerance(double tol);

protected:
    //cleaners and setters
    void clear();
    void setDefaults();

    //execute
    void        execute();
    void        computeConnectivity();
    void        computeWeights();
    /*!Compute dumping function.*/
    virtual void computeDumpingFunction() = 0;
    /*!Propagate field with a "Laplacian smoothing" iterative solver*/
    virtual void solveSmoothing(int nstep) = 0;
    /*!Propagate field solving directly the Laplace equation */ 
    virtual void solveLaplace() = 0;

};

/*!
 * \class PropagateScalarField
 * \brief Executable block that provides the computation of a scalar field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Class/BaseManipulation Object specialization of class PropagateField
 * for the propagation in a volume mesh of a scalar field.
 *
 * Ports available in PropagateScalarField Class :
 *
 *    =========================================================
 *
    |                   Port Input       |||                               |
    |-------|----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER       | setBoundaryConditions        | (MC_MPVECTOR, MD_FLOAT) |


     |             Port Output     |||                                      |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER         | getDeformation    | (MC_MPVECTOR, MD_FLOAT)              |

  Inherited from PropagateField Class :

    |                   Port Input       |||                               |
    |-------|----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setBoundarySurface                    | (MC_SCALAR, MD_MIMMO_) |


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
 * Inherited from PropagateField:
 * - <B>Solver</B>  : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
 * - <B>WeightConstant</B> : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>SmoothingSteps</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>DumpingFactor</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>DumpingRadius</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>Convergence</B> : number of steps the Smoother solver need to perform (1 default);
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateScalarField: public mimmo::PropagateField {

private:

    dmpvector1D   m_bc;              /**< Boundary conditions. */
    dmpvector1D   m_field;           /**< Resulting field.*/

public:

    PropagateScalarField();
    PropagateScalarField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateScalarField();
    PropagateScalarField(const PropagateScalarField & other);
    PropagateScalarField & operator=(PropagateScalarField other);
    void swap(PropagateScalarField & x) noexcept;

    void buildPorts();

    dmpvector1D getField();

    void    setBoundaryConditions(dmpvector1D bc);

private:
    //cleaners and setters
    void clear();
    void setDefaults();

    //execute
    void computeDumpingFunction();
    void solveSmoothing(int nstep);
    void solveLaplace();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    virtual void plotOptionalResults();

};


/*!
 * \class PropagateVectorField
 * \brief Executable block that provides the computation of a 3D array field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Class/BaseManipulation Object specialization of class PropagateField
 * for the propagation in a volume mesh of a 3D array field.
 *
 * Ports available in PropagateVectorField Class :
 *
 *    =========================================================
 *
    |                   Port Input       |||                               |
    |-------|----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS       | setBoundaryConditions        | (MC_MPVECARR3, MD_FLOAT) |


     |             Port Output     |||                                      |
     |-------|----------------|---------------------|--------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS         | getDeformation    | (MC_MPVECARR3, MD_FLOAT)              |

  Inherited from PropagateField Class :

    |                   Port Input       |||                               |
    |-------|----------------|----------------------|-------------------|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setBoundarySurface                    | (MC_SCALAR, MD_MIMMO_) |


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
 * Inherited from PropagateField:
 * - <B>Solver</B>  : 1-true use direct Laplacian solver, 0-false use iterative Smoother (default 0);
 * - <B>WeightConstant</B> : coefficient used to get weights of the stencil points in function of distance (1.0 default);
 * - <B>SmoothingSteps</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>DumpingFactor</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>DumpingRadius</B> : number of steps the Smoother solver need to perform (1 default);
 * - <B>Convergence</B> : number of steps the Smoother solver need to perform (1 default);
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateVectorField: public mimmo::PropagateField {

private:

    dmpvecarr3E   m_bc;              /**< Boundary conditions. */
    dmpvecarr3E   m_field;           /**< Resulting field.*/

public:

    PropagateVectorField();
    PropagateVectorField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateVectorField();
    PropagateVectorField(const PropagateVectorField & other);
    PropagateVectorField & operator=(PropagateVectorField other);
    void swap(PropagateVectorField & x) noexcept;

    void buildPorts();

    dmpvecarr3E getField();

    void    setBoundaryConditions(dmpvecarr3E bc);

private:
    //cleaners and setters
    void clear();
    void setDefaults();

    //execute
    void computeDumpingFunction();
    void solveSmoothing(int nstep);
    void solveLaplace();

    void apply();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

    virtual void plotOptionalResults();

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_FILTER, MC_MPVECTOR, MD_FLOAT,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GDISPLS, MC_MPVECARR3, MD_FLOAT,__PROPAGATEFIELD_HPP__)

REGISTER(BaseManipulation, PropagateScalarField, "mimmo.PropagateScalarField")
REGISTER(BaseManipulation, PropagateVectorField, "mimmo.PropagateVectorField")

};

#endif /* __PROPAGATEFIELD_HPP__ */
