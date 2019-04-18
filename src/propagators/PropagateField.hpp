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
#include "StencilFunctions.hpp"

#if MIMMO_ENABLE_MPI
#include "mimmo_parallel.hpp"
#endif

namespace mimmo{
/*!
 * \class PropagateField
 * \ingroup propagators
 * \brief Executable block that provides the computation of a field
 * over a 3D volume mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Class/BaseManipulation Object managing field defined on the boundaries of a 3D volume mesh.
 * It uses MimmoObject informations as input geometry.
 * The key to handle with constraints is an explicit calculation of the solution of a
 * Laplacian problem.
 * The Laplacian solver employs a Finite Volume Cell based scheme.
 * Quality of the mesh during deformation is controlled by introducing an artificial diffusivity D (DUMPING)
 * in the laplacian calculation: div( D*grad(Phi)) = 0. Such diffusivity can be defined alternatively
 * as variable with distance from a prescribed deforming surface (dumping surface),
 * or with cell volumes distribution, modulated with distance from the dumping surface.
 * In the first case, cells more distant from dumping surfaces cells are forced to move more than nearer ones,
 * in the second case, bigger volume cell far from dumping surface are forced to move more.
 * Result field is stored in m_field member and returned as data field through ports.
 * Boundary info passed on POINT are redistributed on CELL boundary Interfaces by means of
 * an interpolation method.
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
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it.
 * - <B>DumpingType</B> : 0- distance control, 1-volume control.
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>DecayFactor</B>  : exponent to modulate dumping function (as power of), if dumping is active
 * - <B>Tolerance</B> : convergence tolerance for laplacian direct solver;
 * - <B>UpdateThres</B> : lower threshold to internally mark cells whose field norm is above its value, for update purposes
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
template<std::size_t NCOMP>
class PropagateField: public mimmo::BaseManipulation {

protected:

    bitpit::PiercedVector<int>   m_isbp;  /**< list of int to mark volume boundary interfaces.
                                          int >0 marks boundary condition for type. For example 1 is a Dirichlet condition, 2 etc...*/
    MimmoPiercedVector<std::array<double, NCOMP> > m_bc_dir; /**< Dirichlet-type condition values interp on boundaries Interfaces of the target volume mesh */
    MimmoPiercedVector<std::array<double, NCOMP> > m_surface_bc_dir; /**< Dirichlet-type condition value of POINTS on boundary surface*/
    MimmoPiercedVector<std::array<double, NCOMP> > m_field;  /**< Resulting Propagated Field on nodes */

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
    int           m_dumpingType;    /**< 0 distance-control, 1-volume control*/

    double        m_thres;          /**< Lower Threshold to internally mark cells whose solution field norm is above its value. For update purposes only */

    std::unique_ptr<bitpit::SystemSolver> m_solver; /**! linear system solver for laplace */

#if MIMMO_ENABLE_MPI
    std::unique_ptr<GhostCommunicator> m_ghostCommunicator; 			/**<Ghost communicator object */
    int m_ghostTag;														/**< Tag of communicator object*/
    std::unique_ptr<MimmoDataBufferStreamer<NCOMP>> m_ghostStreamer;	/**<Data streamer */

    std::unique_ptr<PointGhostCommunicator> m_pointGhostCommunicator; 			/**<Ghost communicator object */
    int m_pointGhostTag;														/**< Tag of communicator object*/
    std::unique_ptr<MimmoPointDataBufferStreamer<NCOMP>> m_pointGhostStreamer;	/**<Data streamer */
#endif

public:

    PropagateField();
    virtual ~PropagateField();
    PropagateField(const PropagateField & other);
    void swap(PropagateField & x) noexcept;

    //get-set methods
    void    buildPorts();
    void    setGeometry(MimmoObject * geometry_);
    void    setDirichletBoundarySurface(MimmoObject*);
    void    setDumpingBoundarySurface(MimmoObject*);

    void    setDumping(bool flag);
    void    setDumpingOuterDistance(double radius);
    void    setDumpingInnerDistance(double plateau);
    void    setDumpingType( int type=0);
    void    setDecayFactor(double decay);
    void    setTolerance(double tol);
    virtual void    setUpdateThreshold(double thres);

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    //cleaners and setters
    virtual void clear();
    virtual void setDefaults();
    virtual bool checkBoundariesCoherence();
    virtual void distributeBCOnBoundaryInterfaces();

    // core resolution functions.
    virtual void computeDumpingFunction();
    virtual void updateDumpingFunction();
    virtual void initializeLaplaceSolver(FVolStencil::MPVDivergence * laplacianStencils, const liimap & maplocals);
    virtual void updateLaplaceSolver(FVolStencil::MPVDivergence * laplacianStencils, const liimap & maplocals);
    virtual void assignBCAndEvaluateRHS(std::size_t comp, bool unused,
                                        FVolStencil::MPVDivergence * borderLaplacianStencil,
                                        FVolStencil::MPVGradient * borderCCGradientStencil,
                                        const liimap & maplocals,
                                        dvector1D & rhs);
    virtual void solveLaplace(const dvector1D &rhs, dvector1D & result);

    virtual void initializeBoundaryInfo();
    virtual void reconstructResults(const dvector2D & results, const liimap & mapglobals,  livector1D * markedcells = nullptr);

#if MIMMO_ENABLE_MPI
    int createGhostCommunicator(bool continuous);
    int createPointGhostCommunicator(bool continuous);
    void communicateGhostData(MimmoPiercedVector<std::array<double, NCOMP> > *data);
    void communicatePointGhostData(MimmoPiercedVector<std::array<double, NCOMP> > *data);

#endif

};

/*!
 * \class PropagateScalarField
 * \ingroup propagators
 * \brief Executable block that provides the computation of a scalar field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Dirichlet boundary conditions are explicitly provided by the User, identifying boundaries
 * through MimmoObject patches and associating to them the value of the field as Dirichlet condition.
 * A natural zero gradient like condition is automatically provided on unbounded borders.
 *
 * The block can perform multistep evaluation to relax field propagation
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
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | setDirichletBoundarySurface           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM3        | setDumpingBoundarySurface             | (MC_SCALAR, MD_MIMMO_) |
    | M_FILTER       | setDirichletConditions                | (MC_MPVECTOR, MD_FLOAT)|

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER         | getPropagatedField                  | (MC_MPVECTOR, MD_FLOAT) |

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
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it.
 * - <B>DumpingType</B> : 0- distance control, 1-volume control.
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>DecayFactor</B>  : exponent to modulate dumping function (as power of), if dumping is active
 * - <B>Tolerance</B> : convergence tolerance for laplacian direct solver;
 *
 * Proper fo the class:
 * - <B>MultiStep</B> : get field solution in a finite number of substeps;
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateScalarField: public mimmo::PropagateField<1> {

public:

    PropagateScalarField();
    PropagateScalarField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateScalarField();
    PropagateScalarField(const PropagateScalarField & other);
    PropagateScalarField & operator=(PropagateScalarField other);
    void swap(PropagateScalarField & x) noexcept;

    void buildPorts();

    dmpvector1D getPropagatedField();

    void    setDirichletConditions(dmpvector1D bc);
    void    setSolverMultiStep(unsigned int sstep);

    //execute
    void        execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    int           m_nstep;               /**< multistep solver */

    //cleaners and setters
    virtual void clear();
    virtual void setDefaults();
    virtual void plotOptionalResults();

private:
    void    setUpdateThreshold(double thres) override{
        BITPIT_UNUSED(thres);
    };
};


/*!
 * \class PropagateVectorField
 * \ingroup propagators
 * \brief Executable block that provides the computation of a 3D array field
 * over a 3D mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given boundary conditions.
 *
 * Prescribed (Dirichlet) field boundary conditions are explicitly provided by the User,
 * identifying boundaries through MimmoObject patches and associating to them the value of each component
 * of the field as Dirichlet condition.
 * Optionally an inpermeability-like/slip condition (zero vector field normal to boundary surface)
 * can be imposed on chosen boundary patches.
 *
 * The block can perform multistep evaluation to relax field propagation
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
    | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>             |
    | M_GEOM          | setGeometry                 | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM2         | setDirichletBoundarySurface |(MC_SCALAR, MD_MIMMO_)   |
    | M_GEOM3         | setDumpingBoundarySurface   | (MC_SCALAR, MD_MIMMO_)  |
    | M_GDISPLS       | setDirichletConditions      | (MC_MPVECARR3, MD_FLOAT)|
    | M_GEOM4         | setSlipBoundarySurface      | (MC_SCALAR, MD_MIMMO_)  |

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS         | getPropagatedField   | (MC_MPVECARR3, MD_FLOAT) |

 *    =========================================================
 *
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
 * - <B>Dumping</B>         : 1-true activate dumping control, 0-false deactivate it.
 * - <B>DumpingType</B> : 0- distance control, 1-volume control.
 * - <B>DumpingInnerDistance</B> : inner limit of dumping function eta, if dumping is active;
 * - <B>DumpingOuterDistance</B> : outer limit of dumping function eta, if dumping is active;
 * - <B>DecayFactor</B>  : exponent to modulate dumping function (as power of), if dumping is active
 * - <B>Tolerance</B> : convergence tolerance for laplacian  direct solver;
 * - <B>UpdateThres</B> : lower threshold to internally mark cells whose field norm is above its value, for update purposes
 *
 * Proper fo the class:
 * - <B>MultiStep</B> : got deformation in a finite number of substep of solution;
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateVectorField: public mimmo::PropagateField<3> {

protected:

    MimmoObject * m_slipsurface;         /**< MimmoObject boundary patch identifying slip conditions */
    int           m_nstep;               /**< multistep solver */
    MimmoPiercedVector<std::array<double, 3> > m_slip_bc_dir; /**< INTERNAL USE ONLY: Slip-type corrector condition values interp on boundaries Interfaces of the target volume mesh */
    MimmoPiercedVector<std::array<double, 3> > m_surface_slip_bc_dir; /**< INTERNAL USE ONLY: Slip-type corrector condition value of POINTS on boundary surface*/

public:

    PropagateVectorField();
    PropagateVectorField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateVectorField();
    PropagateVectorField(const PropagateVectorField & other);
    PropagateVectorField & operator=(PropagateVectorField other);
    void swap(PropagateVectorField & x) noexcept;

    void buildPorts();

    dmpvecarr3E getPropagatedField();

    void    setSlipBoundarySurface(MimmoObject *);

    void    setDirichletConditions(dmpvecarr3E bc);

    void    setSolverMultiStep(unsigned int sstep);

    //execute
    void        execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    //cleaners and setters
    virtual void clear();
    virtual void setDefaults();

    virtual bool checkBoundariesCoherence();
    virtual void distributeSlipBCOnBoundaryInterfaces();

    void apply();

    virtual void plotOptionalResults();

    virtual void assignBCAndEvaluateRHS(std::size_t comp, bool slipCorrect,
                                FVolStencil::MPVDivergence * borderLaplacianStencil,
                                FVolStencil::MPVGradient * borderCCGradientStencil,
                                const liimap & maplocals,
                                dvector1D & rhs);

    virtual void propagateMaskMovingCells(livector1D & celllist);
    virtual void computeSlipBCCorrector(const MimmoPiercedVector<std::array<double,3> > & guessSolutionOnPoint);
    virtual void subdivideBC();
    virtual void restoreBC();
    void restoreGeometry(bitpit::PiercedVector<bitpit::Vertex> & vertices);
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

#include "PropagateField.tpp"

#endif /* __PROPAGATEFIELD_HPP__ */
