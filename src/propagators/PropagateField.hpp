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
 * over a 3D volume or surface mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given boundary conditions.
 *
 * Class/BaseManipulation Object managing field defined on the boundaries of a
   3D volume or surface mesh. It uses MimmoObject informations as input geometry.
 * The key to handle with constraints is an explicit calculation of the solution of a
 * Laplacian problem.
 * The Laplacian solver employs a node-based scheme (Discrete Graph Laplacian operator).

   To improve quality of the boundary information propagation inside the mesh two
   approaches that directly influence the solution of Laplacian
   \f$\nabla\cdot( \nabla\Phi) = 0\f$ are proposed:

   - Narrow Band Control (NBC) : in a the neighbourhood of some prescribed boundary surfaces
    Laplacian stencils can be altered by using a Weighted Graph Lapalcian operator. The objectove is to
    increase the diffusion length of the boundary conditions inside the bulk mesh.
    This approach requires the User to specify the target boundary surfaces,
    a target distance from them (width of the narrow band) and a relaxation parameter to tune
    the control.

   - Artificial Diffusivity Control (DAMPING function): laplacian solution is manipulated introducing an
    artificial diffusivity D function inside the governing equation \f$\nabla\cdot( D\nabla\Phi) = 0\f$.
    Such diffusivity can be defined alternatively as variable with distance from
    prescribed deforming surfaces (damping surfaces), or with cell volumes distribution,
    even modulated with distance from the damping surface. In the first case, cells more distant
    from damping surfaces allow greater deformations than nearer ones, that have a more stiff behavior.
    In the second case, bigger volume cell far from dumping surface are less stiff and more deformable.
    The effect of artificial diffusivity is confined into a bulk mesh zone starting from the
    reference damping surfaces up to a certain distance from them. This distance
    can be tuned providing the Damping Outer Distance parameter.
    Damping Inner Distance parameter represents an intermediate distance between boundaries and Damping Outer Distance
    so that the damping function decays from inner to outer distance and is set to constant from
    damping surfaces up to inner distance.
   <B>WARNING</B>: Damping Function Control is an EXPERIMENTAL feature, subject
                   of still ongoing investigations. Sometimes it can lead to
                   unpredicatable results.


    <B>Note.</B> Currently, NBC and artificial diffusivity are available only for volume bulk mesh.

 * Result field is stored in m_field member and returned as data field through ports.
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
 * - <B>Damping</B>              : 1-true activate damping control, 0-false deactivate it.
 * - <B>DampingType</B>          : 0- distance control, 1-volume control.
 * - <B>DampingInnerDistance</B> : inner limit of damping function, if damping is active;
 * - <B>DampingOuterDistance</B> : outer limit of damping function, if damping is active;
 * - <B>DecayFactor</B>          : exponent to modulate damping function (as power of),
                                   if damping is active
   - <B>NarrowBand</B>           : 1-true activate narrow band control, 0-false deactivate it.
   - <B>NarrowBandWidth</B>      : size of narrow band (> 0.0), if NarrowBand is
                                   active.
   - <B>NarrowBandRelax</B>      : NBC relaxation parameter within [0,1], where
                                   1 means no relaxation, 0 full relaxation. Meaningful
                                   only if NarrowBand is active.
 * - <B>Tolerance</B>            : convergence tolerance for laplacian solver.
 * - <B>UpdateThres</B>          : lower threshold to internally mark cells whose
                                   field norm is above its value, for update purposes
 * - <B>Print</B>                : print solver debug information, Active only
                                   if MIMMO_ENABLE_MPI is enabled in compilation.
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
template<std::size_t NCOMP>
class PropagateField: public mimmo::BaseManipulation {

protected:
    // General members
    double        m_thres;          /**< Lower Threshold to internally mark cells whose solution field norm is above its value. For update purposes only */
    double        m_tol;            /**< Convergence tolerance. [default tol = 1.0e-12].*/
    bool          m_print;          /**< If true residuals and other info are print during system solving.*/

    std::unique_ptr<bitpit::SystemSolver> m_solver;             /**< linear system solver for Laplace */
    MimmoPiercedVector<std::array<double, NCOMP> > m_field;     /**< Resulting Propagated Field on bulk nodes */

    // Dirichlet surfaces and bcs
    std::unordered_set<MimmoSharedPointer<MimmoObject> >  m_dirichletPatches;              /**<list of MimmoObject boundary patches pointers identifying Dirichlet boundaries.*/
    MimmoPiercedVector<std::array<double, NCOMP> > m_bc_dir;                                /**< Dirichlet-type condition values of POINTS on target volume mesh */
    std::unordered_set<MimmoPiercedVector<std::array<double, NCOMP> >* > m_dirichletBcs;    /**< list of Dirichlet-type conditions on boundary patch POINTS*/

    // Damping surfaces and parameters set
    bool          m_dampingActive;  /**< true the damping control is active, false otherwise.*/
    int           m_dampingType;    /**< 0 distance-control, 1-volume control*/
    double        m_decayFactor;    /**<Damping exponential factor for hybrid diffusivity region.*/
    double        m_radius;         /**<Outer limit distance of damping function. At distance >= m_radius from boundary with bc != 0
                                        the stencil during the laplacian computing is the original one.*/
    double        m_plateau;        /**<Inner limit distance of damping function. At distance <= m_plateau from boundary with bc != 0
                                        the stencil during the laplacian computing account of the maximum artificial diffusivity.*/
    dmpvector1D   m_damping;        /**<Damping field used for weights computing.*/

    std::unordered_set<MimmoSharedPointer<MimmoObject> >  m_dampingSurfaces;    /**<list of MimmoObject boundary patches pointers to identify surface for damping calculation.*/
    MimmoSharedPointer<MimmoObject> m_dampingUniSurface;                        /**<INTERNAL use. Final damping reference surface.*/

    // Narrow Band surfaces and parameters set
    bool          m_bandActive;  /**< true the Narrow Band Control is active, false otherwise.*/
    double        m_bandwidth;   /**< width of the narrow band region.*/
    double        m_bandrelax;   /**< Narrow band relaxation param [0,1]. 1 no relaxing occurs, 0 full relaxation is performed */

    bitpit::PiercedVector<double> m_banddistances; /**< INTERNAL use, list of distances for vertex belonging to narrow band */
    std::unordered_set<MimmoSharedPointer<MimmoObject> >  m_bandSurfaces;   /**<list of MimmoObject boundary patches pointers to identify target baundaries for Narrow Band definition.*/
    MimmoSharedPointer<MimmoObject> m_bandUniSurface; /**< INTERNAL use. Final narrow band reference surface.*/

public:

    PropagateField();
    virtual ~PropagateField();
    PropagateField(const PropagateField & other);
    void swap(PropagateField & x) noexcept;

    void    buildPorts();
    void    setTolerance(double tol);
    virtual void    setUpdateThreshold(double thres);
    void	setPrint(bool print = true);

    void    setGeometry(MimmoSharedPointer<MimmoObject> geometry_);
    BITPIT_DEPRECATED(void    addDirichletBoundarySurface(MimmoSharedPointer<MimmoObject>));
    void    addDirichletBoundaryPatch(MimmoSharedPointer<MimmoObject>);

    void    setNarrowBand(bool flag);
    void    addNarrowBandBoundarySurface(MimmoSharedPointer<MimmoObject>);
    void    setNarrowBandWidth(double width);
    void    setNarrowBandRelaxation(double relax);

    void    setDamping(bool flag);
    void    setDampingType( int type=0);
    void    setDampingDecayFactor(double decay);
    void    addDampingBoundarySurface(MimmoSharedPointer<MimmoObject>);
    void    setDampingInnerDistance(double plateau);
    void    setDampingOuterDistance(double radius);

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    //cleaners and setters
    virtual void clear();
    virtual void setDefaults();
    virtual bool checkBoundariesCoherence();
    virtual void distributeBCOnBoundaryPoints();

    // core resolution functions.
    void initializeUniqueSurface(const std::unordered_set<MimmoSharedPointer<MimmoObject> > & listSurf, MimmoSharedPointer<MimmoObject> & uniSurf);

    //Damping
    virtual void initializeDampingFunction();
    virtual void updateDampingFunction();
    void dampingCellToPoint(MimmoPiercedVector<double> & damping);

    //NarrowBand
    virtual void updateNarrowBand();
    virtual void modifyStencilsForNarrowBand(GraphLaplStencil::MPVStencilUPtr &laplaceStencils );

    // Laplace Solver
    virtual void initializeLaplaceSolver(GraphLaplStencil::MPVStencil * laplacianStencils, const lilimap & maplocals);
    virtual void updateLaplaceSolver(GraphLaplStencil::MPVStencil * laplacianStencils, const lilimap & maplocals);
    virtual void assignBCAndEvaluateRHS(std::size_t comp, bool unused, GraphLaplStencil::MPVStencil * borderLaplacianStencil,
                                        const lilimap & maplocals, dvector1D & rhs);
    virtual void solveLaplace(const dvector1D &rhs, dvector1D & result);

    // Reconstruct final result field
    virtual void reconstructResults(const dvector2D & results, const lilimap & mapglobals,  livector1D * markedcells = nullptr);

};

/*!
 * \class PropagateScalarField
 * \ingroup propagators
 * \brief Executable block that provides the computation of a scalar field
 * over a volume or surface mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given Dirichlet boundary conditions.
 *
 * Dirichlet boundary conditions are explicitly provided by the User,
   identifying boundaries through MimmoObject patches (all types allowed) and associating
   to their points the scalar fields as Dirichlet conditions on each patch.
 * A natural zero gradient like condition is automatically provided on unbounded borders.
 *
 * The block can perform multistep evaluation to relax field propagation
 *
 * Class/BaseManipulation Object specialization of class PropagateField
 * for the propagation in a volume/surface mesh of a scalar field.
 *
 * Ports available in PropagateScalarField Class :
 *
 *    =========================================================
 *
    |Port Input|||
    ||||
    | <B>PortType</B>| <B>variable/function</B>  |<B>DataType</B> |
    | M_GEOM         | setGeometry                           | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM2        | addDirichletBoundaryPatch             | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM3        | addDampingBoundarySurface             | (MC_SCALAR, MD_MIMMO_) |
    | M_GEOM7        | addNarrowBandBoundarySurface          | (MC_SCALAR, MD_MIMMO_) |
    | M_FILTER       | addDirichletConditions                | (MC_SCALAR, MD_MPVECFLOAT_)|

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_FILTER         | getPropagatedField                  | (MC_SCALAR, MD_MPVECFLOAT_) |

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
 * - <B>Damping</B>              : 1-true activate damping control, 0-false deactivate it.
 * - <B>DampingType</B>          : 0- distance control, 1-volume control.
 * - <B>DampingInnerDistance</B> : inner limit of damping function, if damping is active;
 * - <B>DampingOuterDistance</B> : outer limit of damping function, if damping is active;
 * - <B>DecayFactor</B>          : exponent to modulate damping function (as power of),
                                   if damping is active
   - <B>NarrowBand</B>           : 1-true activate narrow band control, 0-false deactivate it.
   - <B>NarrowBandWidth</B>      : size of narrow band (> 0.0), if NarrowBand is
                                   active.
   - <B>NarrowBandRelax</B>      : NBC relaxation parameter within [0,1], where
                                   1 means no relaxation, 0 full relaxation. Meaningful
                                   only if NarrowBand is active.
 * - <B>Tolerance</B>            : convergence tolerance for laplacian solver.
 * - <B>Print</B>                : print solver debug information, Active only
                                   if MIMMO_ENABLE_MPI is enabled in compilation.
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

    dmpvector1D* getPropagatedField();

    void    addDirichletConditions(dmpvector1D *bc);
    void    setSolverMultiStep(unsigned int sstep);

    //cleaners and setters
    virtual void clear();
    virtual void setDefaults();

    //execute
    void        execute();
    void swap(PropagateScalarField & x) noexcept;

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    int           m_nstep;               /**< multistep solver */

    virtual void buildPorts();
    virtual void plotOptionalResults();

private:
    dmpvector1D m_tempfield;             /**< temporary field storage for output purpose of getPropagatedField */
    std::vector<MimmoPiercedVector<std::array<double,1> > > m_tempDirichletBcs; /**< temporary input bc dirichlet storage for input managing purposes */

    //override interface methods, unmeaningful for the current class.
    void    setUpdateThreshold(double thres) override{
        BITPIT_UNUSED(thres);
    };
};


/*!
 * \class PropagateVectorField
 * \ingroup propagators
 * \brief Executable block that provides the computation of a 3D array field
 * over a volume/surface mesh. The field is calculated solving a Laplacian problem over
 * the mesh with given boundary conditions.
 *
 * Dirichlet boundary conditions are explicitly provided by the User,
   identifying boundaries through MimmoObject patches (all types allowed) and associating
   to their points the vector fields as Dirichlet conditions on each patch.

 * Optionally an impermeability-like/slip condition (as a zero vector field normal to
   boundary surface realized with a deformation reprojection onto a reference
   patch) can be imposed on chosen boundary patches. The nodes of the slip
   boundary patch are moved on a reference patch, that is completely independent
   from the slip boundary patches. If not provided by the User
   the reference surface is fixed as a copy of the reconstructed slip boundary patches list.
   <b>Note.</b> Currently, slip conditions are allowed only for bulk volume meshes.

 * The User can force the slip reference surface to be a plane
   by activating the related flag; in this case the mean plane defined over
   the slip boundary patches (the slip reference patch is useless) is used.
   <b>Note.</b> Currently, slip conditions are allowed only for bulk volume meshes.

 * Another option is to set periodic conditions on boundary patches. In this context,
   periodic means that the original boundary shape of the patch remains unaltered,
   but its nodes can move, constrained onto the patch itself. It is a "special" condition
   of slip patches, whose borders are fixed to zero Dirichlet conditions.
   <b>Note.</b> Currently, periodic conditions are allowed only for bulk volume meshes.

 * A natural zero gradient like condition is automatically provided on unbounded borders.
 *
 * The block can perform multistep evaluation to further relax the field propagation. In this case,
   on each step evaluation the vector field is applied on the bulk/boundaries to achieve a
   partial deformation of the mesh.
 *
 * Ports available in PropagateVectorField Class :
 *
 *    =========================================================
 *
    | Port Input|||
    ||||
    | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>             |
    | M_GEOM          | setGeometry                 | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM2         | addDirichletBoundaryPatch   |(MC_SCALAR, MD_MIMMO_)   |
    | M_GEOM3         | addDampingBoundarySurface   | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM4         | addSlipBoundarySurface      | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM5         | addPeriodicBoundarySurface  | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM6         | addSlipReferenceSurface     | (MC_SCALAR, MD_MIMMO_)  |
    | M_GEOM7         | addNarrowBandBoundarySurface| (MC_SCALAR, MD_MIMMO_)  |
    | M_GDISPLS       | addDirichletConditions      | (MC_SCALAR, MD_MPVECARR3FLOAT_)|

    |Port Output|||
    ||||
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS         | getPropagatedField   | (MC_SCALAR, MD_MPVECARR3FLOAT_) |

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
 * - <B>Damping</B>              : 1-true activate damping control, 0-false deactivate it.
 * - <B>DampingType</B>          : 0- distance control, 1-volume control.
 * - <B>DampingInnerDistance</B> : inner limit of damping function, if damping is active;
 * - <B>DampingOuterDistance</B> : outer limit of damping function, if damping is active;
 * - <B>DecayFactor</B>          : exponent to modulate damping function (as power of),
                                   if damping is active
   - <B>NarrowBand</B>           : 1-true activate narrow band control, 0-false deactivate it.
   - <B>NarrowBandWidth</B>      : size of narrow band (> 0.0), if NarrowBand is
                                   active.
   - <B>NarrowBandRelax</B>      : NBC relaxation parameter within [0,1], where
                                   1 means no relaxation, 0 full relaxation. Meaningful
                                   only if NarrowBand is active.
 * - <B>Tolerance</B>            : convergence tolerance for laplacian solver.
 * - <B>UpdateThres</B>          : lower threshold to internally mark cells whose
                                   field norm is above its value, for update purposes
 * - <B>Print</B>                : print solver debug information, Active only
                                   if MIMMO_ENABLE_MPI is enabled in compilation.
 *
 * Proper fo the class:
 * - <B>MultiStep</B> : got deformation in a finite number of substep of solution;
 * - <B>ForcePlanarSlip</B> : (for Quasi-Planar Slip Surface Only)  1- force the
                               class to treat slip surface as plane (without holes),
                               0-use slip reference surface as it is;
 *
 * Geometry, boundary surfaces, boundary condition values
 * for the target geometry have to be mandatorily passed through ports.
 *
 */
class PropagateVectorField: public mimmo::PropagateField<3> {

protected:

    int           m_nstep;  /**< multistep solver steps */

    bool m_forcePlanarSlip; /**< force slip surface to be treated as plane */
    std::unordered_set<MimmoSharedPointer<MimmoObject> > m_slipSurfaces;          /**< list of MimmoObject boundary patches where slip conditions are applied */
    std::unordered_set<MimmoSharedPointer<MimmoObject> > m_slipReferenceSurfaces; /**< list of MimmoObject boundary patches identifying slip reference surface on which the slip nodes of boundary patch are re-projected. */
    MimmoSharedPointer<MimmoObject> m_slipUniSurface;            /**< INTERNAL use. Final slip surface.*/
    MimmoPiercedVector<std::array<double, 3> > m_slip_bc_dir; /**< INTERNAL USE ONLY: Slip-type condition values on POINTS of the target volume mesh */

    std::unordered_set<MimmoSharedPointer<MimmoObject> > m_periodicSurfaces;   /**< MimmoObject boundary patch identifying periodic conditions */
    std::unordered_set<long> m_periodicBoundaryPoints;     /**< list of mesh nodes flagged as periodic */

private:
    std::array<double,3> m_AVGslipNormal;
    std::array<double,3> m_AVGslipCenter;

public:

    PropagateVectorField();
    PropagateVectorField(const bitpit::Config::Section & rootXML);
    virtual ~PropagateVectorField();
    PropagateVectorField(const PropagateVectorField & other);
    PropagateVectorField & operator=(PropagateVectorField other);
    void swap(PropagateVectorField & x) noexcept;

    void buildPorts();

    dmpvecarr3E * getPropagatedField();
    bool        isForcingPlanarSlip();

    void    addSlipBoundarySurface(MimmoSharedPointer<MimmoObject>);
    void    addSlipReferenceSurface(MimmoSharedPointer<MimmoObject>);
    void    addPeriodicBoundarySurface(MimmoSharedPointer<MimmoObject>);

    void    forcePlanarSlip(bool planar);
    void    addDirichletConditions(dmpvecarr3E * bc);

    void    setSolverMultiStep(unsigned int sstep);

    //cleaners and setters
    virtual void setDefaults();
    virtual void clear();

    //execute
    void        execute();

    //XML utilities from reading writing settings to file
    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    virtual bool checkBoundariesCoherence();

    virtual void subdivideBC();
    virtual void restoreBC();
    virtual void apply();
    virtual void restoreGeometry(bitpit::PiercedVector<bitpit::Vertex> & vertices);
    virtual void plotOptionalResults();

    virtual void propagateMaskMovingPoints(livector1D & vertexlist);

    virtual void assignBCAndEvaluateRHS(std::size_t comp, bool slipCorrect,
                                GraphLaplStencil::MPVStencil * borderLaplacianStencil,
                                const lilimap & maplocals,
                                dvector1D & rhs);

    virtual void computeSlipBCCorrector(const MimmoPiercedVector<std::array<double,3> > & guessSolutionOnPoint);

    void initializeSlipSurfaceAsPlane();
};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM2, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM3, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM4, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM5, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM6, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GEOM7, MC_SCALAR, MD_MIMMO_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_FILTER, MC_SCALAR, MD_MPVECFLOAT_,__PROPAGATEFIELD_HPP__)
REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_,__PROPAGATEFIELD_HPP__)

REGISTER(BaseManipulation, PropagateScalarField, "mimmo.PropagateScalarField")
REGISTER(BaseManipulation, PropagateVectorField, "mimmo.PropagateVectorField")

};

#include "PropagateField.tpp"

#endif /* __PROPAGATEFIELD_HPP__ */
