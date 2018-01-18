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

#ifndef __MIMMO_SYSTEM_HPP__
#define __MIMMO_SYSTEM_HPP__

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <bitpit_patchkernel.hpp>
#include <mimmoTypeDef.hpp>
#include <InOut.hpp>
#include <petscksp.h>

namespace mimmo{

typedef std::vector<std::vector<int>>  localivector2D;   /**< mimmo custom typedef*/
typedef std::vector<std::vector<double>>  localdvector2D;   /**< mimmo custom typedef*/
typedef std::vector<double>  localdvector1D;   /**< mimmo custom typedef*/

/*!
 * \struct KSPOptions
 * \ingroup system
 * \brief Info for PETSc linear system solution.
 */
struct KSPOptions {
    PetscInt restart;
    PetscInt levels;
    PetscInt overlap;
    PetscInt sublevels;
    PetscInt maxits;
    PetscScalar rtol;
    PetscScalar subrtol;
    bool nullspace;

    KSPOptions()
    : restart(50), levels(1), overlap(0), sublevels(4),
      maxits(10000), rtol(1.e-13), subrtol(1.e-13), nullspace(false)
    {
    }
};

/*!
 * \struct KSPStatus
 * \ingroup system
 * \brief Info for PETSc linear system execution output.
 */
struct KSPStatus {
    PetscErrorCode error;
    PetscInt its;
    KSPConvergedReason convergence;

    KSPStatus()
    : error(0), its(-1), convergence(KSP_DIVERGED_BREAKDOWN)
    {
    }
};

/*!
 * \class SystemSolver
 * \ingroup system
 * \brief Class handling solution of a linear system with PETSc
 * \TODO to be updated with new bitpit wrapper to Linear System Solver.
 */

class SystemSolver {

public:
    enum PivotType {
        PIVOT_NONE,  // Natural
        PIVOT_ND,    // Nested Dissection
        PIVOT_1WD,   // One-way Dissection
        PIVOT_RCM,   // Reverse Cuthill-McKee
        PIVOT_MD     // Quotient Minimum Degree
    };

    static void addInitOption(std::string options);
    static void addInitOptions(const std::vector<std::string> &options);

#if ENABLE_MPI==1
    SystemSolver(MPI_Comm comm, bool debug = false);
#else
    SystemSolver(bool debug = false);
#endif
    ~SystemSolver();

    void clear();
#if ENABLE_MPI==1
    void initialize(localivector2D &stencils, localdvector2D &weights,
            localdvector1D &rhs, std::unordered_set<long> ghosts, PivotType pivotType = PIVOT_NONE);
#else
    void initialize(localivector2D &stencils, localdvector2D &weights,
            localdvector1D &rhs, PivotType pivotType = PIVOT_NONE);
#endif
    void solve();
    void solve(std::vector<double> &solution, std::vector<double> &rhs);

    void dump(const std::string &directory, const std::string &prefix = "") const;

    PivotType getPivotType();

    KSPOptions & getKSPOptions();
    const KSPOptions & getKSPOptions() const;
    const KSPStatus & getKSPStatus() const;

    double * getRHSRawPtr();
    const double * getRHSRawPtr() const;
    const double * getRHSRawReadPtr() const;
    void restoreRHSRawPtr(double *raw_rhs);
    void restoreRHSRawReadPtr(const double *raw_rhs) const;

    double * getSolutionRawPtr();
    const double * getSolutionRawPtr() const;
    const double * getSolutionRawReadPtr() const;
    void restoreSolutionRawPtr(double *raw_solution);
    void restoreSolutionRawReadPtr(const double *raw_solution) const;

private:
    static int m_nInstances;
    static std::vector<std::string> m_options;

    bool m_initialized;
    PivotType m_pivotType;

    MPI_Comm m_communicator;

    long m_rowGlobalIdOffset;

    Mat m_A;
    std::unordered_map<long, double> m_A_rhs;
    Vec m_rhs;
    Vec m_solution;

    IS m_rpivot;
    IS m_cpivot;

    KSP m_KSP;
    KSPOptions m_KSPOptions;
    KSPStatus m_KSPStatus;

    void matrixInit(localivector2D &stencils);
    void matrixFill(localivector2D &stencils, localdvector2D &weights, localdvector1D &rhs);
    void matrixReorder();

#if ENABLE_MPI == 1
    void vectorsInit(std::unordered_set<long> ghosts);
#else
    void vectorsInit();
#endif
    void vectorsReorder(PetscBool inv);
    void vectorsFill(std::vector<double> &solution, std::vector<double> &rhs);
    void vectorsExport(std::vector<double> &solution);

    void pivotInit(PivotType pivotType);

    void KSPInit();

};

}
#endif
