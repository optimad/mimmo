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
# ifndef __SKDTREEUTILS_HPP__
# define __SKDTREEUTILS_HPP__

# include <patch_skd_tree.hpp>
# include <surfunstructured.hpp>

namespace mimmo{


/*!
 * \brief Utilities employing SkdTree.
 * \ingroup core
 */
namespace skdTreeUtils{

    double distance(const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long &id, double r);
    double signedDistance(const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long &id, std::array<double,3> &normal, double r);
    void distance(std::size_t nP, const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long *id, double *distances, double r);
    void distance(std::size_t nP, const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long *id, double *distances, double *r);
    void signedDistance(std::size_t nP, const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long *ids, double *distances, std::array<double,3> *normals, double r);
    void signedDistance(std::size_t nP, const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, long *ids, double *distances, std::array<double,3> *normals, double *r);
    std::vector<long> selectByPatch(bitpit::PatchSkdTree *selection, bitpit::PatchSkdTree *target, double tol = 1.0e-04);
    void extractTarget(bitpit::PatchSkdTree *target, const std::vector<const bitpit::SkdNode*> & leafSelection, std::vector<long> &extracted, double tol);
    std::array<double,3> projectPoint(const std::array<double,3> *point, const bitpit::PatchSkdTree *tree, double r = std::numeric_limits<double>::max());
    void projectPoint(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, std::array<double,3> *projected_points, long *ids, double r = std::numeric_limits<double>::max());
    void projectPoint(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, std::array<double,3> *projected_points, long *ids, double* r);
    long locatePointOnPatch(const std::array<double, 3> &point, const bitpit::PatchSkdTree *tree);

    std::array<double, 3> computePseudoNormal(const std::array<double,3> &point, const bitpit::SurfUnstructured *surface_mesh, long id);
    bool checkPointBelongsToCell(const std::array<double, 3> &point, const bitpit::SurfUnstructured *surface_mesh, long id);

#if MIMMO_ENABLE_MPI
    void globalDistance(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, double *distances, double r, bool shared = false);
    void globalDistance(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, double *distances, double* r, bool shared = false);
    void signedGlobalDistance(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, std::array<double,3> *normals, double *distances, double r, bool shared = false);
    void signedGlobalDistance(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, std::array<double,3> *normals, double *distances, double* r, bool shared = false);
    void projectPointGlobal(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, std::array<double,3> *projected_points, long *ids, int *ranks, double r = std::numeric_limits<double>::max(), bool shared = false);
    void projectPointGlobal(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, std::array<double,3> *projected_points, long *ids, int *ranks, double* r, bool shared = false);
    void locatePointOnGlobalPatch(std::size_t nP, const std::array<double,3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, bool shared = false);
    std::vector<long> selectByGlobalPatch(bitpit::PatchSkdTree *selection, bitpit::PatchSkdTree *target, double tol = 1.0e-04);
    void extractTarget(bitpit::PatchSkdTree *target, const std::vector<bitpit::SkdBox> &leafSelectionBoxes, std::vector<long> &extracted, double tol);
    void findSharedPointClosestGlobalCell(std::size_t nPoints, const std::array<double, 3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, double *distances, double r = std::numeric_limits<double>::max());
    void findSharedPointClosestGlobalCell(std::size_t nPoints, const std::array<double, 3> *points, const bitpit::PatchSkdTree *tree, long *ids, int *ranks, double *distances, double *r);
#endif

}; //end namespace skdTreeUtils

} //end namespace mimmo

#endif
