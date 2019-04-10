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

# include "bitpit_patchkernel.hpp"
# include "bitpit_surfunstructured.hpp"
# include "surface_skd_tree.hpp"
# include "volume_skd_tree.hpp"

namespace mimmo{


/*!
 * \brief Utilities employing SkdTree.
 * \ingroup core
 */
namespace skdTreeUtils{

    double distance(std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, long &id, double &r);
    double signedDistance(std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, long &id, std::array<double,3> &n, double &r);
    std::vector<long> selectByPatch(bitpit::PatchSkdTree *selection, bitpit::PatchSkdTree *target, double tol = 1.0e-04);
    void extractTarget(bitpit::PatchSkdTree *target, const std::vector<const bitpit::SkdNode*> & leafSelection, std::vector<long> &extracted, double tol);
    std::array<double,3> projectPoint(std::array<double,3> *P_, bitpit::PatchSkdTree *bvtree_, double r_ = 1.0e+18);
    long locatePointOnPatch(const std::array<double, 3> &point, bitpit::PatchSkdTree &tree);
    long closestCellToPoint(const std::array<double, 3> &point, bitpit::PatchSkdTree &tree);
}; //end namespace skdTreeUtils

} //end namespace mimmo

#endif
