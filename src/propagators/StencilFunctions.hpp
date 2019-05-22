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
#ifndef __STENCILFUNCTIONS_HPP__
#define __STENCILFUNCTIONS_HPP__

#include "bitpit_discretization.hpp"
#include "MimmoObject.hpp"
#include "MimmoPiercedVector.hpp"


namespace mimmo{

/*!
 * \brief Collection of methods to precompute stencils of math operators (Laplacian, Gradients, etc...) on mimmo meshes
 * for finite volumes methods
 * \ingroup propagators
 */
namespace FVolStencil{
    typedef MimmoPiercedVector<bitpit::StencilVector> MPVGradient; /**< typedef to save stencil of a Grad stuff */
    typedef MimmoPiercedVector<bitpit::StencilScalar> MPVDivergence; /**< typedef to save stencil of Div stuff */
    typedef std::unique_ptr< MPVGradient > MPVGradientUPtr;  /**< typedef to instantiate unique ptr of Grad stencil */
    typedef std::unique_ptr< MPVDivergence > MPVDivergenceUPtr;/**< typedef to instantiate unique ptr of Div stencil */

    void computeWeightsWLS( const std::vector<std::vector<double>> &P,
                            const std::vector<double> &w,
                            std::vector<std::vector<double>> *weights);

    MPVGradientUPtr   computeFVCellGradientStencil(MimmoObject & geo, const std::vector<long> * updatelist = nullptr);

    MPVGradientUPtr   computeFVFaceGradientStencil(MimmoObject & geo, MPVGradient * cellGradientStencil = nullptr);
    MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, MPVGradient & cellGradientStencil);
    MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, const std::vector<long> & list);

    bitpit::StencilVector computeBorderFaceGradient(const std::array<double,3> & interfaceNormal,
                                                    const bitpit::StencilVector & CCellOwnerStencil);

    bitpit::StencilVector correctionNeumannBCFaceGradient(const double & neuval,
                                                          const std::array<double,3> & interfaceNormal);

    bitpit::StencilVector correctionDirichletBCFaceGradient(const double &dirval,
                                                            const long & ownerID,
                                                            const std::array<double,3> & ownerCentroid,
                                                            const std::array<double,3> & interfaceCentroid,
                                                            const std::array<double,3> & interfaceNormal,
                                                            const double &distD,
                                                            const bitpit::StencilVector & CCellOwnerStencil);

    MPVDivergenceUPtr computeFVLaplacianStencil (MPVGradient & faceGradientStencil, double tolerance = 1.0e-12,
                                                 MimmoPiercedVector<double> * diffusivity = nullptr);


};//end namespace FVolStencil


/*!
 * \brief Collection of methods to precompute stencils of math operators (Laplacian, Gradients, etc...) on mimmo meshes
 * for Graph Laplace deiscretization
 * \ingroup propagators
 */
namespace GraphLaplStencil{
    typedef MimmoPiercedVector<bitpit::StencilScalar> MPVStencil; /**< typedef to save stencil of a Grad stuff */
    typedef std::unique_ptr<MPVStencil> MPVStencilUPtr;  /**< typedef to instantiate unique ptr of Grad stencil */

//    void computeWeightsWLS( const std::vector<std::vector<double>> &P,
//                            const std::vector<double> &w,
//                            std::vector<std::vector<double>> *weights);

//    MPVGradientUPtr   computeStencils(MimmoObject & geo, const std::vector<long> * updatelist = nullptr);

//    MPVGradientUPtr   computeFVFaceGradientStencil(MimmoObject & geo, MPVGradient * cellGradientStencil = nullptr);
//    MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, MPVGradient & cellGradientStencil);
//    MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, const std::vector<long> & list);

//    bitpit::StencilVector computeBorderFaceGradient(const std::array<double,3> & interfaceNormal,
//                                                    const bitpit::StencilVector & CCellOwnerStencil);

//    bitpit::StencilVector correctionNeumannBCFaceGradient(const double & neuval,
//                                                          const std::array<double,3> & interfaceNormal);

//    bitpit::StencilVector correctionDirichletBCFaceGradient(const double &dirval,
//                                                            const long & ownerID,
//                                                            const std::array<double,3> & ownerCentroid,
//                                                            const std::array<double,3> & interfaceCentroid,
//                                                            const std::array<double,3> & interfaceNormal,
//                                                            const double &distD,
//                                                            const bitpit::StencilVector & CCellOwnerStencil);

    MPVStencilUPtr computeLaplacianStencils(MimmoObject & geo, double tolerance = 1.0e-12,
                                                 MimmoPiercedVector<double> * diffusivity = nullptr);

    MPVStencilUPtr computeLaplacianStencils(MimmoObject & geo, std::vector<long>* nodesList, double tolerance,
                                                 MimmoPiercedVector<double> * diffusivity);

};//end namespace stencilFunction

}; //end namespace mimmo

#endif /* __STENCILFUNCTIONS_HPP__ */
