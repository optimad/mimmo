
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

#include "StencilFunctions.hpp"
#include "bitpit_LA.hpp"
#include "bitpit_CG.hpp"

namespace mimmo {

namespace FVolStencil{

    /*!
     * Internal method function -> used by ComputeFVCellStencil.
     * Given a set of points P (already translated w.r.t stencil center)
     * and an initial guess of weights w associated to them
     * evaluate the gradient stencil weights with a Least Squared method.
     *
     * \param[in] P distance set
     * \param[in] w initial  guess weights
     * \param[out] weights evaluated gradient weights
     */
void computeWeightsWLS( const std::vector<std::vector<double>> &P,
                        const std::vector<double> &w,
                        std::vector<std::vector<double>> *gweights)
{
    // Evaluate the coefficient matrix
    int n = P.size() ;
    double epsilon = 1.e-12 ;

    std::array<bool,3> isSingular= {{true,true,true}};
    for(int i=0; i< n; ++i){
        for( int j=0; j<3; ++j){
            isSingular[j] &= std::abs(P[i][j]) < epsilon ;
        }
        if( !isSingular[0] && !isSingular[1] && !isSingular[2] ) break;
    }

    double A[3][3] = {{0.0, 0., 0.}, {0., 0.0, 0.}, {0., 0., 0.0}};

    for (int j = 0; j < 3; ++j) {
        for (int l = 0; l < 3; ++l) {
            for (int i = 0; i < n; ++i) {
                A[j][l] += w[i] * P[i][l] * P[i][j] ;
            }
        }
    }

    // regularize missing directions
    for(int i=0; i<3; ++i){
        A[i][i] += (double) isSingular[i];
    }

    // Invert the matrix
    double A_00_11 = A[0][0] * A[1][1];
    double A_00_12 = A[0][0] * A[1][2];
    double A_00_21 = A[0][0] * A[2][1];
    double A_00_22 = A[0][0] * A[2][2];
    double A_02_11 = A[0][2] * A[1][1];
    double A_02_20 = A[0][2] * A[2][0];
    double A_02_21 = A[0][2] * A[2][1];
    double A_01_12 = A[0][1] * A[1][2];
    double A_01_22 = A[0][1] * A[2][2];
    double A_10_01 = A[1][0] * A[0][1];
    double A_10_02 = A[1][0] * A[0][2];
    double A_10_21 = A[1][0] * A[2][1];
    double A_10_22 = A[1][0] * A[2][2];
    double A_11_20 = A[1][1] * A[2][0];
    double A_11_22 = A[1][1] * A[2][2];
    double A_12_20 = A[1][2] * A[2][0];
    double A_20_01 = A[2][0] * A[0][1];
    double A_20_11 = A[2][0] * A[1][1];
    double A_21_12 = A[2][1] * A[1][2];

    double det = A[0][0] * (A_11_22 - A_21_12) - A[0][1] * (A_10_22 - A_12_20) + A[0][2] * (A_10_21 - A_11_20);
    double inv_det = 1.0 / det;

    std::vector<std::vector<double>> inv_A(3, std::vector<double>(3, 0));
    inv_A[0][0] =   (A_11_22 - A_21_12) * inv_det;
    inv_A[0][1] = - (A_10_22 - A_12_20) * inv_det;
    inv_A[0][2] =   (A_10_21 - A_20_11) * inv_det;
    inv_A[1][0] = - (A_01_22 - A_02_21) * inv_det;
    inv_A[1][1] =   (A_00_22 - A_02_20) * inv_det;
    inv_A[1][2] = - (A_00_21 - A_20_01) * inv_det;
    inv_A[2][0] =   (A_01_12 - A_02_11) * inv_det;
    inv_A[2][1] = - (A_00_12 - A_10_02) * inv_det;
    inv_A[2][2] =   (A_00_11 - A_10_01) * inv_det;

    // Evaluate the weights
    gweights->clear() ;
    *gweights = bitpit::linearalgebra::matmulDiag(bitpit::linearalgebra::matmul(inv_A, bitpit::linearalgebra::transpose(P)), w) ;
}

/*!
 *  Evaluate the gradient Stencil on Center Cell for MimmoObject Volume mesh (type 2).
 *  Internally, it employs a Least Squared method to get the correct weights of the gradient.
 *   The method assumes target mesh has adjacencies already built: no check is performed in this sense.
 *
 *  \param[in] geo target mesh
 *  \param[in] updatelist (optional) if not null, perform computation only on the cell listed,
 *                                   if null compute on all cells.
 *  \return a CELL data pierced vector containing the gradient stencils
 */
MPVGradientUPtr   computeFVCellGradientStencil(MimmoObject & geo, const std::vector<long> * updatelist){

    MPVGradientUPtr result = std::unique_ptr<MPVGradient>(new MPVGradient(&geo, MPVLocation::CELL));
    if(geo.getType() != 2) return result; //this work only on Volume meshes

    livector1D list;
    if(!updatelist){
        list = geo.getCells().getIds();
    }else{
        list.insert(list.end(), updatelist->begin(), updatelist->end());
    }

    if(list.empty()) return result;

    result->reserve(list.size());

    auto patch = geo.getPatch();
    std::vector<double>                   ww;
    std::vector< std::vector<double> >    points;
    std::vector< std::vector<double> >    gweights;
    std::vector< long >                   neighs;
    std::array<double,3>                  dist;
    std::array<double,3>                  targetCentroid;
    long cellID;
    int ncellFaces;

    for(const long &cellID : list){

        targetCentroid = patch->evalCellCentroid(cellID);
        neighs.clear();

        //get all the face neighbours
        patch->findCellNeighs(cellID, i, &neighs); //--> getting only the real neighs. -1 neighs are ignored into the method.
        std::size_t ngsize = neighs.size();
        points.resize(ngsize, std::vector<double>(3));
        ww.resize(ngsize);

        std::vector< std::vector<double> >::iterator pointsItr = points.begin();
        std::vector< double >::iterator              wwItr = ww.begin();

        //get a loop on neighbours. We are assuming no internal boundaries are present in the mesh.
        for( long idN : neighs){
            // get the vector distance from neigh centroid up to current cell centroid.
            dist = patch->evalCellCentroid(idN) - targetCentroid;
            // calculate the initial weight
            *wwItr =  1.0/std::pow( norm2(dist), 1 ); //coefficient to be changed 1,3
            //copy the vector distance in points
            std::copy( dist.begin(), dist.end(), pointsItr->begin());
            //increment the iterators.
            ++wwItr;
            ++pointsItr;
        }

        //calculate the gradient stencil with LeastSquared method.
        computeWeightsWLS(points, ww, &gweights);

        //add an StencilVector Item to MimmoPiercedVector and push gweight into it.
        auto it = result->insert(cellID, bitpit::StencilVector());
        for(int i=0; i<ngsize; ++i){
            it->appendItem(neighs[i], {{gweights[0][i], gweights[1][i], gweights[2][i]}});
        }

        it->addComplementToZero(cellID); // adding the central node weight as minus sum of other weights.
        it->flatten();

    }
    return result;
}

/*!
 *  Evaluate the gradient Stencil on ALL Face Cells for MimmoObject Volume mesh (type 2).
 *  The method assumes target mesh has interfaces already built: no check is performed in this sense.
 *
 * The method first calculate the gradient stencil on center cells with computeFVCellGradientStencil,
 * if not already provided by the User.
 * Then it use these info to distribute them on interface. Face non orthogonality corrections are directly
 * performed on the final stencil (amend and recalculate the gradient normal part to interface).
 *
 * BEWARE: Boundary Interfaces Stencils reports the a Stencil with a Homogeneous Neumann condition on it.
 *         You need to sum the right correction to recover the proper condition.
 *
 * \param[in] geo target mesh
 * \param[in] cellGradientStencil (optional) gradient stencil calculated on ALL cell centers + ghosts.
 * \return a INTERFACE data pierced vector containing the gradient stencil
 */
MPVGradientUPtr   computeFVFaceGradientStencil(MimmoObject & geo, MPVGradient * cellGradientStencil ){

    MPVGradientUPtr calculateCellGradStencil;

    MPVGradient * cgs = nullptr;
    if(!cellGradientStencil){
        calculateCellGradStencil = computeFVCellGradientStencil(geo);
        cgs = calculateCellGradStencil.get();
    }else if(cellGradientStencil->getGeometry() != &geo || cellGradientStencil->getDataLocation() != MPVLocation::CELL) {//check if mpv structure is messed
        calculateCellGradStencil = computeFVCellGradientStencil(geo);
        cgs = calculateCellGradStencil.get();
    }else{
        cgs = cellGradientStencil;
    }
    return updateFVFaceGradientStencil(geo, *cgs);
}

/*!
 *  Re-evaluate the gradient Stencil on Face Cells for a selected list of MimmoObject Volume mesh (type 2) cells.
 *  The method assumes target mesh has interfaces already built: no check is performed in this sense.
 *
 * The method first re-evaluate the center cell gradient stencil on a selected list of cells provided by the User
 * Then it use these info to re-evaluate stencils on all interfaces related to the cells of the list.
 * Face non orthogonality corrections are directly performed on the final stencil (amend and recalculate
 * the gradient normal part to interface).
 *
 * BEWARE: During this update, Boundary Interfaces Stencils (if any) reports the a Stencil with a Homogeneous Neumann condition on it.
 *         You need to sum the right correction to recover the proper condition.
 *
 * \param[in] geo target mesh
 * \param[in] list of cell IDs you need to update.
 * \return a INTERFACE data pierced vector containing gradient stencils at all interfaces relative to the selected cells.
 */
MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, const std::vector<long> & list ){

    MPVGradientUPtr calculateCellGradStencil = computeFVCellGradientStencil(geo, &list);
    return updateFVFaceGradientStencil(geo, *(calculateCellGradStencil.get()));
}

/*!
 *  Re-evaluate the gradient Stencil on Face Cells for a selected list of MimmoObject Volume mesh (type 2) cells.
 *  The method assumes target mesh has interfaces already built: no check is performed in this sense.
 *
 * The method use the already evaluated center cell gradient stencil on a selected list of cells provided by the User
 * Then it re-evaluates stencils on all interfaces related to the cells of the list.
 * Face non orthogonality corrections are directly performed on the final stencil (amend and recalculate
 * the gradient normal part to interface).
 *
 * BEWARE: During this update, Boundary Interfaces Stencils (if any) reports the a Stencil with a Homogeneous Neumann condition on it.
 *         You need to sum the right correction to recover the proper condition.
 *
 * \param[in] geo target mesh
 * \param[in] list of Center Cell Gradient stencils involved into update.
 * \return a INTERFACE data pierced vector containing gradient stencils at all interfaces relative to the selected cells.
 */
MPVGradientUPtr   updateFVFaceGradientStencil(MimmoObject & geo, MPVGradient & cellGradientStencil ){

    MPVGradientUPtr result = std::unique_ptr<MPVGradient>(new MPVGradient(&geo, MPVLocation::INTERFACE));

    if(geo.getType() != 2) return result; //this work only on Volume meshes

    MPVGradient * cgs = &cellGradientStencil;
    //look if you have an empty cc gradient pool
    if(cgs->isEmpty()){
        return result;
    }

    result->reserve(geo.getPatch()->getInterfaceCount());

    //now i have the gradient stencil on cell center.
    bitpit::PiercedVector<bitpit::Interface> & interfaces = geo.getInterfaces();
    bitpit::PiercedVector<bitpit::Cell> & cells = geo.getCells();

    long ownerID, neighID;
    std::array<double,3> ownerCentroid, neighCentroid, interfaceCentroid, interfaceNormal;
    std::array<double,3> ownerProjection, neighProjection;

    bitpit::StencilScalar ownerProjStencil, neighProjStencil;
    bitpit::StencilVector ownerStencil, neighStencil, avgStencil;
    double minDistance;

    livector1D targetlist = geo.getInterfaceFromCellList(cgs->getIds(), false);

    for(long & interfid : targetlist){
        const bitpit::Interface & interface = interfaces.at(interfid);
        ownerID = interface.getOwner();
        neighID = interface.getNeigh();
        interfaceNormal = geo.evalInterfaceNormal(interfid);

        if(neighID < 0){ //interface is a border
            if(!cells.at(ownerID).isInterior()){
                continue; //skip if the owner is a ghost
            }
            // insert the Neumann Homogeneous face gradient @ the boundary.
            result->insert(interfid,
                           computeBorderFaceGradient(interfaceNormal, cgs->at(ownerID)) );
            // you need later to correct this stencil to provide the opportune bc.
        }else{
            if(!cells.at(ownerID).isInterior() && !cells.at(neighID).isInterior()){
                continue; //skip, both owner and neigh are ghosts.
            }

            //clear the stencil but not release their memory reservation. This should be a little more
            // performant, but i cannot say right now. VERIFY-->
            avgStencil.clear(false);
            ownerStencil.clear(false);
            neighStencil.clear(false);
            ownerProjStencil.clear(false);
            neighProjStencil.clear(false);

            //evaluate the gradient stencil on interface, accounting for non-orthogonality correction.
            //get interface centroid
            interfaceCentroid = geo.evalInterfaceCentroid(interfid);

            // get owner and neighbor centroids
            ownerCentroid = geo.getPatch()->evalCellCentroid(ownerID);
            neighCentroid = geo.getPatch()->evalCellCentroid(neighID);

            ownerStencil = cgs->at(ownerID);
            neighStencil = cgs->at(neighID);

            // find projection of owner and neigh cell centroid over the line defined by interface centroid and normal.
            bitpit::CGElem::distancePointLine(ownerCentroid, interfaceCentroid, interfaceNormal, ownerProjection);
            bitpit::CGElem::distancePointLine(neighCentroid, interfaceCentroid, interfaceNormal, neighProjection);

            //evaluate the minimum distances of owner/neighbor from interface centroid
            minDistance = std::min( norm2(ownerProjection - interfaceCentroid),
                                    norm2(neighProjection - interfaceCentroid) );

            // re-evaluate the new positions of owner/neighbor projections
            ownerProjection = interfaceCentroid - minDistance *interfaceNormal;
            neighProjection = interfaceCentroid + minDistance *interfaceNormal;

            // evaluate the scalar stencil of owner and neighbor projection.
            // The meaning of such stencil is the following: this is a stencil
            // that get an approximation of the field A into the new position
            // xxxProjection starting from the Grad(A) defined in the position
            //xxxCentroid.
            ownerProjStencil.appendItem(ownerID, 1.0);
            ownerProjStencil += dotProduct(ownerStencil, ownerProjection-ownerCentroid);
            ownerProjStencil.flatten();

            neighProjStencil.appendItem(neighID, 1.0);
            neighProjStencil += dotProduct(neighStencil, neighProjection-neighCentroid);
            neighProjStencil.flatten();

            //calculate stencil @ interface as the plain average of owner and neigh ccell
            // gradient stencil
            avgStencil = 0.5*(ownerStencil + neighStencil);
            avgStencil.flatten();
            // push it in MPV.
            auto it = result->insert(interfid, avgStencil);

            // amend the normal part, aka the actual stencil projected
            // in the normal direction at interface
            *it -= dotProduct(avgStencil, interfaceNormal) * interfaceNormal;

            // add again the normal part recalculated properly with proj Stencil.
            *it += (neighProjStencil -ownerProjStencil) *
                  (interfaceNormal/(2.0*minDistance)) ;
            //that's it.
        }
    } //end on interface loop

    return result;
}

/*!
 * Evaluate the Gradient Stencil at a generic Boundary Interface using a homogeneous Neumann condition.
 * The method corrects a CenterCell Gradient Stencil referred to the Owner cell of
 * the interface, erasing the cell gradient normal part to the interface and
 * substituting it with a zero onstant part
 *
 * Please note, you can impose your right boundary condition later adding an opportune
 * correction (a real prediction of the normal part). see correctionNeumannBCFaceGradient
 * or correctionDirichletBCFaceGradient
 *
 * \param[in] interfaceNormal normal defined on the BoundaryInterface
 * \param[in] CCellOwnerStencil cell gradient stencil referred to the owner cell
 * \return the homogeneous neumann border stencil.
 */

bitpit::StencilVector computeBorderFaceGradient(const std::array<double,3> & interfaceNormal,
                                                const bitpit::StencilVector & CCellOwnerStencil)
{
    bitpit::StencilVector correction(CCellOwnerStencil);
    //epurate the normal part form the cell gradient
    correction -= dotProduct(CCellOwnerStencil, interfaceNormal) * interfaceNormal;

    return correction;
}

/*!
 * Compute correction to border face gradient to recover a Generic Neumann condition @ the face.
 * Border face Generic Neumann gradient stencil can be decomposed as a sum two parts:
 *
 * - The Homogeneous Neumann part, i.e. the Gradient @ Owner cell centroid
 *   without its contribute along normal interface
 * - The correction part, i.e. the value of the directional derivative along the
 *   interface normal.
 *
 * Please note, this correction is "seen" as a constant part by the stencil.
 *
 * \param[in] neuval value of the Neumann Condition
 * \param[in] interfaceNormal normal defined on the BoundaryInterface
 * \return the neumann correction.
 */

bitpit::StencilVector correctionNeumannBCFaceGradient(const double & neuval,
                                                   const std::array<double,3> & interfaceNormal)
{
    bitpit::StencilVector correction;
    correction.sumConstant(neuval*interfaceNormal);
    return correction;
}
/*!
 * Compute correction to border face gradient to recover a Dirichlet condition @ the face.
 * Border face Dirichlet gradient stencil can be decomposed as a sum two parts:
 *
 * - The Homogeneous Neumann part, i.e. the Gradient @ Owner cell centroid
 *   without its contribute along normal interface
 * - The correction part, i.e. a finite difference estimation of the gradient stencil at
 *  the border interface, where the solution is imposed.
 *
 * The finite-difference correction is [(U_H - U_D) / distD] * interfaceNormal.
 * U_D is the value of the dirichlet condition, H is a point internal at Owner Cell
 * at distance distD form the interface centroid and on the interface normal direction,
 * distD is a charateristic length of the Owner cell (typically the radius of the element
 * in-sphere, incircle), U_H an estimation of field value @ H, evaluated
 * starting from Cell Gradient stencil @ Owner.
 *
 * Please note, the part (U_D/distD)*interfaceNormal of the correction is seen as constant
 * by the stencil.
 *
 * \param[in] dirval value of the Dirichlet Condition U_D
 * \param[in] ownerID id the the ownerCell
 * \param[in] ownerCentroid centroid of the owner of the BoundaryInterface
 * \param[in] interfaceCentroid centroid of the BoundaryInterface
 * \param[in] interfaceNormal normal defined on the BoundaryInterface
 * \param[in] distD charateristic length to pose the internal point H from interface Centroid.
 * \param[in] CCellOwnerStencil cell gradient stencil referred to the owner cell
 * \return the new corrected stencil.
 */

bitpit::StencilVector correctionDirichletBCFaceGradient(const double &dirval,
                                                     const long & ownerID,
                                                     const std::array<double,3> & ownerCentroid,
                                                     const std::array<double,3> & interfaceCentroid,
                                                     const std::array<double,3> & interfaceNormal,
                                                     const double &distD,
                                                     const bitpit::StencilVector & CCellOwnerStencil)
{

    bitpit::StencilVector correction;
    double distance(distD);
    // evaluate the estimation U_H field value obtained used the ccell gradient stencil dot the
    // distance vector from owner Centroid to point H.
    // the point H internal to the Owner is obtained as interfaceCentroid - distD*interfaceNormal
    //(interface normal goes always from Owner to Neighbor).
    bitpit::StencilScalar projpart;
    projpart.appendItem(ownerID, 1.0);
    projpart += dotProduct(CCellOwnerStencil, interfaceCentroid - distance*interfaceNormal - ownerCentroid);

    // build the new normal gradient part.
    // The U_H estimation enter into the stencil
    correction = -1.0 * (projpart/distance) * interfaceNormal;
    // The Dirichlet part enter as constant into the stencil
    correction.sumConstant((dirval/distance) * interfaceNormal);

    //go in peace.
    return correction;
}


/*!
 * The method computes the Laplacian stencils on cells, once the gradient stencils on
 * interfaces and specialized the boundary conditions on boundary interfaces are available.
 * The resulting laplacian stencils will be available on mesh interior cells of the mesh whose
 * interfaces are all reported in the faceGradientStencil list.
 * Providing the right set of interfaces you can use this method both to compute laplacian stencils on
 * the whole mesh or as updater of its subportions.
 * Mesh Interfaces are supposed built,  no control is done in such sense.
 * Only for MimmoObject Volume mesh (type=2).
 *
 * Once computed, the stencils are then divided by their ref cell volume.
 *
 * Both faceGradientStencil and diffusivity(if any) need to be referred to the same mesh.
 * In case of subportions updater usage, be sure diffusivity includes info on all the cells involved.

 * \param[in] faceGradientStencil gradient stencil defined on effective INTERFACES (bc included, no ghost-ghost).
 * \param[in] diffusivity (optional) impose a diffusivity field on center CELLS (provided also on ghosts)
 */

MPVDivergenceUPtr computeFVLaplacianStencil (MPVGradient & faceGradientStencil,
                                             MimmoPiercedVector<double> * diffusivity)
{

    MimmoObject * geo = faceGradientStencil.getGeometry();
    // prepare and allocate the result list for laplacian stencils.
    MPVDivergenceUPtr result = MPVDivergenceUPtr(new MPVDivergence(geo, MPVLocation::CELL));
    if(geo->getType() != 2) return result; //only for Volume Meshes

    result->reserve(geo->getPatch()->getInternalCount());
    for(auto it = geo->getPatch()->internalBegin(); it != geo->getPatch()->internalEnd(); ++it){
        result->insert(it->getId(), bitpit::StencilScalar());
    }

    //loop on interfaces
    double locdiff;
    std::array<double,3> interfaceNormal;
    long ownerID, neighID, interfaceID;
    bitpit::StencilScalar faceGradientNormal;
    bitpit::PiercedVector<bitpit::Cell> & cells = geo->getCells();
    std::unordered_set<long> cellInvolved;

    for(auto it = faceGradientStencil.begin(); it != faceGradientStencil.end(); ++it){

        interfaceID = it.getId();
        const bitpit::Interface & interface = geo->getInterfaces().at(interfaceID);

        interfaceNormal = geo->evalInterfaceNormal(interfaceID);
        faceGradientNormal = dotProduct(*it, interfaceNormal);

        ownerID = interface.getOwner();
        neighID = interface.getNeigh();

        locdiff = geo->evalInterfaceArea(interfaceID);

        if(neighID < 0){
            //borderInterface - owner is an Interior cell for sure
            // so check if you have non unitary diffusivity
            if(diffusivity) locdiff *= diffusivity->at(ownerID);
            // sum up the flux to the owner divergence
            result->at(ownerID) += locdiff * faceGradientNormal;
            cellInvolved.insert(ownerID);
        }else{
            // its an internal interface. Owner ot Neighbor can be ghost.
            // Do not account the flux on the ghost cells.

            // so check if you have non unitary diffusivity
            if(diffusivity){
                locdiff *= 0.5*(diffusivity->at(ownerID) + diffusivity->at(neighID));
            }

            // control and sum up flux to Owner divergence
            if(cells.at(ownerID).isInterior()){
                result->at(ownerID) += locdiff * faceGradientNormal;
                cellInvolved.insert(ownerID);
            }

            // control and substract up flux to Neighbor divergence
            if(cells.at(neighID).isInterior()){
                result->at(neighID) -= locdiff * faceGradientNormal;
                cellInvolved.insert(neighID);
            }
        }
    }

    //clean result and fit its dimension to the only cells involved
    result->squeezeOutExcept(std::vector<long>(cellInvolved.begin(), cellInvolved.end()), true);

    //divide stencils by their cell volume and optimize
    for(MPVDivergence::iterator it = result->begin(); it !=result->end(); ++it){
        *it /= geo->evalCellVolume(it.getId());
        //it->optimize();
    }

    return result;
}


} //end of stencil function namespace

} //end of mimmo namespace
