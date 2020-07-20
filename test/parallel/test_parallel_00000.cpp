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

#include "mimmo_core.hpp"
#include "SkdTreeUtils.hpp"
#include "Partition.hpp"
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;


/*
 * Creating surface triangular mesh and return it in a MimmoObject.
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(MimmoObject * mesh){

    if (mesh->getRank() == 0){
        double dx = 0.25, dy = 0.25;
        int nV, nC;
        //create vertexlist
        dvecarr3E vertex(35,{{0.0,0.0,0.0}});
        livector2D conn(48, livector1D(3));

        for(int i=0; i<7; ++i){
            for(int j=0; j<5; j++){
                nV = 5*i + j;
                vertex[nV][0] = i*dx;
                vertex[nV][1] = j*dy;
            }
        }

        for(int j=0; j<4; ++j){
            for(int i=0; i<3; ++i){
                nC = 8*i + 2*j;

                conn[nC][0] = 5*i + j;
                conn[nC][1] = 5*(i+1) + j;
                conn[nC][2] = 5*i + j+1;

                conn[nC+1][0] = 5*(i+1) + j;
                conn[nC+1][1] = 5*(i+1) + j+1;
                conn[nC+1][2] = 5*i + j+1;
            }
        }

        for(int j=0; j<4; ++j){
            for(int i=3; i<6; ++i){
                nC = 8*i + 2*j;

                conn[nC][0] = 5*i + j;
                conn[nC][1] = 5*(i+1) + j;
                conn[nC][2] = 5*(i+1) + j+1;

                conn[nC+1][0] = 5*i + j;
                conn[nC+1][1] = 5*(i+1) + j+1;
                conn[nC+1][2] = 5*i + j+1;
            }
        }

        mesh->getVertices().reserve(35);
        mesh->getCells().reserve(48);

        //fill the mimmoObject;
        long cV=0;
        for(auto & val: vertex){
            mesh->addVertex(val, cV);
            cV++;
        }

        long cC=0;
        bitpit::ElementType eltype = bitpit::ElementType::TRIANGLE;
        for(auto & val: conn){
            mesh->addConnectedCell(val, eltype, cC, mesh->getRank());
            cC++;
        }
    }

    mesh->buildAdjacencies();

    bool check = true;
    if (mesh->getRank() == 0){
        check = (mesh->getNCells() == 48) && (mesh->getNVertices() == 35);
    }

    mesh->getPatch()->write("support.0");

    return check;
}
// =================================================================================== //
/*
 * Test: testing project point method utility
 */
int test0() {

    // ====================================================================== //
    // Initialize the logger
    // ====================================================================== //
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::MASTER);

    /* Create test mesh */
    MimmoSharedPointer<MimmoObject> mimmo0(new MimmoObject());
    if(!createMimmoMesh(mimmo0.get())) {
        return 1;
    }

    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPlotInExecution(true);
    partition->setGeometry(mimmo0);
    partition->exec();
    mimmo0->getPatch()->write("support.1");

    /*Build the SkdTree of the geometry and retrieve it. */
    mimmo0->buildSkdTree();
    bitpit::PatchSkdTree* tree = mimmo0->getSkdTree();

    /* Initialize points. */
    bitpit::log::cout() << " --- POINTS --- " << std::endl;
    std::size_t np = 5;
    std::vector<std::array<double,3>> points(np);
    for (std::size_t i = 0; i < np; i++){
        std::array<double,3> & point = points[i];
        double x = 1.5/double(np)*double(i);
        double y = 1./double(np)*double(i);
        double z = 1.;
        point = std::array<double,3>({x, y, z});
        bitpit::log::cout()<< point <<std::endl;
    }
    bitpit::log::cout() << std::endl;

    /* Project points on surface. */
    bitpit::log::cout() << " --- PROJECTED POINTS --- " << std::endl;
    std::vector<std::array<double,3>> array_projected_points(np);
    std::vector<long> projectCellIds(np, bitpit::Cell::NULL_ID);
    std::vector<int> ranks(np, -1);
    double radius = std::numeric_limits<double>::max();
    skdTreeUtils::projectPointGlobal(np, points.data(), tree, array_projected_points.data(), projectCellIds.data(), ranks.data(), radius, true);

    for (auto & point : array_projected_points)
        bitpit::log::cout()<< point <<std::endl;
    bitpit::log::cout() << std::endl;


    /* Locate points on surface patch. */
    bitpit::log::cout() << " --- LOCATE POINTS --- " << std::endl;
    std::vector<long> locateCellIds(np, bitpit::Cell::NULL_ID);
    skdTreeUtils::locatePointOnGlobalPatch(np, array_projected_points.data(), tree, locateCellIds.data(), ranks.data(), true);

    for (std::size_t i = 0; i < np; i++){
        bitpit::log::cout()<< " projected cell id : " << projectCellIds[i] << "  locate cell id : " << locateCellIds[i] <<std::endl;
    }

    bool check = true;
    for (std::size_t i = 0; i < np; i++){
        check = check && (std::abs(locateCellIds[i] - projectCellIds[i]) == 0);
    }
    bitpit::log::cout() << std::endl;


    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    MPI_Init(&argc, &argv);

    /**<Calling mimmo Test routines*/
    int val = 1;
    try{
        val = test0() ;
    }
    catch(std::exception & e){
        std::cout<<"test_parallel_00000 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

    MPI_Finalize();

    return val;
}
