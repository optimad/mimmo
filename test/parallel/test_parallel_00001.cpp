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
 * \param[in] origin Origin point of the mesh.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(MimmoObject * mesh, std::array<double,3> origin = std::array<double,3>({0.,0.,0.})){

    if (mesh->getRank() == 0){
        double dx = 0.25, dy = 0.25;
        int nV, nC;
        //create vertexlist
        dvecarr3E vertex(35,{{0.,0.,0.}});
        livector2D conn(48, livector1D(3));

        for(int i=0; i<7; ++i){
            for(int j=0; j<5; j++){
                nV = 5*i + j;
                vertex[nV][0] = origin[0] + i*dx;
                vertex[nV][1] = origin[1] + j*dy;
                vertex[nV][2] = origin[2];
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

    mesh->updateAdjacencies();
    mesh->update();

    bool check = true;
    if (mesh->getRank() == 0){
        check = (mesh->getNCells() == 48) && (mesh->getNVertices() == 35);
    }

    return check;
}
// =================================================================================== //
/*
 * Test: testing project point method utility
 */
int test1() {

    // ====================================================================== //
    // Initialize the logger
    // ====================================================================== //
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    log::manager().initialize(log::COMBINED, true, nProcs, rank);
    log::cout().setVisibility(log::MASTER);

    /* Create target test mesh */
    MimmoSharedPointer<MimmoObject> target(new MimmoObject());
    if(!createMimmoMesh(target.get())) {
        return 1;
    }

    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partitionTarget = new mimmo::Partition();
    partitionTarget->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partitionTarget->setPlotInExecution(true);
    partitionTarget->setGeometry(target);
    partitionTarget->exec();
    target->getPatch()->write("target");

    /*Build the SkdTree of the geometry and retrieve it. */
    target->buildSkdTree();
    bitpit::PatchSkdTree* treeTarget = target->getSkdTree();


    /* Create selection mesh */
    MimmoSharedPointer<MimmoObject> selection(new MimmoObject());
    std::array<double,3> origin({0.6, 0.6, 0.01});
    if(!createMimmoMesh(selection.get(), origin)) {
        return 1;
    }

    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    mimmo::Partition* partitionSelection = new mimmo::Partition();
    partitionSelection->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partitionSelection->setPlotInExecution(true);
    partitionSelection->setGeometry(selection);
    partitionSelection->exec();
    selection->getPatch()->write("selection");

    /*Build the SkdTree of the geometry and retrieve it. */
    selection->buildSkdTree();
    bitpit::PatchSkdTree* treeSelection = selection->getSkdTree();

    /* Select by Patch. */
    std::vector<long> selectedCells;
    selectedCells = skdTreeUtils::selectByGlobalPatch(treeSelection, treeTarget, 2.0e-02);

    std::cout << " --- Selected cells --- " << std::endl;
    std::cout << "#" << rank << " selected cell size : " << selectedCells.size() << std::endl;
    for (long cell : selectedCells){
        std::cout << "#" << rank << " selected cell id : " << cell << std::endl;
    }
    std::cout << std::endl;

    bool check = true;
    std::set<long> result;
    result.insert(20);
    result.insert(21);
    result.insert(22);
    result.insert(23);
    result.insert(28);
    result.insert(29);
    result.insert(30);
    result.insert(31);
    result.insert(36);
    result.insert(37);
    result.insert(38);
    result.insert(39);
    result.insert(44);
    result.insert(45);
    result.insert(46);
    result.insert(47);
    for (long cell : selectedCells){
        check = check && (result.count(cell));
    }

    if (check)
        bitpit::log::cout() << " test passed " << std::endl;

    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    MPI_Init(&argc, &argv);

    /**<Calling mimmo Test routines*/
    int val = 1;
    try{
        val = test1() ;
    }
    catch(std::exception & e){
        std::cout<<"test_parallel_00001 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

    MPI_Finalize();

    return val;
}
