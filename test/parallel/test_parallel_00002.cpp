/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#include "Partition.hpp"
#include <random>
#include <exception>

/*
 * Creating a multi-rank Point Cloud MimmoObject in unit cube.
   Mesh Data will be kept on master rank 0, leaving the other ranks empty.
   \param[in] np number of desired points
   \param[in] log reference to logger
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoPCMesh(int np, bitpit::Logger & log, mimmo::MimmoObject * mesh){

    if (mesh->getRank() == 0){

        dvecarr3E points(np, {{0.0,0.0,0.0}});
        std::random_device rd;
        unsigned int seed = rd(); //fix it up to any unsigned int for reproducibility issue.
        std::mt19937_64 rgen(seed);
        std::uniform_real_distribution<double> distr(0.0, 1.0);

        for(std::array<double,3> & p : points){
            p[0] = distr(rgen);
            p[1] = distr(rgen);
            p[2] = distr(rgen);
        }


        mesh->getVertices().reserve(np);
        mesh->getCells().reserve(np);

        //fill the mimmoObject;
        long cV(0), cC(0);
        bitpit::ElementType eltype = bitpit::ElementType::VERTEX;
        for(std::array<double,3> & val: points){
            mesh->addVertex(val, cV);
            mesh->addConnectedCell(std::vector<long>(1,cV), eltype, cC, mesh->getRank());
            cV++;
            cC++;
        }
    }

    mesh->updateAdjacencies();
    mesh->updateInterfaces();
    mesh->update();

    bool check = true;
    if (mesh->getRank() == 0){
        check = (mesh->getNCells() == np) && (mesh->getNVertices() == np);
    }
    MPI_Allreduce(MPI_IN_PLACE, &check, 1, MPI_C_BOOL, MPI_LAND, mesh->getPatch()->getCommunicator());

    if(!check) log<< "Failed to created point cloud"<<std::endl;

    return check;
}
// =================================================================================== //
/*
 * Test: testing creation and writing of a Point Cloud mesh on multiple ranks.
 */
int test2() {

    // ====================================================================== //
    // Initialize the logger
    // ====================================================================== //
    int nProcs;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    bitpit::log::manager().initialize(bitpit::log::COMBINED, true, nProcs, rank);
    bitpit::log::cout().setVisibility(bitpit::log::MASTER);

    int np = 400;
    double rate = double(np)/double(nProcs);
    /* Create target test mesh */
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> target(new mimmo::MimmoObject(3));
    bool check = createMimmoPCMesh(np, bitpit::log::cout(), target.get());

    std::unordered_map<long, int> partMap;
    if(target->getRank()==0){
        livector1D cellIds = target->getCellsIds();
        for(long id : cellIds){
            partMap[id] = std::min(nProcs-1, int(double(id)/rate));
        }
    }

    std::cout<<target->getRank()<<" "<<partMap.size()<<std::endl;

    mimmo::Partition* partition = new mimmo::Partition();
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partition->setName("testParallel00002_PointCloudDistributed");
    partition->setPlotInExecution(true);
    partition->setGeometry(target);
    partition->setPartition(partMap);
    partition->exec();

    if (check){
        bitpit::log::cout() << " test passed " << std::endl;

    }
    delete partition;
    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

    MPI_Init(&argc, &argv);

    /**<Calling mimmo Test routines*/
    int val = 1;
    try{
        val = test2() ;
    }
    catch(std::exception & e){
        std::cout<<"test_parallel_00002 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

    MPI_Finalize();

    return val;
}
