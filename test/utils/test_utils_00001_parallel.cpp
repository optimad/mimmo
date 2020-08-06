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

#include "mimmo_utils.hpp"
#include "Partition.hpp"
#include <exception>
#include <random>

// =================================================================================== //
/*
 * Test: mirroring a point with data attached w.r.t a plane
 */
int test1() {

    mimmo::Partition * part = new mimmo::Partition();
    bitpit::Logger & log =part->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    //create the mesh
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> pc(new mimmo::MimmoObject(3));
    long np = 20;
    int rank = pc->getRank();
    int nprocs = pc->getProcessorCount();
    double rate = double(np)/double(nprocs);
    std::unordered_map<long, int> partitionMap;

    if(rank == 0){
        std::random_device rd;
        unsigned int seed = rd(); //fix it up to any unsigned int for reproducibility issue.
        std::mt19937_64 rgen(seed);
        std::uniform_real_distribution<double> distr(0.0, 1.0);
        for(long i=0; i<np; ++i){
            pc->addVertex({{distr(rgen),1.,distr(rgen)}}, i);
            pc->addConnectedCell(std::vector<long>(1, i), bitpit::ElementType::VERTEX, 0, i, rank);
            partitionMap[i] = std::min(10, int(double(i)/rate));
        }
    }

    pc->buildAdjacencies();
    pc->setPartitioned();

    //fill and execute partition
    part->setGeometry(pc);
    part->setPartition(partitionMap);
    part->exec();

    //create the data
    mimmo::MimmoPiercedVector<darray3E> data;
    data.initialize(part->getGeometry(),mimmo::MPVLocation::POINT, {{1,2,3}});

//    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"Mirroring points in parallel..."<<std::endl;

    mimmo::SpecularPoints * sp = new mimmo::SpecularPoints();
    sp->setName("test_00001_SpecularPoints_parallel");
    sp->setPointCloud(part->getGeometry());
    sp->setVectorData(&data);
    sp->setPlane({{0,0,0}},{{0,1,0}});
    sp->setPlotInExecution(true);
    sp->exec();

//    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"Mirroring points in parallel... done"<<std::endl;

    bool check = (sp->getMirroredPointCloud()->getNGlobalCells() == 2*np );
    if(!check)  {
        log<<"test failed "<<std::endl;
    }else{
        log<<"test passed "<<std::endl;
    }

    delete sp;
    delete part;
    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
		/**<Calling mimmo Test routines*/
        int val = 1;
        try{
            val = test1();
        }

        catch(std::exception & e){
            std::cout<<"test_utils_00001_parallel exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
