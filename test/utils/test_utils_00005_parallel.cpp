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

#include "mimmo_utils.hpp"
#include "Partition.hpp"
#include <exception>

// =================================================================================== //
/*
 * Test: testing ProjSegmentOnSurface utility in multi procs
 */
int test5() {


    mimmo::MimmoGeometry * reader1 = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::READ);
    bitpit::Logger & log = reader1->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"reading from file "<<std::endl;

    reader1->setReadDir("geodata");
    reader1->setReadFilename("Sphere2Decimated");
    reader1->setReadFileType(FileType::STL);
    reader1->execute();

    log<<"partitioning and distributing geometry "<<std::endl;

    mimmo::Partition * partReader1 = new mimmo::Partition();
    partReader1->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
    partReader1->setGeometry(reader1->getGeometry());
    partReader1->setPlotInExecution(true);
    partReader1->exec();

    log<<"projecting segment "<<std::endl;

    mimmo::ProjSegmentOnSurface * pproj = new mimmo::ProjSegmentOnSurface();
    pproj->setName("test_utils_00005_parallel_ProjectedSegment");
    pproj->setGeometry(partReader1->getGeometry());
    pproj->setSegment({{-0.8,0.6,0.3}}, {{0.8,0.6,0.3}} );
    pproj->setProjElementTargetNCells(320);
    pproj->setPlotInExecution(true);
    pproj->exec();


    int check = 0;
    log<<"test passed "<<std::endl;

    delete reader1;
    delete partReader1;
    delete pproj;

    return check;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif

        int val = 1;

		/**<Calling mimmo Test routines*/
        try{
            val = test5() ;
        }
        catch(std::exception & e){
            std::cout<<"test_utils_00005_parallel exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
