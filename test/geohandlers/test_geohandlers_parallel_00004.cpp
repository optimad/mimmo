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

#include "mimmo_geohandlers.hpp"
#include "Partition.hpp"

// =================================================================================== //
/*!
 * Testing FVSelection classes in parallel.
 */

//creating elementary cube volume test mesh
//nc is the number of cell in x,y,z
mimmo::MimmoSharedPointer<mimmo::MimmoObject> createTestVolumeMesh(const std::array<int,3> & nc ){

    std::array<int,3> np = {{nc[0]+1, nc[1]+1, nc[2]+1}};
    std::array<double,3> dx = {{1./double(nc[0]), 1./double(nc[1]), 1./double(nc[2])}};
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject(2));
    std::vector<std::array<double,3>> points(np[0]*np[1]*np[2]);
    int gentry;
    for(int k=0; k<np[2]; ++k){
        for(int j=0; j<np[1]; ++j){
            for(int i=0; i<np[0]; ++i){
                gentry = np[1]*np[0]*k + np[0]*j + i;
                points[gentry] = std::array<double,3>({{i*dx[0], j*dx[1], k*dx[2]}});
            }
        }
    }
    mesh->getPatch()->reserveVertices(points.size());
    mesh->getPatch()->reserveCells(nc[0]*nc[1]*nc[2]);

    if (mesh->getRank() == 0){
        //push vertices
        long countV(0);
        for(darray3E &pp : points ){
            mesh->addVertex(pp, countV);
            ++countV;
        }

        //push connectivity as HEXA element.
        bitpit::ElementType type = bitpit::ElementType::HEXAHEDRON;
        int rank = -1;
        long PID(1);
        long countC(0);
        livector1D locConn(8);
        for(int k=0; k<nc[2]; ++k){
            for(int j=0; j<nc[1]; ++j){
                for(int i=0; i<nc[0]; ++i){
                    locConn[0] = np[1]*np[0]*k + np[0]*j + i;
                    locConn[1] = np[1]*np[0]*k + np[0]*j + i+1;
                    locConn[2] = np[1]*np[0]*k + np[0]*(j+1) + (i+1);
                    locConn[3] = np[1]*np[0]*k + np[0]*(j+1) + i;
                    locConn[4] = np[1]*np[0]*(k+1) + np[0]*j + i;
                    locConn[5] = np[1]*np[0]*(k+1) + np[0]*j + i+1;
                    locConn[6] = np[1]*np[0]*(k+1) + np[0]*(j+1) + (i+1);
                    locConn[7] = np[1]*np[0]*(k+1) + np[0]*(j+1) + i;

                    mesh->addConnectedCell(locConn, type, PID, countC, rank);
                    ++countC;
                }
            }
        }
    }
    mesh->update();
    return mesh;
}


//version 1 of the test
int testcore1() {

    mimmo::Partition * part = new mimmo::Partition();
    bitpit::Logger & log = part->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    //create bulk mesh on master rank and extract its boundary master rank
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> bulk = createTestVolumeMesh({{20,18,16}});
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> boundary = bulk->extractBoundaryMesh();

    log<<"test v1: created mesh and boundary"<<std::endl;

    //partition it
    part->setName("tg4_v1Partition");
    part->setGeometry(bulk);
    part->setBoundaryGeometry(boundary);
    part->setPartitionMethod(1);
    part->setPlotInExecution(true);
    part->exec();

    log<<"test v1: partitioned"<<std::endl;

    //use sphere selector for bulk+boundary compound.
    mimmo::FVSelectionBySphere * fvsel = new mimmo::FVSelectionBySphere();
    fvsel->setName("tg4_v1FVSelection");
    fvsel->setOrigin({{0.0,0.5,0.5}});
    fvsel->setSpan({{0.5, 2.0*BITPIT_PI, BITPIT_PI}});
    fvsel->setGeometry(part->getGeometry());
    fvsel->setBoundaryGeometry(part->getBoundaryGeometry());
    fvsel->setPlotInExecution(true);
    fvsel->exec();

    log<<"test v1: selected by sphere"<<std::endl;

    log<<"test v1 passed :"<<std::endl;

    delete fvsel;
    delete part;

    return 0;
}


/*
//version 2 of the test
int testcore2() {

    mimmo::Partition * part = new mimmo::Partition();
    bitpit::Logger & log = part->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    //create bulk mesh on master rank and extract its boundary master rank
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> bulk = createTestVolumeMesh({{30,20,16}});
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> boundary = bulk->extractBoundaryMesh();

    log<<"test v2: created mesh and boundary"<<std::endl;

    //use sphere selector for bulk+boundary compound (return selection on master rank).
    mimmo::FVSelectionBySphere * fvsel = new mimmo::FVSelectionBySphere();
    fvsel->setName("tg4_v2FVSelection");
    fvsel->setOrigin({{0.0,0.5,0.5}});
    fvsel->setSpan({{0.5, 2.0*BITPIT_PI, BITPIT_PI}});
    fvsel->setGeometry(bulk);
    fvsel->setBoundaryGeometry(boundary);
    fvsel->setPlotInExecution(true);
    fvsel->exec();

    log<<"test v2: selected mesh and boundaries"<<std::endl;

    //partition bulk and internal boundaries
    part->setName("tg4_v2Partition");
    part->setGeometry(fvsel->getVolumePatch());
    part->setBoundaryGeometry(fvsel->getInternalBoundaryPatch());
    part->setPartitionMethod(1);
    part->setPlotInExecution(false);
    part->exec();

    log<<"test v2: partition of bulk selection and internal boundaries"<<std::endl;

    part->getGeometry()->getPatch()->write("tg4_v2FVSelection_Volume_Patch_Partitioned");
    part->getBoundaryGeometry()->getPatch()->write("tg4_v2FVSelection_InternalBoundary_Patch_Partitioned");

    log<<"test v2 passed :"<<std::endl;

    delete fvsel;
    delete part;

    return 0;
}
*/
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
        val = testcore1();
        //val = std::max(testcore1(), val);
	}
	catch(std::exception & e){
		std::cout<<"test_geohandlers_parallel_00004 exited with an error of type : "<<e.what()<<std::endl;
		return 1;
	}
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
