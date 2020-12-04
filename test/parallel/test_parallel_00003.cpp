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
 #include "Partition.hpp"

// =================================================================================== //
/*!
 * Creating a partitioned boundary volume and extract the boundary from MimmoObject in
   distributed mode
 */

//creating elementary cube volume test mesh on master rank 0
//nc is the number of cell in x,y,z
mimmo::MimmoSharedPointer<mimmo::MimmoObject> createTestVolumeMesh(const std::array<int,3> & nc ){

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject(2));
    if(mesh->getRank() == 0){

        std::array<int,3> np = {{nc[0]+1, nc[1]+1, nc[2]+1}};
        std::array<double,3> dx = {{1./double(nc[0]), 1./double(nc[1]), 1./double(nc[2])}};

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


//proper core of the test
int testcore() {

    //create bulk mesh and extract its boundary
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> bulk = createTestVolumeMesh({{10,9,8}});

    mimmo::Partition * part = new mimmo::Partition();
    part->setGeometry(bulk);
    part->setPartitionMethod(1);
    part->exec();

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> boundary = part->getGeometry()->extractBoundaryMesh();

    part->getGeometry()->getPatch()->write("tp00003_originalBulk");
    boundary->getPatch()->write("tp00003_originalBoundary");


    bitpit::Logger & log = boundary->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);
    log<<"test 1 passed :"<<std::endl;

    delete part;

    return 0;
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
		val = testcore();
	}
	catch(std::exception & e){
		std::cout<<"test_parallel_00003 exited with an error of type : "<<e.what()<<std::endl;
		return 1;
	}
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
