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
#include "mimmo_parallel.hpp"
#include "mimmo_iogeneric.hpp"
#include <iostream>
#include <chrono>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

// =================================================================================== //
/*
	example parallel_example_00003.cpp

	brief Example of usage of Serialization of a partioned geometry.

	Block used: Partition and Serialization.

	<b>To run</b>: mpirun -np X parallel_example_00003 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //
//creating a mesh on master rank 0 and fill on it lists of marked boundary vertices
mimmo::MimmoSharedPointer<mimmo::MimmoObject> createTestVolumeMesh(int rank, std::vector<long> &bcdir1_vertlist, std::vector<long> &bcdir2_vertlist){

	std::array<double,3> center({{0.0,0.0,0.0}});
	double radiusin(2.0), radiusout(5.0);
	double azimuthin(0.0), azimuthout(0.5*BITPIT_PI);
	double heightbottom(-1.0), heighttop(1.0);
	int nr(20), nt(20), nh(20);

	double deltar = (radiusout - radiusin)/ double(nr);
	double deltat = (azimuthout - azimuthin)/ double(nt);
	double deltah = (heighttop - heightbottom)/ double(nh);

	std::vector<std::array<double,3> > verts ((nr+1)*(nt+1)*(nh+1));

	int counter = 0;
	for(int k=0; k<=nh; ++k){
		for(int j=0; j<=nt; ++j){
			for(int i=0; i<=nr; ++i){
				verts[counter][0] =(radiusin + i*deltar)*std::cos(azimuthin + j*deltat);
				verts[counter][1] =(radiusin + i*deltar)*std::sin(azimuthin + j*deltat);
				verts[counter][2] =(heightbottom + k*deltah);
				++counter;
			}
		}
	}

	//create the volume mesh mimmo.
	mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh(new mimmo::MimmoObject(2));

	if (rank == 0){

		mesh->getPatch()->reserveVertices((nr+1)*(nt+1)*(nh+1));
		mesh->getPatch()->reserveCells(nr*nt*nh);

		//pump up the vertices
		for(const auto & vertex : verts){
			mesh->addVertex(vertex); //automatic id assigned to vertices.
		}

		//create connectivities for hexa elements
		std::vector<long> conn(8,0);
		for(int k=0; k<nh; ++k){
			for(int j=0; j<nt; ++j){
				for(int i=0; i<nr; ++i){
					conn[0] = (nr+1)*(nt+1)*k + (nr+1)*j + i;
					conn[1] = (nr+1)*(nt+1)*k + (nr+1)*j + i+1;
					conn[2] = (nr+1)*(nt+1)*k + (nr+1)*(j+1) + i+1;
					conn[3] = (nr+1)*(nt+1)*k + (nr+1)*(j+1) + i;
					conn[4] = (nr+1)*(nt+1)*(k+1) + (nr+1)*j + i;
					conn[5] = (nr+1)*(nt+1)*(k+1) + (nr+1)*j + i+1;
					conn[6] = (nr+1)*(nt+1)*(k+1) + (nr+1)*(j+1) + i+1;
					conn[7] = (nr+1)*(nt+1)*(k+1) + (nr+1)*(j+1) + i;

					mesh->addConnectedCell(conn, bitpit::ElementType::HEXAHEDRON);
				}
			}
		}

		bcdir1_vertlist.clear();
		bcdir2_vertlist.clear();
		bcdir1_vertlist.reserve(mesh->getNVertices());
		bcdir2_vertlist.reserve(mesh->getNVertices());

		for(int k=0; k<=nh; ++k){
			for(int i=0; i<=nr; ++i){
				bcdir1_vertlist.push_back((nr+1)*(nt+1)*k + i);
				bcdir2_vertlist.push_back((nr+1)*(nt+1)*k + (nr+1)*nt + i);
			}
		}
	}

	mesh->updateInterfaces();
	mesh->update();

	return mesh;
}

// =================================================================================== //

int test00003() {

	// Initialize mpi
	int nProcs;
	int    rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<long> bc1list, bc2list;
	mimmo::MimmoSharedPointer<mimmo::MimmoObject> mesh = createTestVolumeMesh(rank, bc1list, bc2list);
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> boundary= mesh->extractBoundaryMesh();
    //pidding the boundary
    livector1D pidlist1 = boundary->getCellFromVertexList(bc1list, true);
    livector1D pidlist2 = boundary->getCellFromVertexList(bc2list, true);
    //pid the cells.
    for(long cellid : pidlist1){
        boundary->setPIDCell(cellid, long(1));
    }
    for(long cellid : pidlist2){
        boundary->setPIDCell(cellid, long(2));
    }

    //create a subselection of the boundary
    mimmo::MimmoSharedPointer<mimmo::MimmoObject> bDirMesh(new mimmo::MimmoObject(boundary->getType()));

    if(rank == 0){   //create a subselection of the whole boundary
        //vertices
        std::unordered_set<long> potvertices(bc1list.begin(), bc1list.end());
        potvertices.insert(bc2list.begin(), bc2list.end());
        for(long idv: potvertices){
            bDirMesh->addVertex(boundary->getVertexCoords(idv), idv);
        }
        //cell
        //mark to delete all cells not in the pidlist1 and 2.
        std::unordered_set<long> preservedcells (pidlist1.begin(), pidlist1.end());
        preservedcells.insert(pidlist2.begin(), pidlist2.end());

        for(long idc: preservedcells){
            bDirMesh->addCell(boundary->getPatch()->getCell(idc), idc, rank);
        }
    }
    bDirMesh->updateInterfaces();
    bDirMesh->update();

	/* Instantiation of a Partition object with default patition method space filling curve.
	 * Plot Optional results during execution active for Partition block.
	 */
	mimmo::Partition* partition = new mimmo::Partition();
	partition->setPlotInExecution(true);
	partition->setGeometry(mesh);
	partition->setBoundaryGeometry(bDirMesh);
    partition->setPartitionMethod(mimmo::PartitionMethod::PARTGEOM);
	auto t1 = Clock::now();
	if (rank == 0)
		std::cout << "#" << rank  << " Start Partition mesh " << std::endl;
	partition->exec();
	auto t2 = Clock::now();
	if (rank == 0){
		std::cout << "#" << rank << " Partition mesh execution time: "
				<< std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
				<< " seconds" << std::endl;
	}
	/* Creation of mimmo containers. MimmoGeometry used to dump partitioned mesh
	 */
	mimmo::MimmoGeometry * mimmoVolumeOut = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::WRITE);
	mimmoVolumeOut->setWriteDir("./");
	mimmoVolumeOut->setWriteFileType(FileType::MIMMO);
	mimmoVolumeOut->setWriteFilename("parallel_example_00003.volume.partitioned");
	mimmoVolumeOut->setGeometry(mesh);
	mimmoVolumeOut->exec();

	/* Creation of mimmo containers. MimmoGeometry used to restore partitioned mesh
	 */
	mimmo::MimmoGeometry * mimmoVolumeIn = new mimmo::MimmoGeometry(mimmo::MimmoGeometry::IOMode::CONVERT);
	mimmoVolumeIn->setReadDir("./");
	mimmoVolumeIn->setReadFileType(FileType::MIMMO);
	mimmoVolumeIn->setReadFilename("parallel_example_00003.volume.partitioned");
	mimmoVolumeIn->setWriteDir("./");
	mimmoVolumeIn->setWriteFileType(FileType::VOLVTU);
	mimmoVolumeIn->setWriteFilename("parallel_example_00003.volume.restored");
	mimmoVolumeIn->exec();

	/* Instantiation of a Partition object with serialize partition method.
	 * Plot Optional results during execution active for Partition block.
	 */
	mimmo::Partition* serialize = new mimmo::Partition();
	serialize->setName("mimmo.Serialization");
	serialize->setPlotInExecution(true);
	serialize->setGeometry(mesh);//mimmoVolumeIn->getGeometry());
	serialize->setBoundaryGeometry(bDirMesh);
	serialize->setPartitionMethod(mimmo::PartitionMethod::SERIALIZE);

	t1 = Clock::now();
	if (rank == 0)
		std::cout << "Start Serialize mesh " << std::endl;
	serialize->exec();
	t2 = Clock::now();
	if (rank == 0){
		std::cout << "Serialize mesh execution time: "
		<< std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
		<< " seconds" << std::endl;
	}

	bool error = false;

	delete partition;
	delete serialize;
	return error;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
	/**<Calling mimmo Test routines*/

	int val = test00003() ;

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
