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
#include "mimmo_propagators.hpp"
#include "mimmo_iogeneric.hpp"
#include <iostream>
#include <chrono>
#include <fstream>

typedef std::chrono::high_resolution_clock Clock;

// =================================================================================== //
/*
	example parallel_example_00004.cpp

	brief Example of usage of Serialization of a partioned geometry.

	Block used: Partition and Serialization.

	<b>To run</b>: mpirun -np X parallel_example_00004 \n

	<b> visit</b>: <a href="http://optimad.github.io/mimmo/">mimmo website</a> \n
 */


// =================================================================================== //

std::unique_ptr<mimmo::MimmoObject> createTestVolumeMesh(int rank, std::vector<bitpit::Vertex> &bcdir1_vertlist, std::vector<bitpit::Vertex> &bcdir2_vertlist){

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
	std::unique_ptr<mimmo::MimmoObject> mesh = std::unique_ptr<mimmo::MimmoObject>(new mimmo::MimmoObject(2));

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
		std::vector<long> list1, list2;
		list1.reserve(mesh->getNVertices());
		list2.reserve(mesh->getNVertices());

		for(int k=0; k<=nh; ++k){
			for(int i=0; i<=nr; ++i){
				bcdir1_vertlist.push_back(mesh->getPatch()->getVertex((nr+1)*(nt+1)*k + i));
				bcdir2_vertlist.push_back(mesh->getPatch()->getVertex((nr+1)*(nt+1)*k + (nr+1)*nt + i));
				list1.push_back((nr+1)*(nt+1)*k + i);
				list2.push_back((nr+1)*(nt+1)*k + (nr+1)*nt + i);
			}
		}
	}

	mesh->buildAdjacencies();
	mesh->buildInterfaces();

	return mesh;
}

// =================================================================================== //

int test00004() {

	// Initialize mpi
	int nProcs;
	int    rank;
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::vector<bitpit::Vertex> bc1list, bc2list;
	std::unique_ptr<mimmo::MimmoObject> mesh = createTestVolumeMesh(rank, bc1list, bc2list);

	std::vector<long> bc1list_, bc2list_;
	for (auto v : bc1list)
		bc1list_.push_back(v.getId());
	for (auto v : bc2list)
		bc2list_.push_back(v.getId());

	livector1D cellInterfaceList1 = mesh->getInterfaceFromVertexList(bc1list_, true, true);
	livector1D cellInterfaceList2 = mesh->getInterfaceFromVertexList(bc2list_, true, true);

	//create the portion of boundary mesh carrying Dirichlet conditions
	std::unique_ptr<mimmo::MimmoObject> bdirMesh = std::unique_ptr<mimmo::MimmoObject>(new mimmo::MimmoObject(1));
	if (rank == 0){
		bdirMesh->getPatch()->reserveVertices(bc1list.size()+bc2list.size());
		bdirMesh->getPatch()->reserveCells(cellInterfaceList1.size()+cellInterfaceList2.size());

		for(auto & val : bc1list_){
			bdirMesh->addVertex(mesh->getVertexCoords(val), val);
		}
		for(auto & val : bc2list_){
			bdirMesh->addVertex(mesh->getVertexCoords(val), val);
		}
		for(auto & val : cellInterfaceList1){
			int sizeconn =mesh->getInterfaces().at(val).getConnectSize();
			long * conn = mesh->getInterfaces().at(val).getConnect();
			bdirMesh->addConnectedCell(std::vector<long>(&conn[0], &conn[sizeconn]),
					bitpit::ElementType::QUAD, val);
			bdirMesh->getPatch()->getCell(val).setPID(1);
		}
		for(auto & val : cellInterfaceList2){
			int sizeconn =mesh->getInterfaces().at(val).getConnectSize();
			long * conn = mesh->getInterfaces().at(val).getConnect();
			bdirMesh->addConnectedCell(std::vector<long>(&conn[0], &conn[sizeconn]),
					bitpit::ElementType::QUAD, val);
			bdirMesh->getPatch()->getCell(val).setPID(2);
		}
	}
	bdirMesh->buildAdjacencies();
	bdirMesh->buildInterfaces();

	/* Instantiation of a Partition object with default patition method space filling curve.
	 * Plot Optional results during execution active for Partition block.
	 */
	mimmo::Partition* partition = new mimmo::Partition();
	partition->setPlotInExecution(true);
	partition->setGeometry(mesh.get());
	partition->setBoundaryGeometry(bdirMesh.get());
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
	mimmo::MimmoGeometry * mimmoVolumeOut = new mimmo::MimmoGeometry();
	mimmoVolumeOut->setIOMode(IOMode::WRITE);
	mimmoVolumeOut->setWriteDir("./");
	mimmoVolumeOut->setWriteFileType(FileType::MIMMO);
	mimmoVolumeOut->setWriteFilename("parallel_example_00004.volume.partitioned");
	mimmoVolumeOut->setGeometry(mesh.get());
	mimmoVolumeOut->exec();

	/* Creation of mimmo containers. MimmoGeometry used to restore partitioned mesh
	 */
	mimmo::MimmoGeometry * mimmoVolumeIn = new mimmo::MimmoGeometry();
	mimmoVolumeIn->setIOMode(IOMode::CONVERT);
	mimmoVolumeIn->setReadDir("./");
	mimmoVolumeIn->setReadFileType(FileType::MIMMO);
	mimmoVolumeIn->setReadFilename("parallel_example_00004.volume.partitioned");
	mimmoVolumeIn->setWriteDir("./");
	mimmoVolumeIn->setWriteFileType(FileType::VOLVTU);
	mimmoVolumeIn->setWriteFilename("parallel_example_00004.volume.restored");
	mimmoVolumeIn->exec();

	/* Instantiation of a Partition object with serialize partition method.
	 * Plot Optional results during execution active for Partition block.
	 */
	mimmo::Partition* serialize = new mimmo::Partition();
	serialize->setName("mimmo.Serialization");
	serialize->setPlotInExecution(true);
	serialize->setGeometry(mimmoVolumeIn->getGeometry());
	serialize->setBoundaryGeometry(bdirMesh.get());
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

	int val = test00004() ;

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
