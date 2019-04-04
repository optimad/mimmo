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
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;

#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

// =================================================================================== //

std::unique_ptr<MimmoObject> createTestVolumeMesh( std::vector<bitpit::Vertex> &bcdir1_vertlist, std::vector<bitpit::Vertex> &bcdir2_vertlist){

    std::array<double,3> center({{0.0,0.0,0.0}});
    double radiusin(2.0), radiusout(5.0);
    double azimuthin(0.0), azimuthout(0.5*BITPIT_PI);
    double heightbottom(-1.0), heighttop(1.0);
    int nr(50), nt(100), nh(50);

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
    std::unique_ptr<MimmoObject> mesh = std::unique_ptr<MimmoObject>(new MimmoObject(2));
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
    mesh->buildAdjacencies();
    mesh->buildInterfaces();

    bcdir1_vertlist.clear();
    bcdir2_vertlist.clear();
    bcdir1_vertlist.reserve(mesh->getNVertices());
    bcdir2_vertlist.reserve(mesh->getNVertices());

    for(int k=0; k<=nh; ++k){
        for(int i=0; i<=nr; ++i){
            bcdir1_vertlist.push_back(mesh->getPatch()->getVertex((nr+1)*(nt+1)*k + i));
            bcdir2_vertlist.push_back(mesh->getPatch()->getVertex((nr+1)*(nt+1)*k + (nr+1)*nt + i));
        }
    }
    return mesh;
}


// =================================================================================== //

int test00003() {

    std::vector<bitpit::Vertex> bc1list, bc2list;
    std::unique_ptr<MimmoObject> mesh = createTestVolumeMesh(bc1list, bc2list);

    std::vector<long> bc1list_, bc2list_;
    for (auto v : bc1list)
    	bc1list_.push_back(v.getId());
    for (auto v : bc2list)
    	bc2list_.push_back(v.getId());

    livector1D cellInterfaceList1 = mesh->getInterfaceFromVertexList(bc1list_, true, true);
    livector1D cellInterfaceList2 = mesh->getInterfaceFromVertexList(bc2list_, true, true);

    //create the portion of boundary mesh carrying Dirichlet conditions
    std::unique_ptr<MimmoObject> bdirMesh = std::unique_ptr<MimmoObject>(new MimmoObject(1));
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
    }
    for(auto & val : cellInterfaceList2){
        int sizeconn =mesh->getInterfaces().at(val).getConnectSize();
        long * conn = mesh->getInterfaces().at(val).getConnect();
        bdirMesh->addConnectedCell(std::vector<long>(&conn[0], &conn[sizeconn]),
                                    bitpit::ElementType::QUAD, val);
    }

    bdirMesh->buildAdjacencies();


    /* Instantiation of a Partition object with default patition method space filling curve.
     * Plot Optional results during execution active for Partition block.
     */
    Partition* partition = new Partition();
    partition->setPlotInExecution(true);
    partition->setGeometry(mesh.get());
    partition->setBoundaryGeometry(bdirMesh.get());
    auto t1 = Clock::now();
    std::cout << "Start Partition mesh " << std::endl;
    partition->exec();
    auto t2 = Clock::now();
    std::cout << "Partition mesh execution time: "
              << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
              << " seconds" << std::endl;

    // and the field of Dirichlet values on its nodes.
    MimmoPiercedVector<std::array<double,3>> bc_surf_field;
    bc_surf_field.setGeometry(partition->getBoundaryGeometry());
    bc_surf_field.setDataLocation(MPVLocation::POINT);
    bc_surf_field.reserve(partition->getBoundaryGeometry()->getNVertices());

    partition->getBoundaryGeometry()->buildKdTree();
    bitpit::KdTree<3, bitpit::Vertex, long> *  tree = partition->getBoundaryGeometry()->getKdTree();

    for(auto & val : bc1list){
    	int node = tree->exist(&val);
    	if (node > -1){
    		bc_surf_field.insert(tree->nodes[node].object_->getId(), {{10., 5., 1.}});
    	}
    }
    for(auto & val : bc2list){
    	int node = tree->exist(&val);
    	if (node > -1){
    		bc_surf_field.insert(tree->nodes[node].object_->getId(), {{0., 0., 0.}});
    	}
    }

    // Now create a PropagateScalarField and solve the laplacian.
    PropagateVectorField * prop = new PropagateVectorField();
    prop->setGeometry(mesh.get());
    prop->setDirichletBoundarySurface(partition->getBoundaryGeometry());
    prop->setDirichletConditions(bc_surf_field);
    prop->setDumping(false);
    prop->setPlotInExecution(true);

    t1 = Clock::now();
    std::cout << "Start Propagator vector field " << std::endl;
    prop->exec();
    t2 = Clock::now();
    std::cout << "Propagator vector field execution time: "
              << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()
              << " seconds" << std::endl;


//    auto values = prop->getPropagatedField();
//    long targetNode = (10 +1)*(6+1)*3 + (6+1)*5 + 3;
//    bool check = std::abs(values.at(targetNode)-5.0) > 1.0E-6;

    bool error = false;

    delete partition;
    delete prop;
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
