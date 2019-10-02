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

#include "mimmo_propagators.hpp"
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;

// =================================================================================== //

std::unique_ptr<MimmoObject> createTestVolumeMesh(std::vector<long> &bcdir1_vertlist, std::vector<long> &bcdir2_vertlist){

    std::array<double,3> center({{0.0,0.0,0.0}});
    double radiusin(2.0), radiusout(5.0);
    double azimuthin(0.0), azimuthout(0.5*BITPIT_PI);
    double heightbottom(-1.0), heighttop(1.0);
    int nr(6), nt(10), nh(6);

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

    //pump up the vertices
    for(const auto & vertex : verts){
        mesh->addVertex(vertex); //automatic id assigned to vertices.
    }
    mesh->getPatch()->reserveCells(nr*nt*nh);

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
            bcdir1_vertlist.push_back((nr+1)*(nt+1)*k + i);
            bcdir2_vertlist.push_back((nr+1)*(nt+1)*k + (nr+1)*nt + i);
        }
    }
    return mesh;
}


// =================================================================================== //

int test1() {

    std::vector<long> bc1list, bc2list;
    std::unique_ptr<MimmoObject> mesh = createTestVolumeMesh(bc1list, bc2list);

    livector1D cellInterfaceList1 = mesh->getInterfaceFromVertexList(bc1list, true, true);
    livector1D cellInterfaceList2 = mesh->getInterfaceFromVertexList(bc2list, true, true);

    //create the portion of boundary mesh carrying Dirichlet conditions
    std::unique_ptr<MimmoObject> bdirMesh = std::unique_ptr<MimmoObject>(new MimmoObject(1));
    bdirMesh->getPatch()->reserveVertices(bc1list.size()+bc2list.size());
    bdirMesh->getPatch()->reserveCells(cellInterfaceList1.size()+cellInterfaceList2.size());

    for(auto & val : bc1list){
        bdirMesh->addVertex(mesh->getVertexCoords(val), val);
    }
    for(auto & val : bc2list){
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
    bool check = false;
    long targetNode =  (10 +1)*(6+1)*3 + (6+1)*5 + 3;

//TESTING THE SCALAR PROPAGATOR //////
    // and the scalar field of Dirichlet values on its nodes.
    MimmoPiercedVector<double> bc_surf_field;
    bc_surf_field.setGeometry(bdirMesh.get());
    bc_surf_field.setDataLocation(MPVLocation::POINT);
    bc_surf_field.reserve(bdirMesh->getNVertices());
    for(auto & val : bc1list){
        bc_surf_field.insert(val, 10.0);
    }
    for(auto & val : bc2list){
        bc_surf_field.insert(val, 0.0);
    }


    // Now create a PropagateScalarField and solve the laplacian.
    PropagateScalarField * prop = new PropagateScalarField();
    prop->setName("test00001_PropagateScalarField");
    prop->setGeometry(mesh.get());
    prop->setDirichletBoundarySurface(bdirMesh.get());
    prop->setDirichletConditions(bc_surf_field);
    prop->setDumping(false);
    prop->setPlotInExecution(true);

    prop->exec();

    auto values = prop->getPropagatedField();

    check = check || (std::abs(values.at(targetNode)-5.0) > 1.0E-6);


//TESTING THE VECTOR PROPAGATOR //////
    // and the scalar field of Dirichlet values on its nodes.
    MimmoPiercedVector<std::array<double,3>> bc_surf_3Dfield;
    bc_surf_3Dfield.setGeometry(bdirMesh.get());
    bc_surf_3Dfield.setDataLocation(MPVLocation::POINT);
    bc_surf_3Dfield.reserve(bdirMesh->getNVertices());
    for(auto & val : bc1list){
        bc_surf_3Dfield.insert(val, {{10.0, 7.0, -4.0}});
    }
    for(auto & val : bc2list){
        bc_surf_3Dfield.insert(val, {{0.0,0.0,0.0}});
    }

    // Now create a PropagateScalarField and solve the laplacian.
    PropagateVectorField * prop3D = new PropagateVectorField();
    prop3D->setName("test00001_PropagateVectorField");
    prop3D->setGeometry(mesh.get());
    prop3D->setDirichletBoundarySurface(bdirMesh.get());
    prop3D->setDirichletConditions(bc_surf_3Dfield);
    prop3D->setDumping(true);
    prop3D->setDumpingType(1);
    prop3D->setDecayFactor(1.0);
    prop3D->setDumpingInnerDistance(0.5);
    prop3D->setDumpingOuterDistance(3.5);

    prop3D->setPlotInExecution(true);

    prop3D->exec();

    auto values3D = prop3D->getPropagatedField();
    check = check || (norm2(values3D.at(targetNode)-std::array<double,3>({{5.0,3.5,-2.0}})) > 1.0E-6);


    delete prop3D;
    delete prop;

    return check;
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
            val = test1() ;
        }
        catch(std::exception & e){
            std::cout<<"test_propagators_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
