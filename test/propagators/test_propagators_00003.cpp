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

std::unique_ptr<MimmoObject> createTestVolumeMesh(std::unique_ptr<MimmoObject> & boundary){

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

    //decompose manually the six boundary faces.
    std::vector<std::vector<long>> bverts(6, std::vector<long>());

    //vertex at azimuthal limits.
    for(int k=0; k<=nh; ++k){
        for(int i=0; i<=nr; ++i){
            bverts[0].push_back((nr+1)*(nt+1)*k + i);
            bverts[1].push_back((nr+1)*(nt+1)*k + (nr+1)*nt + i);
        }
    }

    //vertex at quote limits.
    for(int j=0; j<=nt; ++j){
        for(int i=0; i<=nr; ++i){
            bverts[2].push_back((nr+1)*j + i);
            bverts[3].push_back((nr+1)*(nt+1)*nh + (nr+1)*j + i);
        }
    }

    //vertex at radius limits.
    for(int k=0; k<=nh; ++k){
        for(int j=0; j<=nt; ++j){
            bverts[4].push_back((nr+1)*(nt+1)*k + (nr+1)*j);
            bverts[5].push_back((nr+1)*(nt+1)*k + (nr+1)*j + nr);
        }
    }

    //put together
    MimmoPiercedVector<long> pidfaces;
    long pidder = 1;
    for(livector1D & vert : bverts ){
        livector1D boundaryfaces = mesh->getInterfaceFromVertexList(vert, true, true);
        for(long id : boundaryfaces){
            pidfaces.insert(id, pidder);
        }
        ++pidder;
    }
    bverts.clear();
    std::vector<long> boundaryverts = mesh->extractBoundaryVertexID();

    //create the boundary meshes
    auto orinterfaces = mesh->getInterfaces();
    auto orcells = mesh->getCells();

    boundary = std::unique_ptr<MimmoObject>(new MimmoObject(1));
    boundary->getPatch()->reserveVertices(boundaryverts.size());
    boundary->getPatch()->reserveCells(pidfaces.size());

    //push in verts
    for(long id : boundaryverts){
        boundary->addVertex(mesh->getVertexCoords(id), id);
    }

    //push in interfaces ad 2D cells
    for(auto it=pidfaces.begin(); it!=pidfaces.end(); ++it){
        bitpit::Interface & ii = orinterfaces.at(it.getId());
        long * conn = ii.getConnect();
        std::size_t connsize = ii.getConnectSize();

        bitpit::ElementType et = orcells.at(ii.getOwner()).getFaceType(ii.getOwnerFace());

        boundary->addConnectedCell(std::vector<long>(&conn[0], &conn[connsize]), et, *it,it.getId());
    }
    boundary->buildAdjacencies();
    boundary->getPatch()->write("piddedBoundary_test3");

    return mesh;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
std::unique_ptr<MimmoObject> pidExtractor(MimmoObject* boundary, const livector1D & pids){
    if(!boundary) return nullptr;

    livector1D cellsID = boundary->extractPIDCells(pids);
    livector1D vertID = boundary->getVertexFromCellList(cellsID);

    std::unique_ptr<MimmoObject> extracted(new MimmoObject(1));
    extracted->getPatch()->reserveVertices(vertID.size());
    extracted->getPatch()->reserveCells(cellsID.size());

    for(long id: vertID){
        extracted->addVertex(boundary->getVertexCoords(id), id);
    }
    for(long id: cellsID){
        extracted->addCell(boundary->getCells().at(id), id);
    }
    extracted->buildAdjacencies();
    return extracted;
}

// =================================================================================== //

int test1() {

    std::unique_ptr<MimmoObject> bDirMesh;
    std::unique_ptr<MimmoObject> mesh = createTestVolumeMesh(bDirMesh);


    std::unique_ptr<MimmoObject> bDIR = pidExtractor(bDirMesh.get(), {{1,2}});
    std::unique_ptr<MimmoObject> bSLIP = pidExtractor(bDirMesh.get(), {{5,6}});

//TESTING THE VECTOR PROPAGATOR //////
    // and the scalar field of Dirichlet values on pidded surface 0 and 1 nodes.
    MimmoPiercedVector<std::array<double,3>> bc_surf_3Dfield;
    bc_surf_3Dfield.setGeometry(bDIR.get());
    bc_surf_3Dfield.setDataLocation(MPVLocation::POINT);
    bc_surf_3Dfield.reserve(bDIR->getNVertices());

    std::vector<long> d1 = bDIR->getVertexFromCellList(bDIR->extractPIDCells(1));
    for(long id : d1){
        bc_surf_3Dfield.insert(id, {{0.0, 0.0, 0.0}});
    }
    //apply a rotation of alpha around z centered in the origin 0,0,0, so that deformation is
    //x-directed and equal to tan(alpha) * y-coord of the point.
    //BEWARE IF THE ROTATION ARC IS similar to a linear segment, this deformation works,
    //for large rotation, where the arc path is very different from a different segment, this
    // could lead to not properly expected results.
    double alpha = 9.0*BITPIT_PI/180.0;
    std::vector<long> d2 = bDIR->getVertexFromCellList(bDIR->extractPIDCells(2));
    for(long id : d2){
        darray3E coord = bDIR->getVertexCoords(id);
        double normc = std::sqrt(coord[0]*coord[0] + coord[1]*coord[1]);
        bc_surf_3Dfield.insert(id, {{std::sin(alpha)*normc,
                                     (std::cos(alpha)-1.0)*normc,
                                      0.0}}  );
    }

    // Now create a PropagateScalarField and solve the laplacian.
    PropagateVectorField * prop3D = new PropagateVectorField();
    prop3D->setName("test00003_PropagateVectorField");
    prop3D->setGeometry(mesh.get());
    prop3D->setDirichletBoundarySurface(bDIR.get());
    prop3D->setDirichletConditions(bc_surf_3Dfield);
    prop3D->setSlipBoundarySurface(bSLIP.get());

    prop3D->setMethod(PropagatorMethod::GRAPHLAPLACE);

    prop3D->setDumping(true);
    prop3D->setDumpingType(0);
    prop3D->setDecayFactor(2.0);
    prop3D->setDumpingInnerDistance(0.5);
    prop3D->setDumpingOuterDistance(2.5);

    prop3D->setSolverMultiStep(3);

    prop3D->setPlotInExecution(true);

    prop3D->exec();

    //deform the mesh and write deformation
    dmpvecarr3E values3D = prop3D->getPropagatedField();
    darray3E work;
    for(auto it=values3D.begin(); it!=values3D.end(); ++it){
        work = mesh->getVertexCoords(it.getId());
        mesh->modifyVertex(work + *it, it.getId());
    }

    mesh->getPatch()->write("test00003_deformedMesh");


    bool check = false;
    long targetNode =  (10 +1)*(6+1)*3 + (6+1)*5 + 3;

    std::cout<<values3D.at(targetNode)<<std::endl;

    check = check || (norm2(values3D.at(targetNode)-std::array<double,3>({{0.210569, -0.0519423, -1.40474e-18}})) > 1.0E-5);


    delete prop3D;
    return check;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routines*/

        int val = test1() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return val;
}
