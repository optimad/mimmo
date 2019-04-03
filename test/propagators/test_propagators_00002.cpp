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
/*!
 * Create bar grid with prismatic layer + tetrahedral bulk and an intermediate ref level between them
 * \param[out]boundary layer surface grid pidded in 3, 0 front face prismatic, 1 back bulk bar faces, 2 other.
 * \return the volume grid
 */
std::unique_ptr<MimmoObject> createTestVolumeMesh( std::unique_ptr<MimmoObject> &boundary)
{

    std::array<double, 3> origin({{0.0,0.0,0.0}});
    double width(0.2), height(0.2);
    double prismdepth(0.4), bulkdepth(0.6);

    int nw(10), nh(10), npd(20), nbd(30);

    double deltaw = width/ double(nw);
    double deltah = height/ double(nh);
    double deltapd = prismdepth/ double(npd);
    double deltabd = bulkdepth/ double(nbd);

    std::vector<std::array<double,3> > verts ((nw+1)*(nh+1)*(npd+nbd+1));

    int counter = 0;
    //prism layer verts
    for(int k=0; k<npd; ++k){
        for(int j=0; j<=nh; ++j){
            for(int i=0; i<=nw; ++i){
                verts[counter][0] = origin[0] + deltaw  * i;
                verts[counter][1] = origin[1] + deltah  * j;
                verts[counter][2] = origin[2] + deltapd * k;
                ++counter;
            }
        }
    }
    //bulk layer verts
    for(int k=0; k<=nbd; ++k){
        for(int j=0; j<=nh; ++j){
            for(int i=0; i<=nw; ++i){
                verts[counter][0] = origin[0] + deltaw  * i;
                verts[counter][1] = origin[1] + deltah  * j;
                verts[counter][2] = origin[2] + deltapd * npd + deltabd * k;
                ++counter;
            }
        }
    }

    //create the volume mesh mimmo.
    std::unique_ptr<MimmoObject> mesh = std::unique_ptr<MimmoObject>(new MimmoObject(2));
    mesh->getPatch()->reserveVertices(counter);

    //pump up the vertices
    for(const auto & vertex : verts){
        mesh->addVertex(vertex); //automatic id assigned to vertices.
    }
    //using Kuhn decomposition in tetrahedra from hexahedron, cut on diagonal for prism layer cells from hexahedron
    mesh->getPatch()->reserveCells(2*nw*nh*npd + 6*nw*nh*nbd);

    //create connectivities for wedge elements  of prism layer
    std::vector<long> conn(6,0);
    for(int k=0; k<npd; ++k){
        for(int j=0; j<nh; ++j){
            for(int i=0; i<nw; ++i){
                //FIRST 0-2-1-4-6-7 from elemental hexa
                conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                conn[1] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i+1;
                conn[2] = (nw+1)*(nh+1)*k + (nw+1)*j + i+1;
                conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i;
                conn[4] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;
                conn[5] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i+1;

                mesh->addConnectedCell(conn, bitpit::ElementType::WEDGE);
                //SECOND 0-3-2-4-7-6 from elemental hexa
                conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                conn[1] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i;
                conn[2] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i+1;
                conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i;
                conn[4] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i;
                conn[5] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                mesh->addConnectedCell(conn, bitpit::ElementType::WEDGE);
            }
        }
    }

    //create connectivities for tetra elements: USING KUHN
    //BEWARE, the coarser the tetrahedral kuhn adapted part the worse the non-orthogonality of cells.
    // This impacts on LAPLACIAL SOLUTION. //TODO Use a better tetrahedral element generation.

    conn.clear();
    conn.resize(4,0);
    for(int k=npd; k<(npd+nbd); ++k){
        for(int j=0; j<nh; ++j){
            for(int i=0; i<nw; ++i){
                    //FIRST 0-1-2-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*k + (nw+1)*j + i+1;
                    conn[2] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i+1;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                    //SECOND 0-5-1-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i+1;
                    conn[2] = (nw+1)*(nh+1)*k + (nw+1)*j + i+1;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                    //THIRD 0-2-3-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i+1;
                    conn[2] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                    //FIRST 0-4-5-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i;
                    conn[2] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i+1;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                    //SECOND 0-7-4-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i;
                    conn[2] = (nw+1)*(nh+1)*(k+1) + (nw+1)*j + i;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                    //THIRD 0-3-7-6 from elemenhal hexa
                    conn[0] = (nw+1)*(nh+1)*k + (nw+1)*j + i;
                    conn[1] = (nw+1)*(nh+1)*k + (nw+1)*(j+1) + i;
                    conn[2] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i;
                    conn[3] = (nw+1)*(nh+1)*(k+1) + (nw+1)*(j+1) + i+1;

                    mesh->addConnectedCell(conn, bitpit::ElementType::TETRA);

                }
            }
        }

    mesh->buildAdjacencies();
    mesh->buildInterfaces();

    livector1D borderverts = mesh->extractBoundaryVertexID();

    //get the subset of ids nearest at z=0 and z= prismdepth+bulkdepth;
    double zin(0.0), zout(prismdepth+bulkdepth), zwork;
    std::vector<long> boundary1verts, boundary2verts;
    for(long idv : borderverts){
        zwork = mesh->getVertexCoords(idv)[2];
        if(std::abs(zwork - zin) < 0.5*deltapd){
            boundary1verts.push_back(idv);
        }
        if(std::abs(zwork - zout) < 0.5*deltabd){
            boundary2verts.push_back(idv);
        }
    }
    livector1D boundary1faces = mesh->getInterfaceFromVertexList(boundary1verts, true, true);
    livector1D boundary2faces = mesh->getInterfaceFromVertexList(boundary2verts, true, true);

    //put together
    MimmoPiercedVector<long> pidfaces;
    pidfaces.reserve(boundary1faces.size()+boundary2faces.size());
    for(long id : boundary1faces){
        pidfaces.insert(id,0);
    }
    for(long id : boundary2faces){
        pidfaces.insert(id,1);
    }

    std::vector<long> boundaryverts(boundary1verts.begin(), boundary1verts.end());
    boundaryverts.insert(boundaryverts.end(), boundary2verts.begin(), boundary2verts.end());

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

    //DEBUG
    // std::unique_ptr<MimmoObject> pippo = std::unique_ptr<MimmoObject>(new MimmoObject(1));
    // livector1D allinterf = mesh->getInterfaceFromVertexList(borderverts, true, true);
    // pippo->getPatch()->reserveVertices(borderverts.size());
    // pippo->getPatch()->reserveCells(allinterf.size());
    //
    // for(long id: borderverts){
    //     pippo->addVertex(mesh->getVertexCoords(id), id);
    // }
    // dvecarr3E normals;
    // normals.reserve(allinterf.size());
    //
    // for(long id: allinterf){
    //     normals.push_back(mesh->evalInterfaceNormal(id));
    //     bitpit::Interface & ii = orinterfaces.at(id);
    //     long * conn = ii.getConnect();
    //     std::size_t connsize = ii.getConnectSize();
    //
    //     bitpit::ElementType et = orcells.at(ii.getOwner()).getFaceType(ii.getOwnerFace());
    //
    //     pippo->addConnectedCell(std::vector<long>(&conn[0], &conn[connsize]), et);
    // }
    //
    // pippo->getPatch()->getVTK().addData("normals", bitpit::VTKFieldType::VECTOR, bitpit::VTKLocation::CELL, normals);
    // pippo->getPatch()->write("interfacenormal");

    return mesh;
}


// =================================================================================== //

int test1() {

    std::unique_ptr<MimmoObject> boundary;
    std::unique_ptr<MimmoObject> mesh = createTestVolumeMesh(boundary);


    bool check = false;
    long targetNode =  (11*11)*25 + (11)*5 + 5;

//TESTING THE SCALAR PROPAGATOR //////
    // and the scalar field of Dirichlet values on its nodes.
    MimmoPiercedVector<double> bc_surf_field;
    bc_surf_field.setGeometry(boundary.get());
    bc_surf_field.setDataLocation(MPVLocation::POINT);
    bc_surf_field.reserve(boundary->getNVertices());
    for(long val : boundary->getVertices().getIds()){
        bc_surf_field.insert(val, 0.0);
    }

    //extract vertex of boundary cells pidded 0 and assign to them a condition 10.0
    std::vector<long> listPP = boundary->getVertexFromCellList(boundary->extractPIDCells(0));
    for(long id : listPP){
        bc_surf_field.at(id) = 10.0;
    }
    // dvector1D data = bc_surf_field.getDataAsVector();
    // boundary->getPatch()->getVTK().addData("dirichlet", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, data);
    // boundary->getPatch()->write("boundary");
    // boundary->getPatch()->getVTK().removeData("dirichlet");
    //
    // Now create a PropagateScalarField and solve the laplacian.
    PropagateScalarField * prop = new PropagateScalarField();
    prop->setName("mimmo.PropagateScalarField_test2");
    prop->setGeometry(mesh.get());
    prop->setDirichletBoundarySurface(boundary.get());
    prop->setDirichletConditions(bc_surf_field);
    prop->setDumping(false);
    prop->setPlotInExecution(true);

    prop->exec();

    auto values = prop->getPropagatedField();
    std::cout<< values.at(targetNode) <<std::endl;
    check = check || (std::abs(values.at(targetNode)-4.99937) > 1.0E-3);

    delete prop;

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
