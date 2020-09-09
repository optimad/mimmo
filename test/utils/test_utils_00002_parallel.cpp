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

/*
 * Creating surface triangular mesh and return it in a MimmoObject.
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(mimmo::MimmoObject * mesh){
    if(mesh->getRank() == 0){
        double dx = 0.25, dy = 0.25;
        int nV, nC;
        //create vertexlist
        dvecarr3E vertex(35,{{0.0,0.0,0.0}});
        livector2D conn(48, livector1D(3));

        for(int i=0; i<7; ++i){
            for(int j=0; j<5; j++){
                nV = 5*i + j;
                vertex[nV][0] = i*dx;
                vertex[nV][1] = j*dy;
            }
        }

        for(int j=0; j<4; ++j){
            for(int i=0; i<3; ++i){
                nC = 8*i + 2*j;

                conn[nC][0] = 5*i + j;
                conn[nC][1] = 5*(i+1) + j;
                conn[nC][2] = 5*i + j+1;

                conn[nC+1][0] = 5*(i+1) + j;
                conn[nC+1][1] = 5*(i+1) + j+1;
                conn[nC+1][2] = 5*i + j+1;
            }
        }

        for(int j=0; j<4; ++j){
            for(int i=3; i<6; ++i){
                nC = 8*i + 2*j;

                conn[nC][0] = 5*i + j;
                conn[nC][1] = 5*(i+1) + j;
                conn[nC][2] = 5*(i+1) + j+1;

                conn[nC+1][0] = 5*i + j;
                conn[nC+1][1] = 5*(i+1) + j+1;
                conn[nC+1][2] = 5*i + j+1;
            }
        }

        mesh->getVertices().reserve(35);
        mesh->getCells().reserve(48);

        //fill the mimmoObject;
        long cV=0;
        for(auto & val: vertex){
            mesh->addVertex(val, cV);
            cV++;
        }

        long cC=0;
        bitpit::ElementType eltype = bitpit::ElementType::TRIANGLE;
        for(auto & val: conn){
            mesh->addConnectedCell(val, eltype, 0, cC, 0);
            cC++;
        }
    }

    mesh->updateAdjacencies();
    mesh->update();

    bool check = (mesh->getNGlobalCells() == 48) && (mesh->getNGlobalVertices() == 35);
    return check;
}
// =================================================================================== //
/*
 * Test: testing CreateSeedsOnSurface utility
 */
int test2() {

    mimmo::Partition * part = new mimmo::Partition();
    bitpit::Logger & log = part->getLog();
    log.setPriority(bitpit::log::Priority::NORMAL);

    mimmo::MimmoSharedPointer<mimmo::MimmoObject> m1(new mimmo::MimmoObject());
    if(!createMimmoMesh(m1.get())) {
        return 1;
    }

    part->setName("support_t2p_utils");
    part->setGeometry(m1);
    part->setPlotInExecution(true);
    part->exec();


    darray3E normal = {{0,0,1}};
    mimmo::CreateSeedsOnSurface * cseed = new mimmo::CreateSeedsOnSurface();
    cseed->setGeometry(part->getGeometry());
    cseed->setSeed({{0.0,0.0,0.0}});
    cseed->setNPoints(16);
    cseed->setEngine(0); //random engine
    cseed->setRandomFixed(true);
    cseed->setRandomSignature(123456789);
    cseed->setPlotInExecution(true);

    cseed->exec();

    //std::cout<<cseed->getRank()<<"----"<<cseed->getPoints()<<std::endl;

    bool check = std::abs(dotProduct(cseed->getPoints()[1], normal)) < 1.e-18;


    delete cseed;
    delete part;
    log<<"test passed :" <<check<<std::endl;
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
            val = test2() ;
        }
        catch(std::exception & e){
            std::cout<<"test_utils_00002_parallel exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
