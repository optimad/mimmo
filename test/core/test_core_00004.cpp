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

#include "mimmo_core.hpp"

/*
 * Test 00004
 * Testing BasicMeshes and Lattice construction
 */

// =================================================================================== //

int test4() {

	mimmo::Lattice * mesh = new mimmo::Lattice();
	mesh->setShape(mimmo::ShapeType::CUBE);
	mesh->setOrigin({{0.5,0.5,0.5}});
	mesh->setSpan(1.0,1.0,1.0);
	iarray3E dim = {{3,2,3}};
	mesh->setDimension(dim);
	mesh->execute();

	bool check = true;

	check = ((int)mesh->getGlobalCellCentroids().size() == 4);

	dvector1D temp(2, 0.25);
	temp[1] +=0.5;

	darray3E vec;

	if(check){
		for(int i=0; i<2; ++i){
			for(int j=0; j<1; ++j){
				for(int k=0; k<2; ++k){
					vec[0] = temp[i];
					vec[1] = 0.5;
					vec[2] = temp[k];
					check =  check && (norm2(vec - mesh->getGlobalCCell(i,j,k)) <= 1.e-15);
				}
			}
		}
	}

	if(!check){
		std::cout<<"ERROR.Not able to build Lattice Structured Mesh"<<std::endl;
		delete mesh;
        return 1;
	}else{
		mesh->plotGrid(".","lattice_t4", 0, 0 );
		std::cout<<"Successfully built Lattice Structured Mesh and written to file lattice_t4.0000.vtu"<<std::endl;
	}

	delete mesh;

	mimmo::Lattice * mesh2 = new mimmo::Lattice();
    mesh2->setShape(mimmo::ShapeType::CYLINDER);
    mesh2->setOrigin({{-1920.0, 0.0, 0.0}});
    mesh2->setSpan({{560.0, 6.28, 310.0}});
    mesh2->setInfLimits({{220.0, 0.0, 0.0}});
    iarray3E dim2 = {{2, 4, 2}};
    mesh2->setDimension(dim2);
    mesh2->setPlotInExecution(true);
    mesh2->exec();

    delete mesh2;

    return 0;
}

// =================================================================================== //

int test5() {

	mimmo::Lattice * mesh = new mimmo::Lattice();
    mesh->setShape(mimmo::ShapeType::WEDGE);
    mesh->setOrigin({{0.0,0.0,1.5}});
    mesh->setSpan(1.2,0.5,3.0);
    iarray3E dim = {{2,2,2}};
    mesh->setDimension(dim);
    mesh->execute();

    bool check = true;

    check = ((int)mesh->getGlobalCellCentroids().size() == 1);

    if(!check){
        std::cout<<"ERROR.Not able to build Lattice Structured Mesh"<<std::endl;
        delete mesh;
        return 1;
    }else{
        mesh->plotGrid(".","lattice_t4", 1, 0 );
        std::cout<<"Successfully built Lattice Structured Mesh and written to file lattice_t4.0001.vtu"<<std::endl;
    }

    delete mesh;
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
        val = test4() ;
        std::max(val, test5());
    }
    catch(std::exception & e){
        std::cout<<"test_core_00004 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
