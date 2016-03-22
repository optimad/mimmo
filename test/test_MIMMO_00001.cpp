/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "MiMMO.hpp"

using namespace std;
using namespace bitpit;


#include <chrono>

using namespace std::chrono;

// =================================================================================== //

void test0001() {



	//Creation of MiMMO container.
	MimmoObject mimmo0;
	//Input triangulation
	int		np,	nt;
	darray3E point;
	{
		//Import STL
//		STLObj stl("geo_data/placca.stl", true);
//		STLObj stl("geo_data/placca0.stl", true);
//		STLObj stl("geo_data/sphere.stl", true);
		STLObj stl("geo_data/sphere2.stl", true);
//		STLObj stl("geo_data/cad.stl", true);

		dvector2D V,N;
		ivector2D T;
		stl.load(np, nt, V, N, T);
		
		for (long ip=0; ip<np; ip++){
			point = conArray<double,3>(V[ip]);
			mimmo0.setVertex(point);
		}
		mimmo0.setConnectivity(&T);
		mimmo0.cleanGeometry();
	}
	//Write undeformed geometry
	string filename = "mimmo0";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();
	
	//Instantiation of a FFDobject (and Input object).
	FFDLattice* lattice = new FFDLattice();
	//Set lattice

	//sphere2
	darray3E origin = {-0.0, -0.0,-0.0};
	darray3E span;
	span[0]= 1.2;
	span[1]= 1.2;
	span[2]= 1.2;


	//Set Lattice dimensions and degree
	ivector1D dim(3), deg(3);
	dim[0] = 20;
	dim[1] = 20;
	dim[2] = 20;
	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;


	lattice->setMesh(origin,span,BasicShape::ShapeType::CUBE,dim, deg);

	//Set geometry
	lattice->setGeometry(&mimmo0);

	//Set Input with Init Displacements
	np = lattice->getNNodes();
	dvecarr3E displ(np);
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<np; i++){
		for (int j=0; j<3; j++){
			displ[i][j] = 0.15*( (double) (rand()) / RAND_MAX );
		}
	}

	GenericInput* input = new GenericInput();
	input->setInput(displ);

	//create applier
	Apply* applier = new Apply();
	applier->setGeometry(&mimmo0);


	//Create PINS
	addPin(input, lattice, &GenericInput::getResult<dvecarr3E>, &FFDLattice::setDisplacements);
	addPin(lattice, applier, &FFDLattice::getResult<dvecarr3E>, &Apply::setInput<dvecarr3E>);

	//Create chain
	Chain ch0;
	cout << "add input" << endl;
	ch0.addObject(input);
	cout << "add lattice" << endl;
	ch0.addObject(lattice);
	cout << "add applier" << endl;
	ch0.addObject(applier);

	//Execution of chain
	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	ch0.exec();
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//Plot results
	lattice->plotGrid("./", "lattice", 0, false, false);
	lattice->plotGrid("./", "lattice", 1, false, true);
	filename = "mimmo1";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	//Delete and nullify pointer
	delete lattice, applier, input;
	lattice = NULL;
	applier = NULL;
	input 	= NULL;

	//Print execution time
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

    return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

        test0001() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
