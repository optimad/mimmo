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
#include "customOperators.hpp"
using namespace std;
using namespace bitpit;


#include <chrono>

using namespace std::chrono;
using namespace mimmo;
using namespace mimmo::pin;

// =================================================================================== //

void test0001() {


	//Creation of MiMMO container.
	MimmoGeometry * mimmo0 = new MimmoGeometry();

	mimmo0->setRead(true);
	mimmo0->setReadDir("geo_data");
	mimmo0->setReadFileType(FileType::STL);
	mimmo0->setReadFilename("sphere2");
	mimmo0->setBuildBvTree(true);

	mimmo0->setWrite(true);
	mimmo0->setWriteDir("./");
	mimmo0->setWriteFileType(FileType::STL);
	mimmo0->setWriteFilename("mimmo_0001.0000");

	mimmo0->execute();
	MimmoObject * object = new MimmoObject();
	object->setHARDCopy(mimmo0->getGeometry());

	delete mimmo0; mimmo0=NULL;

	//Instantiation of a FFDobject (and Input object).
	FFDLattice* lattice = new FFDLattice();
	//Set lattice
	darray3E origin = {-0.0, -0.0,-0.0};
	darray3E span;
	span[0]= 1.2;
	span[1]= 1.2;
	span[2]= 1.2;

	//Set Lattice dimensions and degree
	iarray3E dim, deg;
	dim[0] = 20;
	dim[1] = 20;
	dim[2] = 20;
	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

	lattice->setLattice(origin,span,ShapeType::CUBE,dim, deg);
	lattice->setGeometry(object);

	GenericInput* input = new GenericInput();
	input->setReadFromFile(true);
	input->setFilename("input/input_MIMMO_00001.txt");

	//create applier
	Apply* applier = new Apply();
	applier->setGeometry(object);

	cout << "create pin " << endl;
	//Create PINS
	//Add pin with port ID
	cout << "add pin info : " << boolalpha << addPin(input, lattice, 10, 10) << endl;
	//Add pin with port TAG
	cout << "add pin info : " << boolalpha << addPin(lattice, applier, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;

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
	lattice->plotGrid("./", "lattice_0001", 0, true, false);
	lattice->plotGrid("./", "lattice_0001", 1, false, true);

	mimmo0 = new MimmoGeometry();

	mimmo0->setWrite(true);
	mimmo0->setWriteDir("./");
	mimmo0->setWriteFileType(FileType::STL);
	mimmo0->setWriteFilename("mimmo_0001.0001");
	mimmo0->setGeometry(object);
	mimmo0->execute();

	//Delete and nullify pointer
	delete lattice, applier, input, mimmo0, object;
	lattice = NULL;
	applier = NULL;
	input 	= NULL;
	object  = NULL;
	mimmo0  = NULL;

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
