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
#include <functional>
#include "customOperators.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test0004() {

	//Creation of MiMMO container.
	MimmoGeometry * mimmo0 = new MimmoGeometry();
	
	mimmo0->setRead(true);
	mimmo0->setReadDir("geo_data");
	mimmo0->setReadFileType(mimmo::FileType::STL);
	mimmo0->setReadFilename("sphere2");
	
	mimmo0->setWrite(true);
	mimmo0->setWriteDir(".");
	mimmo0->setWriteFileType(mimmo::FileType::STL);
	mimmo0->setWriteFilename("mimmo_0004.0000");
	
	mimmo0->execute();
	MimmoObject * object = mimmo0->getGeometry();

	//Instantiation of a FFDobject (and Input object).
	FFDLattice* lattice = new FFDLattice();

	//Set lattice
	lattice->setGeometry(object);

	//create aux lattice for mesh and nodes coordinates
	Lattice* mesh = new Lattice();

	//Set Inputs with Shape and Mesh Info
	darray3E origin = {0.0, 0.0, 0.0};
	darray3E span;
	span[0]= 1.2;
	span[1]= 1.2;
	span[2]= 1.2;

	iarray3E dim, deg;
	dim[0] = 20;
	dim[1] = 20;
	dim[2] = 20;
	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;


	//TODO Ports for generic input
	int t = 0;
	GenericInput* inputshapet = new GenericInput();
	inputshapet->setInput(t);
	mesh->setShape(t);

	GenericInput* inputorig = new GenericInput();
	inputorig->setInput(origin);
	mesh->setOrigin(origin);

	GenericInput* inputspan = new GenericInput();
	inputspan->setInput(span);
	mesh->setSpan(span);

	GenericInput* inputdim = new GenericInput();
	inputdim->setInput(dim);
	mesh->setDimension(dim);

	GenericInput* inputdeg = new GenericInput();
	inputdeg->setInput(deg);
	lattice->setDegrees(deg);

	GenericInput* inputname = new GenericInput();
	string name = "test_MIMMO_0004.out";
	inputname->setInput(name);

	GenericOutput* output = new GenericOutput();
	output->setFilename(name);

	//Set Input with Init Displacements
	int ndeg = (dim[0])*(dim[1])*(dim[2]);
	dvecarr3E displ(ndeg);
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
		for (int j=0; j<3; j++){
			displ[i][j] = 0.25*( (double) (rand()) / RAND_MAX - 0.5);
		}
	}
	GenericInput* input = new GenericInput();
	input->setInput(displ);

	//create Mask
	Mask* mask = new Mask();
	mask->setThresholds({{-10.0,0}}, 0);
	mask->setThresholds({{-10.0,10.0}}, 1);
	mask->setThresholds({{-10.0,10.0}}, 2);
	mask->setInside(0, true);
	mask->setInside(1, true);
	mask->setInside(2, true);

	//create Bend
	Bend* bend = new Bend();
	uint32_t	degree = 2;
	dvector1D	coeffs(degree+1);
	coeffs[0] = 0;
	coeffs[1] = 0;
	coeffs[2] = 1;
	bend->setDegree(2,0,degree);
	bend->setCoeffs(2,0,coeffs);

	//create applier
	Apply* applier = new Apply();
	applier->setGeometry(object);

	//Set PINS
	cout << "set pins" << endl;

	//TODO Ports for generic input
//	addPin(inputshapet, mesh, &GenericInput::getResult<int>, &Lattice::setShape);
//	addPin(inputorig, mesh, &GenericInput::getResult<darray3E>, &Lattice::setOrigin);
//	addPin(inputspan, mesh, &GenericInput::getResult<darray3E>, &Lattice::setSpan);
//	addPin(inputdim, mesh, &GenericInput::getResult<iarray3E>, &Lattice::setDimension);

	addPin(mesh, mask, GLOBAL, COORDS);
	addPin(input, mask, DISPLS, DISPLS);

	addPin(mask, bend, COORDS, COORDS);
	addPin(mask, bend, DISPLS, DISPLS);

	//TODO Ports for generic input
//	addPin(inputshapet, lattice, &GenericInput::getResult<int>, &FFDLattice::setShape);
//	addPin(inputorig, lattice, &GenericInput::getResult<darray3E>, &FFDLattice::setOrigin);
//	addPin(inputspan, lattice, &GenericInput::getResult<darray3E>, &FFDLattice::setSpan);
//	addPin(inputdim, lattice, &GenericInput::getResult<iarray3E>, &FFDLattice::setDimension);
//	addPin(inputdeg, lattice, &GenericInput::getResult<iarray3E>, &FFDLattice::setDegrees);

//	addPin(inputname, output, &GenericInput::getResult<string>, &GenericOutput::setFilename);

	addPin(bend, output, DISPLS, DISPLS);

	addPin(bend, lattice, DISPLS, DISPLS);
	addPin(lattice, applier, GDISPLS, GDISPLS);

	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs " << endl;
	ch0.addObject(inputorig);
	ch0.addObject(inputshapet);
	ch0.addObject(inputspan);
	ch0.addObject(inputdim);
	ch0.addObject(inputdeg);
	ch0.addObject(inputname);
	ch0.addObject(input);
	cout << "add mesh" << endl;
	ch0.addObject(mesh);
	cout << "add mask" << endl;
	ch0.addObject(mask);
	cout << "add bend" << endl;
	ch0.addObject(bend);
	cout << "add lattice" << endl;
	ch0.addObject(lattice);
	cout << "add output" << endl;
	ch0.addObject(output);
	cout << "add applier" << endl;
	ch0.addObject(applier);

	//Execution of chain
	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	ch0.exec();
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//Plot results
	lattice->plotGrid("./", "lattice_0004", 0, false, false);
	lattice->plotGrid("./", "lattice_0004", 1, false, true);
	
	mimmo0->setRead(false);
	mimmo0->setWriteFilename("mimmo_0004.0001");
	mimmo0->execute();

	//Delete and nullify pointer
	delete lattice;
	delete applier;
	delete input;
	delete inputorig;
	delete inputspan;
	delete inputshapet;
	delete inputdim;
	delete inputdeg;
	delete mimmo0;

	lattice 	= NULL;
	applier 	= NULL;
	inputorig 	= NULL;
	inputspan 	= NULL;
	inputshapet	= NULL;
	inputdim 	= NULL;
	inputdeg 	= NULL;
	input 		= NULL;
	mimmo0 = NULL;
	object = NULL;

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

		test0004() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
