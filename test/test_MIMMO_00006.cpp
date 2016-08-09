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

void test0006() {

	//Creation of MiMMO container.
	MimmoGeometry * mimmoP = new MimmoGeometry();
	MimmoGeometry * mimmoD = new MimmoGeometry();
	
	mimmoP->setRead(true);
	mimmoP->setReadDir("geo_data");
	mimmoP->setReadFileType(FileType::STL);
	mimmoP->setReadFilename("plane");
	
	mimmoP->setWrite(true);
	mimmoP->setWriteDir(".");
	mimmoP->setWriteFileType(FileType::STL);
	mimmoP->setWriteFilename("mimmo_0006p.0000");
	
	mimmoP->execute();
	MimmoObject * objectPlane = mimmoP->getGeometry();

	mimmoD->setRead(true);
	mimmoD->setReadDir("geo_data");
	mimmoD->setReadFileType(FileType::STL);
	mimmoD->setReadFilename("disk");
	
	mimmoD->setWrite(true);
	mimmoD->setWriteDir(".");
	mimmoD->setWriteFileType(FileType::STL);
	mimmoD->setWriteFilename("mimmo_0006d.0000");
	
	mimmoD->execute();
	MimmoObject * objectDisk = mimmoD->getGeometry();
	

	//Instantiation of a FFDobject (and Input object).
	FFDLattice* lattice = new FFDLattice();
	//Set lattice
	lattice->setGeometry(objectDisk);
	lattice->setDisplGlobal(true);


	//Set Inputs with Shape and Mesh Info
	darray3E origin = {0.05, 0.05, 0.0};
	darray3E span;
	span[0]= 0.1;
	span[1]= 2*M_PI;
	span[2]= 0.1;

	iarray3E dim, deg;
	dim[0] = 2;
	dim[1] = 40;
	dim[2] = 3;
	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

	int t = 1;
	GenericInput* inputshapet = new GenericInput();
	inputshapet->setInput(t);
	inputshapet->setName("MiMMO.InputShape");

	GenericInput* inputorig = new GenericInput();
	inputorig->setInput(origin);
	inputorig->setName("MiMMO.InputOrigin");

	GenericInput* inputspan = new GenericInput();
	inputspan->setInput(span);
	inputspan->setName("MiMMO.InputSpan");

	GenericInput* inputdim = new GenericInput();
	inputdim->setInput(dim);
	inputdim->setName("MiMMO.InputDim");

	GenericInput* inputdeg = new GenericInput();
	inputdeg->setInput(deg);
	inputdeg->setName("MiMMO.InputDeg");

	//Set Input with Init Displacements
	int ndeg = (dim[0])*(dim[1])*(dim[2]);
	dvecarr3E displ(ndeg);
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
			displ[i][0] = 0.0*( (double) (rand()) / RAND_MAX - 0.5);
	}
	GenericInput* input = new GenericInput();
	input->setInput(displ);
	input->setName("MiMMO.InputDispl");

	//create aux lattice for mesh and nodes coordinates
	Lattice* mesh = new Lattice();

	//create Bend
	Bend* bend = new Bend();
	uint32_t	degree = 2;
	dvector1D	coeffs(degree+1);
	coeffs[0] = -0.0625;
	coeffs[1] = 0;
	coeffs[2] = 0;
	bend->setDegree(2,1,degree);
	bend->setCoeffs(2,1,coeffs);

	//create applier
	Apply* applier = new Apply();
	applier->setGeometry(objectDisk);

	//createRBF
	MRBF* mrbf = new MRBF();
	mrbf->setMode(MRBFSol::NONE);
	mrbf->setGeometry(objectPlane);
	mrbf->setSupportRadius(-1.0);
	mrbf->addNode(objectDisk);

	//create applier
	Apply* applier2 = new Apply();
	applier2->setGeometry(objectPlane);

	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(inputshapet, mesh, PortType::M_VALUEI, PortType::M_SHAPEI) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(inputorig, mesh, PortType::M_POINT, PortType::M_POINT) << endl;
	cout << "add pin info 3 : " << boolalpha << addPin(inputspan, mesh, PortType::M_SPAN, PortType::M_SPAN) << endl;
	cout << "add pin info 4 : " << boolalpha << addPin(inputdim, mesh, PortType::M_DIMENSION, PortType::M_DIMENSION) << endl;

	cout << "add pin info 5 : " << boolalpha << addPin(input, bend, PortType::M_DISPLS, PortType::M_DISPLS) << endl;
	cout << "add pin info 6 : " << boolalpha << addPin(mesh, bend, PortType::M_GLOBAL, PortType::M_COORDS) << endl;

	cout << "add pin info 7 : " << boolalpha << addPin(bend, lattice, PortType::M_DISPLS, PortType::M_DISPLS) << endl;

	cout << "add pin info 8 : " << boolalpha << addPin(inputshapet, lattice, PortType::M_VALUEI, PortType::M_SHAPEI) << endl;
	cout << "add pin info 9 : " << boolalpha << addPin(inputorig, lattice, PortType::M_POINT, PortType::M_POINT) << endl;
	cout << "add pin info 10 : " << boolalpha << addPin(inputspan, lattice, PortType::M_SPAN, PortType::M_SPAN) << endl;
	cout << "add pin info 11 : " << boolalpha << addPin(inputdim, lattice, PortType::M_DIMENSION, PortType::M_DIMENSION) << endl;
	cout << "add pin info 12 : " << boolalpha << addPin(inputdeg, lattice, PortType::M_DEG, PortType::M_DEG) << endl;

	cout << "add pin info 13 : " << boolalpha << addPin(lattice, applier, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;

	cout << "add pin info 14 : " << boolalpha << addPin(lattice, mrbf, PortType::M_GDISPLS, PortType::M_DISPLS) << endl;
	cout << "add pin info 15 : " << boolalpha << addPin(mrbf, applier2, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;

	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(inputorig);
	ch0.addObject(inputshapet);
	ch0.addObject(inputspan);
	ch0.addObject(inputdim);
	ch0.addObject(inputdeg);
	ch0.addObject(input);
	ch0.addObject(mesh);
	ch0.addObject(bend);
	ch0.addObject(lattice);
	ch0.addObject(applier);
	ch0.addObject(mrbf);
	ch0.addObject(applier2);

	//Execution of chain
	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	ch0.exec(true);
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//Plot results
	lattice->plotGrid("./", "lattice_0006", 0, false, false);
	lattice->plotGrid("./", "lattice_0006", 1, false, true);
	
	
	mimmoP->setWriteFilename("mimmo_0006p.0001");
	mimmoP->write();
	mimmoD->setWriteFilename("mimmo_0006d.0001");
	mimmoD->write();

	//Delete and nullify pointer
	delete lattice;
	delete mesh;
	delete bend;
	delete applier;
	delete mrbf;
	delete applier2;
	delete input;
	delete inputorig;
	delete inputspan;
	delete inputshapet;
	delete inputdim;
	delete inputdeg;
	delete mimmoP;
	delete mimmoD;
	
	lattice 	= NULL;
	mesh 		= NULL;
	bend	 	= NULL;
	applier 	= NULL;
	mrbf 		= NULL;
	applier2 	= NULL;
	inputorig 	= NULL;
	inputspan 	= NULL;
	inputshapet	= NULL;
	inputdim 	= NULL;
	inputdeg 	= NULL;
	input 		= NULL;
	mimmoP 		= NULL;
	mimmoD 		= NULL;
	objectPlane = NULL;
	objectDisk 	= NULL;
	
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

		test0006() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
