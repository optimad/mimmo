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
	MimmoObject mimmoPlane, mimmoDisk;
	//Input triangulation
	int		np,	nt;
	darray3E point;

	{
		//Import STL
		STLObj stl("geo_data/plane.stl", true);
		dvector2D V,N;
		ivector2D T;
		stl.load(np, nt, V, N, T);
		for (long ip=0; ip<np; ip++){
			point = conArray<double,3>(V[ip]);
			mimmoPlane.setVertex(point);
		}
		mimmoPlane.setConnectivity(&T);
		mimmoPlane.cleanGeometry();
	}
	string filename = "mimmo_0006p.0000";
	mimmoPlane.m_geometry->setName(filename);
	mimmoPlane.m_geometry->write();


	{
		//Import STL
		STLObj stl("geo_data/disk.stl", true);
		dvector2D V,N;
		ivector2D T;
		stl.load(np, nt, V, N, T);
		for (long ip=0; ip<np; ip++){
			point = conArray<double,3>(V[ip]);
			mimmoDisk.setVertex(point);
		}
		mimmoDisk.setConnectivity(&T);
		mimmoDisk.cleanGeometry();
	}
	//Write undeformed geometry
	filename = "mimmo_0006d.0000";
	mimmoDisk.m_geometry->setName(filename);
	mimmoDisk.m_geometry->write();

	//Instantiation of a FFDobject (and Input object).
	FFDLattice* lattice = new FFDLattice();
	//Set lattice
	lattice->setGeometry(&mimmoDisk);
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
	coeffs[0] = -0.25;
//	coeffs[1] = -1;
//	coeffs[2] = 15;
	coeffs[1] = 0;
	coeffs[2] = 0;
//	bend->setDegree(0,1,degree);
//	bend->setCoeffs(0,1,coeffs);
	bend->setDegree(2,1,degree);
	bend->setCoeffs(2,1,coeffs);

	//create applier
	Apply* applier = new Apply();
	applier->setGeometry(&mimmoDisk);

	//createRBF
	MRBF* mrbf = new MRBF();
	mrbf->setGeometry(&mimmoPlane);
	mrbf->addNode(&mimmoDisk);
	mrbf->setTol(0.00001);

	//create applier
	Apply* applier2 = new Apply();
	applier2->setGeometry(&mimmoPlane);


	//Set PINS
	cout << "set pins" << endl;

	addPin(inputshapet, mesh, &GenericInput::getResult<int>, &Lattice::setShape);
	addPin(inputorig, mesh, &GenericInput::getResult<darray3E>, &Lattice::setOrigin);
	addPin(inputspan, mesh, &GenericInput::getResult<darray3E>, &Lattice::setSpan);
	addPin(inputdim, mesh, &GenericInput::getResult<iarray3E>, &Lattice::setDimension);

	addPin(input, bend, &GenericInput::getResult<dvecarr3E>, &Bend::setInput<dvecarr3E>);
	addPin(mesh, bend, &Lattice::getGlobalCoords, &Bend::setCoords);

	addPin(bend, lattice, &Bend::getResult<dvecarr3E>, &FFDLattice::setDisplacements);

	addPin(inputshapet, lattice, &GenericInput::getResult<int>, &FFDLattice::setShape);
	addPin(inputorig, lattice, &GenericInput::getResult<darray3E>, &FFDLattice::setOrigin);
	addPin(inputspan, lattice, &GenericInput::getResult<darray3E>, &FFDLattice::setSpan);
	addPin(inputdim, lattice, &GenericInput::getResult<iarray3E>, &FFDLattice::setDimension);
	addPin(inputdeg, lattice, &GenericInput::getResult<iarray3E>, &FFDLattice::setDegrees);

	addPin(lattice, applier, &FFDLattice::getResult<dvecarr3E>, &Apply::setInput<dvecarr3E>);

	addPin(lattice, mrbf, &FFDLattice::getResult<dvecarr3E>, &MRBF::setDisplacements);
	addPin(mrbf, applier2, &MRBF::getResult<dvecarr3E>, &Apply::setInput<dvecarr3E>);

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
	ch0.exec();
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//Plot results
	lattice->plotGrid("./", "lattice_0006", 0, false, false);
	lattice->plotGrid("./", "lattice_0006", 1, false, true);
	filename = "mimmo_0006d.0001";
	mimmoDisk.m_geometry->setName(filename);
	mimmoDisk.m_geometry->write();

	filename = "mimmo_0006p.0001";
	mimmoPlane.m_geometry->setName(filename);
	mimmoPlane.m_geometry->write();

	//Delete and nullify pointer
	delete lattice;
	delete applier;
	delete mrbf;
	delete applier2;
	delete input;
	delete inputorig;
	delete inputspan;
	delete inputshapet;
	delete inputdim;
	delete inputdeg;

	lattice 	= NULL;
	applier 	= NULL;
	mrbf 		= NULL;
	applier2 	= NULL;
	inputorig 	= NULL;
	inputspan 	= NULL;
	inputshapet	= NULL;
	inputdim 	= NULL;
	inputdeg 	= NULL;
	input 		= NULL;

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
