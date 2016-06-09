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
using namespace mimmo;


#include <chrono>

using namespace std::chrono;

// =================================================================================== //

void test0003() {

	//Creation of MiMMO container.
	MimmoGeometry * mimmo0 = new MimmoGeometry();
	
	mimmo0->setRead(true);
	mimmo0->setReadDir("geo_data");
	mimmo0->setReadFileType(mimmo::FileType::STL);
	mimmo0->setReadFilename("ball");
	
	mimmo0->setWrite(true);
	mimmo0->setWriteDir(".");
	mimmo0->setWriteFileType(mimmo::FileType::STL);
	mimmo0->setWriteFilename("mimmo_0003.0000");
	
	mimmo0->execute();
	MimmoObject * object = mimmo0->getGeometry();

	//********************************************************************************************
	// 	//CREATING LATTICE
	//Instantiation of a FFDobject of spherical shape 
	FFDLattice* lattice = new FFDLattice();

	//Set spherical lattice
	darray3E origin = {0.0, 0.0,0.0};
	darray3E span;
	span[0]= 3.01;
	span[1]= 2*M_PI;
	span[2]= M_PI;

	iarray3E dim, deg;
	dim[0] = 30;
	dim[1] = 30;
	dim[2] = 30;

	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

	//set lattice
	lattice->setLattice(origin,span,ShapeType::SPHERE,dim, deg);
	lattice->setRefSystem(2, darray3E{0,1,0});
	lattice->setCoordType(CoordType::CLAMPED, 2);
	lattice->setGeometry(object);

	//Set Input with Init Displacements
	int ndeg = lattice->getNNodes();
	dvecarr3E displ(ndeg, darray3E{0,0,0});
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
		int l1,l2,l3;
		int index = lattice->accessGridFromDOF(i);
		lattice->accessPointIndex(index,l1,l2,l3);
		if(l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] < M_PI){
			displ[i][0] = 1.0*( (double) (rand()) / RAND_MAX );
		}
		if( (l1 > 0 && lattice->getLocalPoint(l1,l2,l3)[1] >= M_PI)
				|| lattice->getLocalPoint(l1,l2,l3)[1] == 0){
			displ[i][0] = 1.5;
		}

	}

	//********************************************************************************************
	//	CREATING INPUT
	cout << "input setup" << endl;
	GenericInput* input = new GenericInput();
	input->setInput(displ);
	cout << "input setup done" << endl;

	//********************************************************************************************
	//CREATE APPLIER
	cout << "applier setup" << endl;
	Apply* applier = new Apply();
	applier->setGeometry(object);
	cout << "applier setup done" << endl;

	//********************************************************************************************
	//CREATE FILTER MASK

	//********************************************************************************************
	//CREATE BENDER-WRAPPER

	//********************************************************************************************
	//SETUP PINS
	addPin(input, lattice, M_DISPLS, M_DISPLS);
	addPin(lattice, applier, M_GDISPLS, M_GDISPLS);

	//********************************************************************************************
	//Creating ELEMENT chain
	Chain ch0;
	ivector1D chain_pos;
	chain_pos.push_back(ch0.addObject(input));
	chain_pos.push_back(ch0.addObject(applier));
	chain_pos.push_back(ch0.addObject(lattice));

	//********************************************************************************************
	//Executing CHAIN

	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();

	ch0.exec();

	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//********************************************************************************************
	//PLOT RESULTS

	lattice->plotGrid("./", "lattice_0003", 0, false, false);
	lattice->plotGrid("./", "lattice_0003", 1, false, true);

	mimmo0->setRead(false);
	mimmo0->setWriteFilename("mimmo_0003.0001");
	mimmo0->execute();

	//********************************************************************************************
	//clean up & exit;
	delete lattice, applier, input, mimmo0;

	lattice = NULL;
	applier = NULL;
	input 	= NULL;
	object = NULL;
	mimmo0 = NULL;
	

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO Deformation execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
}

// =================================================================================== //

int	main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/
		test0003() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return(0);
}
