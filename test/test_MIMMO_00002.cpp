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

void test0002() {

	//Creation of MiMMO container.
	MimmoObject mimmo0;

	int		np,	nt;
	darray3E point;
	{
		//Import STL
		STLObj stl("geo_data/pipe.stl", true);
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

	//writing undeformed geometry stock in mimmo container
	string filename = "mimmo_pipe0";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	//********************************************************************************************
	// 	//CREATING GLOBAL DISPLACEMENTS LATTICE
	//Instantiation of a FFDobject of cylindrical shape
	FFDLattice* lattice = new FFDLattice();

	//Set cylindrical lattice
	darray3E origin = {0.0, 0.0,0.0};
	darray3E span;
	span[0]= 0.51;
	//span[1]= 6*std::atan(1.0);
	span[1]= 2*M_PI;
	span[2]= 8.51;

	ivector1D dim(3), deg(3);
	dim[0] = 2;
	dim[1] = 15;
	dim[2] = 40;

	deg[0] = 1;
	deg[1] = 3;
	deg[2] = 8;

// 	dim[0] = 2;
// 	dim[1] = 5;
// 	dim[2] = 2;
// 	deg[0] = 1;
// 	deg[1] = 1;
// 	deg[2] = 1;
	
	
	//set lattice
	lattice->setMesh(origin,span,BasicShape::ShapeType::CYLINDER,dim, deg);
	lattice->setRefSystem(2, darray3E{0,1,0});
	//lattice->setInfLimits(1.0*M_PI,1);
	//Set geometry
	lattice->setGeometry(&mimmo0);

	// 	//Set release Info
	lattice->setReleaseInfo(false);

	//Set Input with Init Displacements
	int ndof = lattice->getNDeg();
	dvecarr3E displ(ndof, darray3E{0,0,0});
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndof; i++){
		int l1,l2,l3;
		int index = lattice->accessGridFromDOF(i);
		lattice->accessPointIndex(index,l1,l2,l3);
 		if(l1>0){
				displ[i][0] = 0.5*( (double) (rand()) / RAND_MAX );
 		}
		// 			if(l2 == nnn[1]-1){
		// 				displ[i][0] = 1.0;
		// 			}

		//	}
	}

	lattice->setDisplacements(displ);
	lattice->setDisplGlobal(false);

	//********************************************************************************************
	// 	//CREATING LOCAL DISPLACEMENTS LATTICE
	//Instantiation of a FFDobject of cylindrical shape 
	FFDLattice* lattice2 = new FFDLattice();

	//Set cylindrical lattice
	span[0]= 0.51;
	//span[1]= 6*std::atan(1.0);
	span[1]= 2*M_PI;
	span[2]= 8.51;

	dim[0] = 10;
	dim[1] = 10;
	dim[2] = 40;

	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

	//set lattice
	lattice2->setMesh(origin,span,BasicShape::ShapeType::CYLINDER,dim, deg);
	lattice2->setRefSystem(2, darray3E{0,1,0});
	//Set geometry
	lattice2->setGeometry(&mimmo0);

	// 	//Set release Info
	lattice2->setReleaseInfo(true);

	//Set Input with Init Displacements
	ndof = lattice2->getNDeg();
	dvecarr3E displ2(ndof, darray3E{0,0,0});
	Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndof; i++){
		int l1,l2,l3;
		lattice2->accessPointIndex(i,l1,l2,l3);
		if(l1>0){
			displ2[i][0] = 0.0*( (double) (rand()) / RAND_MAX );
		}
		// 			if(l2 == nnn[1]-1){
		// 				displ[i][0] = 1.0;
		// 			}

		//	}
	}

	lattice2->setDisplacements(displ2);
	lattice2->setDisplGlobal(false);


	//********************************************************************************************
	//CREATING INPUT	
	// 	InputDoF* input = new InputDoF(ndeg, displ);
	// 	string file = "input.txt";
	// 	InputDoF* input0 = new InputDoF(file);
	//
	// 	cout << "input setup done" << endl;
	//********************************************************************************************
	//CREATE APPLIERS
	cout << "applier setup" << endl;
	Apply* applier = new Apply(&mimmo0);
	Apply* applier2 = new Apply(&mimmo0);

	lattice2->addChild(applier2);
	lattice->addChild(applier);

	cout << "applier setup done" << endl;

	//********************************************************************************************
	//CREATE FILTER MASK
	Mask* mask = new Mask();
	mask->setNDeg(lattice->getNDeg());
	darray3E thres;
	thres[0] = -10.0;
	thres[1] = 0.0;
	thres[2] = -10.0;
	mask->setThresholds(thres);
	mask->setForward(0,false);
	mask->setForward(1,false);
	mask->setForward(2,false);
	//set filter to lattice
	mask->addChild(lattice);

	//********************************************************************************************
	//CREATE BENDER-WRAPPER
	//********************************************************************************************
	Bend* bend = new Bend();
	dvecarr3E degree(3);
	degree[0][1] = 2;
	degree[1][0] = 2;
	bend->setDegree(degree);
	dvector3D coeffs(3, vector<vector<double> >(3) );
	coeffs[0][1].resize(degree[0][1]+1);
	coeffs[0][1][0] = 0.0;
	coeffs[0][1][1] = 0.0;
	coeffs[0][1][2] = 0.2;
	coeffs[1][0].resize(degree[1][0]+1);
	coeffs[1][0][0] = 0.0;
	coeffs[1][0][1] = -0.5;
	coeffs[1][0][2] = 0.0;
	bend->setCoeffs(coeffs);
	//set bend to lattice
	bend->addChild(mask);
	bend->setDisplacements(displ);

	//********************************************************************************************
	//CREATE BENDER-WRAPPER
	//********************************************************************************************
	Bend* bend2 = new Bend();
	degree.clear();
	degree.resize(3);
	degree[0][1] = 2;
	bend2->setDegree(degree);
	coeffs.clear();
	coeffs.resize(3, vector<vector<double> >(3));
	coeffs[0][1].resize(degree[0][1]+1);
	coeffs[0][1][0] = 0.0;
	coeffs[0][1][1] = 0.1;
	coeffs[0][1][2] = 0.05;
	bend2->setCoeffs(coeffs);
	//set bend to lattice
	bend2->addChild(lattice2);
	bend2->setDisplacements(displ2);

	//create output
	cout << "output setup" << endl;
	OutputDoF* output = new OutputDoF();
	lattice->addChild(output);
	cout << "output setup done" << endl;

	//********************************************************************************************

	//Creating ELEMENT chain
	Chain ch0;
	//ch0.addObject(lattice2);
	ch0.addObject(output);
	ch0.addObject(lattice);
//	ch0.addObject(mask);
//	ch0.addObject(bend);
	//ch0.addObject(applier2);
	ch0.addObject(applier);

	//********************************************************************************************
	//Executing CHAIN

	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();

	ch0.exec();

	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//********************************************************************************************
	//PLOT RESULTS

	lattice->plotGrid("./", "lattice_pipe", 0, false, false);
	lattice->plotGrid("./", "lattice_pipe", 1, false, true);
// 	lattice2->plotGrid("./", "lattice2_pipe", 0, false, false);
// 	lattice2->plotGrid("./", "lattice2_pipe", 1, false, true);

	filename = "mimmo_pipe1";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	//********************************************************************************************
	//clean up & exit;
	delete lattice, applier, output, applier2, lattice2, mask, bend;

	lattice = NULL;
	applier = NULL;
	lattice2 = NULL;
	applier2 = NULL;
	mask = NULL;
	bend = NULL;

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
		test0002() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return(0);
}
