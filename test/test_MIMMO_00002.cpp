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
	//CREATING LATTICE
	//Instantiation of a FFDobject of cylindrical shape 
	FFDLattice* lattice = new FFDLattice();

	//Set cylindrical lattice
	darray3E origin = {0.0, 0.0,0.0};
	darray3E span;
	span[0]= 0.6;
	span[1]= 8*std::atan(1.0);
	span[2]= 8.6;

	ivector1D dim(3), deg(3);
	dim[0] = 2;
	dim[1] = 20;
	dim[2] = 12;

	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

// 	//Set cylindrical lattice
// 	darray3E origin = {-0.51, -4.26,-0.51};
// 	darray3E span;
// 	span[0]= 1.02;
// 	span[1]= 8.52;
// 	span[2]= 1.02;
// 	
// 	ivector1D dim(3), deg(3);
// 	dim[0] = 2;
// 	dim[1] = 4;
// 	dim[2] = 2;
// 	
// 	deg[0] = 1;
// 	deg[1] = 1;
// 	deg[2] = 1;
	
	//set lattice
	lattice->setMesh(origin,span,BasicShape::ShapeType::CYLINDER,dim, deg);
	lattice->getShape()->setRefSystem(2, darray3E{0,1,0});	
	//Set geometry
	lattice->setGeometry(&mimmo0);

// 	//Set release Info
// 	lattice->setReleaseInfo(true);

	//Set Input with Init Displacements
	int ndeg = lattice->getNDeg();
	dvecarr3E displ(ndeg, darray3E{0,0,0});
	time_t Time = time(NULL);
	srand(Time);
	ivector1D nnn = lattice->getDimension();
	for (int i=0; i<ndeg; i++){
			int l1,l2,l3;
			lattice->accessPointIndex(i,l1,l2,l3);
			if(l1>0){
				displ[i][0] = 0.5;
			//	displ[i][0] = 0.15*( (double) (rand()) / RAND_MAX );
			}
// 			if(l2 == nnn[1]-1){
// 				displ[i][0] = 1.0;
// 			}
			
	//	}	
	}
	
	
	
 	lattice->setDisplacements(displ);
	
// 	
// 	lattice->plotGrid("./", "lattice_pipe", 0, false, false);
// 	lattice->plotGrid("./", "lattice_pipe", 1, false, true);
// 	exit(1);
//********************************************************************************************	
	//CREATING INPUT	
// 	InputDoF* input = new InputDoF(ndeg, displ);
// 	string file = "input.txt";
// 	InputDoF* input0 = new InputDoF(file);
// 
// 	cout << "input setup done" << endl;
//********************************************************************************************	
	//CREATE APPLIER
	cout << "applier setup" << endl;
	Apply* applier = new Apply(&mimmo0);

	lattice->addChild(applier);

	cout << "applier setup done" << endl;
//********************************************************************************************	
	//CREATE FILTER MASK
// 	Mask* mask = new Mask();
// 	dvecarr3E coords(ndeg),coords2(ndeg);
// 	for (int i=0; i<dim[0]; i++){
// 		for (int j=0; j<dim[1]; j++){
// 			for (int k=0; k<dim[2]; k++){
// 				coords[lattice->accessPointIndex(i,j,k)] = origin + lattice->getLocalPoint(i,j,k);
// 				coords2[lattice->accessPointIndex(i,j,k)] = {{0.0, 0.0, 0.0}};
// 			}
// 		}
// 	}
// //	mask->setCoords(coords2);
// 	darray3E thres;
// //	thres[0] = 0.5;
// //	thres[1] = -10.0;
// //	thres[2] = -10.0;
// 	thres[0] = -0.0;
// 	thres[1] = -0.0;
// 	thres[2] = -0.0;
// 	mask->setThresholds(thres);
// 	mask->setForward(0,false);
// 	mask->setForward(1,false);
// 	mask->setForward(2,false);
// 
// 	//set filter to lattice
// 	mask->addChild(lattice);
// 
// 	cout << "mask setup done" << endl;
//********************************************************************************************	
	//CREATE BENDER-WRAPPER
// 	Bend* bend = new Bend();
// 	bend->setCoords(coords);
// 	dvecarr3E degree(3);
// 	degree[2][0] = 2;
// 	bend->setDegree(degree);
// 	dvector3D coeffs(3, vector<vector<double> >(3) );
// 
// 	coeffs[2][0].resize(degree[2][0]+1);
// 	coeffs[2][0][0] = 0.195;
// 	coeffs[2][0][1] = -0.8;
// 	coeffs[2][0][2] = 0.8;
// 
// 	bend->setCoeffs(coeffs);
// 	//set bend to lattice
// 	bend->addChild(mask);
// 	cout << "bend setup done" << endl;
//********************************************************************************************	
	//create output
	cout << "output setup" << endl;
	OutputDoF* output = new OutputDoF();
	lattice->addChild(output);
	cout << "output setup done" << endl;
//********************************************************************************************	

//Creating ELEMENT chain
	Chain ch0;
	ivector1D chain_pos;
// 		input0->addChild(bend);
// 		cout << "input" << endl;
// 		cout << ch0.addObject(input0) << endl;
// 		input->addChild(bend);
// 		cout << "input" << endl;
// 		cout << ch0.addObject(input) << endl;

	chain_pos.push_back(ch0.addObject(output));
// 	cout << ch0.addObject(mask) << endl;
	chain_pos.push_back(ch0.addObject(applier));
	chain_pos.push_back(ch0.addObject(lattice));
// 	cout << ch0.addObject(bend) << endl;
	cout<<chain_pos<<endl;
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

	filename = "mimmo_pipe1";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

//********************************************************************************************	
	//clean up & exit;
	delete lattice, applier, output;

	lattice = NULL;
	applier = NULL;
	output 	= NULL;

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
