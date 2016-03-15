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

void test0003() {

	//Creation of MiMMO container.
	MimmoObject mimmo0;

	int		np,	nt;
	darray3E point;
	{
		//Import STL
		STLObj stl("geo_data/ball.stl", true);
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
	string filename = "mimmo_ball0";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

//********************************************************************************************	
// 	//CREATING LATTICE
	//Instantiation of a FFDobject of spherical shape 
	FFDLattice* lattice = new FFDLattice();

	//Set cylindrical lattice
	darray3E origin = {0.0, 0.0,0.0};
	darray3E span;
	span[0]= 3.01;
	span[1]= 2*M_PI;
	span[2]= M_PI;

	ivector1D dim(3), deg(3);
	dim[0] = 4;
	dim[1] = 31;
	dim[2] = 12;
// 	dim[0] = 3;
// 	dim[1] = 5;
// 	dim[2] = 3;
	
	deg[0] = 3;
	deg[1] = 5;
	deg[2] = 3;

	//set lattice
	lattice->setMesh(origin,span,BasicShape::ShapeType::SPHERE,dim, deg);
	lattice->setRefSystem(2, darray3E{0,1,0});
//	lattice->setInfLimits(0.75*M_PI,2);
//	lattice->setInfLimits(std::atan(1.0),1);
//	lattice->setCoordType(BasicShape::CoordType::CLAMPED, 2);
	//Set geometry
	lattice->setGeometry(&mimmo0);
	
// 	//Set release Info
// 	lattice->setReleaseInfo(true);

	//Set Input with Init Displacements
	int ndeg = lattice->getNDeg();
	dvecarr3E displ(ndeg, darray3E{0,0,0});
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
			int l1,l2,l3;
			int index = lattice->accessGridFromDOF(i);
			lattice->accessPointIndex(index,l1,l2,l3);
			if(l1>0){
// 				displ[i][0] = 0.5;
// 				displ[i][0] = 0.5*( (double) (rand()) / RAND_MAX );
// 				displ[i][2] = 0.5*( (double) (rand()) / RAND_MAX );
				
// 				displ[i][0] = 1.0*( (double) (rand()) / RAND_MAX );
			}
			
		}	
	
	for(int k=1; k<dim[2]-1; ++k){
		int indGrid  = lattice->accessPointIndex(dim[0]-1,0,k);
		int indDof = lattice->accessDOFFromGrid(indGrid);
		displ[indDof][0] = 1.0;
	}
	
	lattice->setDisplacements(displ);
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
//********************************************************************************************	
	//CREATE BENDER-WRAPPER
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

	chain_pos.push_back(ch0.addObject(output));
	chain_pos.push_back(ch0.addObject(applier));
	chain_pos.push_back(ch0.addObject(lattice));
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

	lattice->plotGrid("./", "lattice_ball", 0, false, false);
	lattice->plotGrid("./", "lattice_ball", 1, false, true);

	filename = "mimmo_ball1";
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
	  test0003() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif

	return(0);
}
