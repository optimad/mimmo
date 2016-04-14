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

void test0002() {

	//Creation of MiMMO container.
	MimmoObject mimmo0;
	int	np = 0;
	int	nt = 0;
	darray3E point;
	{
		//Import STL
		STLObj stl("geo_data/catpipe.stl", true);
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
	string filename = "mimmo_0002.0000";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();


	//********************************************************************************************
	// 	//CREATING LOCAL DISPLACEMENTS LATTICE
	//Instantiation of a FFDobject of cylindrical shape
	FFDLattice* lattice = new FFDLattice();

	//Set cylindrical lattice
	darray3E origin = {-1537.5, -500.0, 3352.5};
	darray3E span;
	span[0]= 100.0;
	span[1]= 2*M_PI;
	span[2]= 1000.0;

	iarray3E dim, deg;
	dim[0] = 2;
	dim[1] = 15;
	dim[2] = 20;

	deg[0] = 1;
	deg[1] = 2;
	deg[2] = 10;

	//set lattice
	lattice->setLattice(origin,span,ShapeType::CYLINDER,dim, deg);
	lattice->setRefSystem(2, darray3E{0,-1,0});
	lattice->setGeometry(&mimmo0);

	//Set Input with Init Displacements
	int ndof = lattice->getNNodes();
	dvecarr3E displ(ndof, darray3E{0,0,0});
	double lt = 2*dim[2]/3;
	double lx;
	double a, b, c, max;
	c = 6;
	b = -15;
	a = 10;
	max = 100;
	for (int i=0; i<ndof; i++){
		int l1,l2,l3;
		int index = lattice->accessGridFromDOF(i);
		lattice->accessPointIndex(index,l1,l2,l3);
		if (l3 < int(lt) && l1 > 0){
			lx = double(l3)/lt;
			displ[i][0] = max*(c*pow(lx,5) + b*pow(lx,4) + a*pow(lx,3));
		}
		else if (l1 > 0){
			displ[i][0] = max;
		}
	}

	lattice->setDisplacements(displ);
	lattice->setDisplGlobal(false);

	//********************************************************************************************
	//********************************************************************************************
	// 	//CREATING LOCAL DISPLACEMENTS LATTICE
	//Instantiation of a FFDobject of box shape
	FFDLattice* lattice2 = new FFDLattice();

	//Set box lattice
	origin = {-110.0, 592.0, 3420.0};
	span[0]= 1000.0;
	span[1]= 400.0;
	span[2]= 150.0;

	dim[0] = 10;
	dim[1] = 2;
	dim[2] = 2;

	deg[0] = 5;
	deg[1] = 2;
	deg[2] = 2;

	//set lattice
	lattice2->setLattice(origin,span,ShapeType::CUBE,dim, deg);
	//Set geometry
	lattice2->setGeometry(&mimmo0);

	//Set Input with Init Displacements
	ndof = lattice2->getNNodes();
	displ.clear();
	displ.resize(ndof, darray3E{0,0,0});
	for (int i=0; i<ndof; i++){
		int l1,l2,l3;
		int index = lattice2->accessGridFromDOF(i);
		lattice2->accessPointIndex(index,l1,l2,l3);
		if (l1 > 3 && l1 < 6){
			if (l2 == 0){
				displ[i][1] = 100;
			}
			if (l2 == 2){
				displ[i][1] = -100;
			}
		}
	}

	lattice2->setDisplacements(displ);
	lattice2->setDisplGlobal(true);


	//********************************************************************************************
	//CREATE BENDER-WRAPPER
	//Bend* bend = new Bend();


	//********************************************************************************************

	//CREATE APPLIERS
	cout << "applier setup" << endl;
	Apply* applier = new Apply();
	applier->setGeometry(&mimmo0);
	Apply* applier2 = new Apply();
	applier2->setGeometry(&mimmo0);
	cout << "applier setup done" << endl;


	//********************************************************************************************
	//Creating PINS
	addPin(lattice, applier, &FFDLattice::getResult<dvecarr3E>, &Apply::setInput<dvecarr3E>);
	addPin(lattice2, applier2, &FFDLattice::getResult<dvecarr3E>, &Apply::setInput<dvecarr3E>);


	//********************************************************************************************

	//Creating ELEMENT chain
	Chain ch0;
	ch0.addObject(lattice);
	ch0.addObject(applier);
	ch0.addObject(lattice2);
	ch0.addObject(applier2);

	//********************************************************************************************
	//Executing CHAIN

	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();

	ch0.exec();

	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//********************************************************************************************
	//PLOT RESULTS

	lattice->plotGrid("./", "lattice_0002", 0, false, false);
	lattice->plotGrid("./", "lattice_0002", 1, false, true);

	lattice2->plotGrid("./", "lattice2_0002" , 0, false, false);
	lattice2->plotGrid("./", "lattice2_0002" , 1, false, true);

	filename = "mimmo_0002.0001";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	//********************************************************************************************
	//clean up & exit;
	delete lattice, applier;
	delete lattice2, applier2;

	lattice = NULL;
	applier = NULL;
	lattice2 = NULL;
	applier2 = NULL;

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
