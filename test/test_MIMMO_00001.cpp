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
	//placca
//	darray3E origin = {-0.01, -0.01,-0.01};
//	darray3E span;
//	span[0]= 1.02;
//	span[1]= 0.12;
//	span[2]= 0.02;
	//sphere2
	darray3E origin = {-0.0, -0.0,-0.0};
	darray3E span;
	span[0]= 1.2;
	span[1]= 1.2;
	span[2]= 1.2;
//	//cadstl
//	darray3E origin = {-0.9, -1.0, 0.03};
//	darray3E span;
//	span[0]= 5.1;
//	span[1]= 2;
//	span[2]= 1.31;

	//Set Lattice dimensions and degree
	ivector1D dim(3), deg(3);
	dim[0] = 10;
	dim[1] = 10;
	dim[2] = 10;
	deg[0] = 9;
	deg[1] = 9;
	deg[2] = 9;

//	dim[0] = 20;
//	dim[1] = 8;
//	dim[2] = 6;
//
//	deg[0] = 2;
//	deg[1] = 2;
//	deg[2] = 2;

	lattice->setMesh(origin,span,BasicShape::ShapeType::CUBE,dim, deg);

	//Set geometry
	lattice->setGeometry(&mimmo0);

	//Set release Info
	lattice->setReleaseInfo(true);


	//Set Input with Init Displacements
	int ndeg = lattice->getNDeg();
	dvecarr3E displ(ndeg);
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
		for (int j=0; j<3; j++){
			displ[i][j] = 0.15*( (double) (rand()) / RAND_MAX );
//			displ[i][j] = 0.1+0.1*j;
		}
	}
	InputDoF* input = new InputDoF(ndeg, displ);
	string file = "input/inputMIMMO_00001.txt";
	InputDoF* input0 = new InputDoF(file);


	//create applier
	Apply* applier = new Apply(&mimmo0);
	lattice->addChild(applier);

	//create filter mask
	Mask* mask = new Mask();
	dvecarr3E coords(ndeg),coords2(ndeg);
	for (int i=0; i<dim[0]; i++){
		for (int j=0; j<dim[1]; j++){
			for (int k=0; k<dim[2]; k++){
				coords[lattice->accessPointIndex(i,j,k)] = origin + lattice->getLocalPoint(i,j,k);
				coords2[lattice->accessPointIndex(i,j,k)] = {{0.0, 0.0, 0.0}};
			}
		}
	}
//	mask->setCoords(coords2);
	darray3E thres;
//	thres[0] = 0.5;
//	thres[1] = -10.0;
//	thres[2] = -10.0;
	thres[0] = -5.0;
	thres[1] = -5.0;
	thres[2] = -5.0;
	mask->setThresholds(thres);
	mask->setForward(0,false);
	mask->setForward(1,false);
	mask->setForward(2,false);
	//set filter to lattice
	mask->addChild(lattice);

	//create bend
	Bend* bend = new Bend();
	bend->setCoords(coords);
	dvecarr3E degree(3);
	degree[2][0] = 2;
	bend->setDegree(degree);
	dvector3D coeffs(3, vector<vector<double> >(3) );
	coeffs[2][0].resize(degree[2][0]+1);
	coeffs[2][0][0] = 0.195;
	coeffs[2][0][1] = -0.8;
	coeffs[2][0][2] = 0.8;
	bend->setCoeffs(coeffs);
	//set bend to lattice
	bend->addChild(mask);

	//create output
	OutputDoF* output = new OutputDoF();
	lattice->addChild(output);

	//create translation
	TranslationBox* transl = new TranslationBox();
	transl->setDirection({ {0.5, 0.25, 1} });
	transl->setTranslation(0.0);
	//set translation
	transl->addChild(lattice);

	//create rotation
	RotationBox* rotation = new RotationBox();
	rotation->setDirection({ {1.0, 1.0, 0.0} });
	rotation->setOrigin({ {0.0, 0.0, 0.0} });
//	rotation->setRotation(1.5707963267/2);
	rotation->setRotation(0.0);
	//set translation
	rotation->addChild(lattice);

	//Create chain
	Chain ch0;
	int inp;
	cout << "input zero (0) or random (1)?" << endl;
	cin >> inp;
	if (inp==0){
		input0->addChild(bend);
		cout << "add input 0" << endl;
		ch0.addObject(input0);
	}else{
		input->addChild(bend);
		cout << "add input 1" << endl;
		ch0.addObject(input);
	}
	cout << "add translation" << endl;
	ch0.addObject(transl);
	cout << "add rotation" << endl;
	ch0.addObject(rotation);
	cout << "add output" << endl;
	ch0.addObject(output);
	cout << "add mask" << endl;
	ch0.addObject(mask);
	cout << "add applier" << endl;
	ch0.addObject(applier);
	cout << "add lattice" << endl;
	ch0.addObject(lattice);
	cout << "add bend" << endl;
	ch0.addObject(bend);

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
	delete lattice, applier, mask, bend, input, input0, output;
	lattice = NULL;
	applier = NULL;
	mask 	= NULL;
	bend 	= NULL;
	input 	= NULL;
	input0 	= NULL;
	output 	= NULL;

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
