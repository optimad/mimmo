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

// =================================================================================== //

void test0001() {

	//Creation of MiMMO container.
	MimmoObject mimmo0;

	int		np,	nt;
	darray3E point;
	{
		//Import STL
		STLObj stl("placca.stl", true);
		dvector2D V,N;
		ivector2D T;
		stl.load(np, nt, V, N, T);
		
		for (long ip=0; ip<np; ip++){
			point = conArray<double,3>(V[ip]);
			mimmo0.setVertex(ip, point);
		}
		mimmo0.setConnectivity(&T);
	}

	string filename = "mimmo0";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	//Instantiation of a FFDobject (and Input object).
	FFDLatticeBox* lattice = new FFDLatticeBox();

	//Set lattice
	darray3E origin = {{-0.5, -0.5, -0.5}};
	darray3E span;
	span.fill(1);
	//Set Lattice
	ivector1D dim(3,7), deg(3);
	lattice->setMesh(origin, span[0], span[1], span[2], dim[0], dim[1], dim[2]);

	deg[0] = dim[0]-3;
	deg[1] = dim[1]-3;
	deg[2] = dim[2]-3;
	//Set number of nodes (and degrees of curves)
	lattice->setDimension(dim, deg);

	//Set geometry
	lattice->setGeometry(&mimmo0);

	//Init Displacements
	int ndeg = lattice->getNDeg();
	dvecarr3E displ(ndeg);
	time_t Time = time(NULL);
	srand(Time);
	for (int i=0; i<ndeg; i++){
		for (int j=0; j<3; j++){
			displ[i][j] = 0.5*( (double) (rand()) / RAND_MAX );
		}
	}
	lattice->setDisplacements(displ);

	//create applier
	Apply* applier = new Apply(&mimmo0, lattice);

	//create filter mask
	MaskFilter* mask = new MaskFilter();
	dvecarr3E coords(ndeg);
	for (int i=0; i<dim[0]; i++){
		for (int j=0; j<dim[1]; j++){
			for (int k=0; k<dim[2]; k++){
				coords[lattice->accessPointData(i,j,k)] = origin + lattice->getGridPoint(i,j,k);
			}
		}
	}
	mask->setCoords(coords);
	darray3E thres;
	thres.fill(-0.25);
	mask->setThresholds(thres);
	mask->setForward(true);

	//set filter to lattice
	lattice->setFilter(mask);

	//Create execution chain
	vector<BaseManipulation*> chain;
	chain.push_back(lattice);
	chain.push_back(applier);

	for (int i=0; i<chain.size(); i++){
		chain[i]->exec();
	}

	//Plot results
	lattice->plotGrid("./", "lattice", 0, false, false);
	lattice->plotGrid("./", "lattice", 1, false, true);
	filename = "mimmo1";
	mimmo0.m_geometry->setName(filename);
	mimmo0.m_geometry->write();

	delete lattice, applier;
	lattice = NULL;
	applier = NULL;


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
