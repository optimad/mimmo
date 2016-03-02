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

void test0002() {

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
	FFDLattice* lattice = new FFDLattice();

	//Set lattice
	darray3E origin = {{0.5, 0.05, 0.0}};
	dmatrix32E limits;
	limits[0][0]= limits[1][0] = limits[2][0] = -0.01;
	limits[0][1]= 1.01;
	limits[1][1]= 0.11;
	limits[2][1]= 0.01;
	//Set Lattice
	ivector1D dim(3), deg(3);
	dim[0] = 21;
	dim[1] = 7;
	dim[2] = 7;
	
	deg[0] = 2;
	deg[1] = 2;
	deg[2] = 2;

	lattice->setMesh(origin, limits,BasicShape::ShapeType::CUBE,dim, deg);

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
			displ[i][j] = 0.025*( (double) (rand()) / RAND_MAX );
		}
	}
	lattice->setDisplacements(displ);
	
	//do not deform;
	lattice->execute();
	lattice->plotGrid("./", "lattice", 0, false, false);
	lattice->plotGrid("./", "lattice", 1, false, true);

	string filename2 = "mimmo1";
	mimmo0.m_geometry->setName(filename2);
	mimmo0.m_geometry->write();
	
	return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

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
}
