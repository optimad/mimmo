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

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00013() {

	//Creation of MiMMO container.
	MultipleMimmoGeometries * g1 = new MultipleMimmoGeometries(1,FileType::STL);
	
	g1->setRead(true);
	g1->setAddReadAbsolutePathFile("geo_data", "sphere");
	g1->setAddReadAbsolutePathFile("geo_data", "stanfordBunny_red");
	
	MultipleMimmoGeometries * g2 = new MultipleMimmoGeometries(1,FileType::STVTU);
	
	g2->setWrite(true);
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(g1, g2, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(g1, g2, PortType::M_MAPGEOM, PortType::M_MAPGEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(g1, g2, PortType::M_FINFO, PortType::M_FINFO2) << endl;
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(g1);
	ch0.addObject(g2);
	
	duration<double> time_span;
	steady_clock::time_point t1,t2;
	//Execution of chain
	cout << "execution start" << endl;
	t1 = steady_clock::now();
		ch0.exec(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cout << "execution done" << endl;

	//Delete and nullify pointer
	delete g1;
	delete g2;
	
	g1		= NULL;
	g2 		= NULL;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test00013() ;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
