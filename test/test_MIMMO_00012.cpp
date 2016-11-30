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

void test00012() {

	//Creation of MiMMO container.
	MimmoGeometry * g1 = new MimmoGeometry();
	
	g1->setRead(true);
	g1->setReadDir("geo_data");
	g1->setReadFileType(FileType::STL);
	g1->setReadFilename("plane");
	
	MimmoGeometry * g2 = new MimmoGeometry();
	
	g2->setWrite(true);
	g2->setWriteDir("./");
	g2->setWriteFileType(FileType::PCVTU);
	g2->setWriteFilename("PointCloudPlane");
	
	MimmoGeometry * g3 = new MimmoGeometry();
	
	g3->setRead(true);
	g3->setReadDir("./");
	g3->setReadFileType(FileType::PCVTU);
	g3->setReadFilename("PointCloudPlane");
	
	MimmoGeometry * g4 = new MimmoGeometry();
	
	g4->setWrite(true);
	g4->setWriteDir("./");
	g4->setWriteFileType(FileType::OFP);
	g4->setWriteFilename("PointCloudPlane");
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(g1, g2, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 1 : " << boolalpha << addPin(g3, g4, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0,ch1;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(g1);
	ch0.addObject(g2);
	ch1.addObject(g3);
	ch1.addObject(g4);
	
	duration<double> time_span;
	steady_clock::time_point t1,t2;
	//Execution of chain
	cout << "execution start" << endl;
	t1 = steady_clock::now();
		ch0.exec(true);
		ch1.exec(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cout << "execution done" << endl;

	//Delete and nullify pointer
	delete g1;
	delete g2;
	delete g3;
	delete g4;
	
	g1		= NULL;
	g2 		= NULL;
	g3 		= NULL;
	g4		= NULL;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test00012() ;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
