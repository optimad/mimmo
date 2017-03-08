/*---------------------------------------------------------------------------*\
 * 
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
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
 \ *---------------------------------------------------------------------------*/

#include "bitpit.hpp"
#include "MiMMO.hpp"
#include "customOperators.hpp"
#include "MeshSelection.hpp"
#include "BvTree.hpp"
#include <functional>

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00015(MimmoObject * selection, MimmoObject * selection2) {

	steady_clock::time_point t1,t2;
	duration<double> time_span;
	
	SelectionByBox  	* boxSel = new SelectionByBox();
	SelectionBySphere   * sphSel = new SelectionBySphere();
	SelectionByMapping  * mapSel = new SelectionByMapping();
	SelectionByPID		* pidSel = new SelectionByPID();

	MimmoGeometry * service = new MimmoGeometry();
	service->setWrite(true);
	service->setRead(false);
	service->setWriteDir("./");
	service->setWriteFileType(FileType::STVTU);
	
		
	boxSel->setOrigin({{-0.5,-0.5,0.2}});
	boxSel->setSpan(0.6,0.6,0.6);
	boxSel->setGeometry(selection);
	boxSel->setDual(true);
	
	sphSel->setOrigin({{-0.5, 0.5,0.2}});
	sphSel->setSpan(0.34, 2*M_PI, M_PI);
	sphSel->setGeometry(selection);
	
	mapSel->setGeometry(selection);
	mapSel->addFile(std::make_pair("geo_data/presa1.stl", FileType::STL));
	mapSel->addFile(std::make_pair("geo_data/presa2.stl", FileType::STL));

	pidSel->setPID(2);
	pidSel->setPID(3);
	pidSel->setPID(4);
	pidSel->setGeometry(selection2);
	
	t1 = steady_clock::now();
	boxSel->exec();
	t2 = steady_clock::now();
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "elementary box selection execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	service->setGeometry(boxSel->getPatch());
	service->setWriteFilename("boxSelection");
	service->execute();
	
	t1 = steady_clock::now();
	sphSel->exec();
	t2 = steady_clock::now();
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "elementary sphere selection execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	service->setGeometry(sphSel->getPatch());
	service->setWriteFilename("sphereSelection");
	service->execute();
	
	t1 = steady_clock::now();
	mapSel->exec();
	t2 = steady_clock::now();
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "mapping selection execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	service->setGeometry(mapSel->getPatch());
	service->setWriteFilename("mapSelection");
	service->execute();
	
	t1 = steady_clock::now();
	mapSel->setDual(true);
	mapSel->exec();
	t2 = steady_clock::now();
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "mapping dual selection execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	service->setGeometry(mapSel->getPatch());
	service->setWriteFilename("mapDualSelection");
	service->execute();

	t1 = steady_clock::now();
	pidSel->exec();
	t2 = steady_clock::now();
	//Print execution time
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "pid selection execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	service->setGeometry(pidSel->getPatch());
	service->setWriteFilename("pidSelection");
	service->execute();
	

	delete boxSel;		boxSel = NULL;
	delete sphSel;		sphSel = NULL;
	delete mapSel;		mapSel = NULL;
	delete pidSel;		pidSel = NULL;
	delete service; 	service = NULL;
		
	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routines*/


		MimmoGeometry * mimmo0 = new MimmoGeometry();
		mimmo0->setRead(true);
		mimmo0->setReadDir("geo_data");
		mimmo0->setReadFileType(FileType::STL);
		mimmo0->setReadFilename("drivAerBin2");
		mimmo0->setWrite(true);
		mimmo0->setWriteDir("./");
		mimmo0->setWriteFileType(FileType::STL);
		mimmo0->setWriteFilename("root1");
		mimmo0->setBuildBvTree(true);
		mimmo0->execute();

		MimmoGeometry * mimmo1 = new MimmoGeometry();
		mimmo1->setRead(true);
		mimmo1->setReadDir("geo_data");
		mimmo1->setReadFileType(FileType::STL);
		mimmo1->setReadFilename("ahmed");
		mimmo1->setWrite(true);
		mimmo1->setWriteDir("./");
		mimmo1->setWriteFileType(FileType::STL);
		mimmo1->setWriteFilename("root2");
		mimmo1->setBuildBvTree(true);
		mimmo1->execute();
		mimmo1->setRead(false);
		mimmo1->setWrite(true);
		mimmo1->setWriteDir("./");
		mimmo1->setWriteFileType(FileType::NAS);
		mimmo1->setWriteFilename("root3");
		mimmo1->setBuildBvTree(false);
		mimmo1->execute();
		
		
		test00015(mimmo0->getGeometry(), mimmo1->getGeometry());
		
		delete mimmo0;
		mimmo0 = NULL;

		delete mimmo1;
		mimmo1 = NULL;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

