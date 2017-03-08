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
#include "CreateSeedsOnSurface.hpp"
#include "MimmoGeometry.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00016(MimmoObject * geo) {

	steady_clock::time_point t1,t2;
	duration<double> time_span;
	
	CreateSeedsOnSurface * cso =  new CreateSeedsOnSurface();
	
	cso->setGeometry(geo);
	cso->setNPoints(20);
	
 	cso->setEngineENUM(CSeedSurf::LEVELSET);
	cso->setMassCenterAsSeed(true);

	cout << "execution LEVELSET engine start" << endl;
	t1 = steady_clock::now();
 	cso->solve(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "LEVEL SET engine took " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	
	std::cout<<"minimum distance allowed is "<<cso->getMinDistance()<<std::endl;
 	cso->plotCloud("./", "distLS", 0, false);
 	
	
	cso->setEngineENUM(CSeedSurf::CARTESIANGRID);
	cso->setMassCenterAsSeed(true);

	cout << "execution CARTESIANGRID engine start" << endl;
	t1 = steady_clock::now();
	cso->solve(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "CARTESIANGRID engine took " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cso->plotCloud("./", "distCG", 0, false);
	
	cso->setRandomFixed(12141516);
	cso->setEngineENUM(CSeedSurf::RANDOM);
	cso->setMassCenterAsSeed(true);
	cso->setPlotInExecution(true);	
	
	cout << "execution RANDOM engine start" << endl;
	t1 = steady_clock::now();
	cso->solve(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "RANDOM engine took " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cso->plotCloud("./", "distRM", 0, false);
	
	delete cso;
	cso = NULL;
		
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

		
		test00016(mimmo0->getGeometry());
		
		delete mimmo0;
		mimmo0 = NULL;

		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

