/*---------------------------------------------------------------------------*\
 *
 *  mimmino
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  mimmino is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmino is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmino. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "bitpit.hpp"
#include "CreateSeedsOnSurface.hpp"
#include "MimmoGeometry.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmino;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test0009(MimmoObject * geo) {

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
		/**<Calling mimmino Test routines*/


		MimmoGeometry * mimmino0 = new MimmoGeometry();
		mimmino0->setRead(true);
		mimmino0->setReadDir("geo_data");
		mimmino0->setReadFileType(FileType::STL);
		mimmino0->setReadFilename("drivAerBin2");
		mimmino0->setWrite(true);
		mimmino0->setWriteDir("./");
		mimmino0->setWriteFileType(FileType::STL);
		mimmino0->setWriteFilename("root1");
		mimmino0->setBuildBvTree(true);
		mimmino0->execute();

		
		test0009(mimmino0->getGeometry());
		
		delete mimmino0;
		mimmino0 = NULL;

		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

