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
#include "MRBF.hpp"
#include "Apply.hpp"
#include "ClipGeometry.hpp"
#include "SpecularPoints.hpp"
#include "Chain.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00017(MimmoGeometry * geo) {

	ClipGeometry * clip = new ClipGeometry();
	CreateSeedsOnSurface * seeder0 = new CreateSeedsOnSurface();
	SpecularPoints * mirror = new SpecularPoints();
	MRBF * morpher = new MRBF();
	Apply * applier = new Apply();
	MimmoGeometry* writerGeo = new MimmoGeometry();
	
	MimmoGeometry * geoCopy = new MimmoGeometry();
	geoCopy->setHARDCopy(geo);
	
	//customize blocks
	//seeder0
	seeder0->setNPoints(4);
	seeder0->setClassCounter(0);
	seeder0->setEngineENUM(CSeedSurf::RANDOM);
	seeder0->setMassCenterAsSeed(true);
	seeder0->setRandomFixed(12131415);
	seeder0->setPlotInExecution(true);
	seeder0->setOutputPlot(".");
	
	clip->setClipPlane({{0.0,1.0,0.0,0.0}});	
	clip->setPlotInExecution(true);
	clip->setOutputPlot(".");
	clip->setClassCounter(1);
	clip->setGeometry(geo->getGeometry());
	
	dvecarr3E displacements(4);
	dvecarr2E limits(3);
	{
		time_t Time = time(NULL);
		srand(Time);
		for (auto & val:displacements){
			val[0] = 0.03*( (double) (rand()) / RAND_MAX );
			val[1] = 0.06*( (double) (rand()) / RAND_MAX );
			val[2] = 0.02*( (double) (rand()) / RAND_MAX );
		}
	}

	{
		limits[0][0] = -1.0;	 	limits[0][1] = 1.0;
		limits[1][0] = -0.005;	limits[1][1] = 0.05;
		limits[2][0] = 0;		limits[2][1] = 0.7;
	}
	
	mirror->setPlane({{0.0,1.0,0.0,0.0}});	
	mirror->setPlotInExecution(true);
	mirror->setOutputPlot(".");
	mirror->setClassCounter(2);
	mirror->setVectorData(displacements);
	mirror->setGeometry(geo->getGeometry());
	
	
	morpher->setSupportRadius(-1.0);
	morpher->setClassCounter(3);
	morpher->setGeometry(geo->getGeometry());
	
	applier->setClassCounter(4);
	applier->setGeometry(geo->getGeometry());
	
	writerGeo->setWriteDir(".");
	writerGeo->setWriteFileType(FileType::STL);
	writerGeo->setWriteFilename("helmet_deformed");
	writerGeo->setWrite(true);
	writerGeo->setGeometry(geo->getGeometry());
	
	//Set PINS
	cout << "set pins" << endl;
	std::cout << " add pin " <<  boolalpha << addPin(clip, seeder0, PortType::M_GEOM, PortType::M_GEOM) << std::endl;
	std::cout << " add pin " <<  boolalpha << addPin(seeder0, mirror, PortType::M_COORDS, PortType::M_COORDS) << std::endl;
	std::cout << " add pin " <<  boolalpha << addPin(mirror, morpher, PortType::M_COORDS, PortType::M_COORDS) << std::endl;
	std::cout << " add pin " <<  boolalpha << addPin(mirror, morpher, PortType::M_DISPLS, PortType::M_DISPLS) << std::endl;
	std::cout << " add pin " <<  boolalpha << addPin(morpher, applier, PortType::M_GDISPLS, PortType::M_GDISPLS) << std::endl;
	
	cout << "set pins done" << endl;
	
	// 	//Create chain
	mimmo::Chain ch0, ch1;
	ch0.addObject(clip);
	ch0.addObject(seeder0);
	ch0.addObject(mirror);
	ch0.addObject(morpher);
	ch0.addObject(applier);
	
	
	ch1.addObject(writerGeo);
	
	//Execution of chain
	cout << "execution start" << endl;
	ch0.exec(true);
	ch1.exec(true);
	cout << "execution done" << endl;
	
	delete seeder0;
	delete clip;
	delete morpher;
	delete writerGeo;
	delete mirror;
	delete applier;
	delete geoCopy;
	
	seeder0 = NULL;
	clip = NULL;
	morpher = NULL;
	writerGeo = NULL;
	applier = NULL;
	mirror = NULL;
	geoCopy = NULL;
	return;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);
	
#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmo Test routines*/


		MimmoGeometry * mimmo0 = new MimmoGeometry();
		mimmo0->setRead(true);
		mimmo0->setReadDir("geo_data");
		mimmo0->setReadFileType(FileType::STL);
		mimmo0->setReadFilename("helmet");
		mimmo0->setWrite(false);
		mimmo0->setBuildBvTree(true);
		mimmo0->execute();

		
		test00017(mimmo0);
		
		delete mimmo0;
		mimmo0 = NULL;

		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

