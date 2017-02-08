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
// #include <valgrind/callgrind.h>
using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00014() {

	//Creation of MiMMO container.
	MimmoGeometry * g1 = new MimmoGeometry();
	ProjSegmentOnSurface * psegment= new ProjSegmentOnSurface();
	MRBFStyleObj * rbf = new MRBFStyleObj();
	Apply * applier = new Apply();
	
	//settings
	g1->setRead(true);
	g1->setReadDir("geo_data");
	g1->setReadFilename("plane");
	g1->setReadFileType(FileType::STL);
	
	psegment->setSegment({{-0.5,-0.17,1.0}}, {{0.5,-0.17,1.0}});
	psegment->setBuildBvTree(true);
	psegment->setProjElementTargetNCells(100);
	psegment->setPlotInExecution(true);
	
	dvecarr3E displ(1);
	displ[0] = {{0.0,0.0,0.2}};
	rbf->setSupportRadiusValue(0.2);
	rbf->setDisplacements(displ);
	//rbf->setFunction(bitpit::RBFBasisFunction::LINEAR);
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(g1, psegment, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(g1, rbf, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 3 : " << boolalpha << addPin(psegment, rbf, PortType::M_GEOM, PortType::M_GEOM2)<<endl;
	cout << "add pin info 4 : " << boolalpha << addPin(g1, applier, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 5 : " << boolalpha << addPin(rbf, applier, PortType::M_GDISPLS, PortType::M_GDISPLS)<<endl;
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(g1);
	ch0.addObject(psegment);
	ch0.addObject(rbf);
	ch0.addObject(applier);

	duration<double> time_span;
	steady_clock::time_point t1,t2;
	//Execution of chain
	cout << "execution start" << endl;
	t1 = steady_clock::now();
		
//		CALLGRIND_START_INSTRUMENTATION;
		ch0.exec(true);
// 		CALLGRIND_STOP_INSTRUMENTATION;
// 		CALLGRIND_DUMP_STATS;
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cout << "execution done" << endl;

// 	auto map1 = stitcher->getCellDivisionMap();
// 	auto map2 = stitcher->getVertDivisionMap();
// 	
// 	std::cout<<"=========MAPCELL======================="<<endl;
// 	for(auto & v: map1){
// 		std::cout<<v.first<<'\t'<<v.second.first<<'\t'<<v.second.second<<std::endl;
// 	}
// 	
// 	std::cout<<"=========MAPVERT======================="<<endl;
// 	for(auto & v: map2){
// 		std::cout<<v.first<<'\t'<<v.second.first<<'\t'<<v.second.second<<std::endl;
// 	}
// 	
	
	g1->getGeometry()->getPatch()->write("result");
	
	//Delete and nullify pointer
	delete g1;
	delete psegment;
	delete rbf;
	delete applier;
	
	g1		= NULL;
	psegment = NULL;
	rbf = NULL;
	applier = NULL;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test00014() ;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
