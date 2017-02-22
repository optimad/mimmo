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
 \ *---------------------------------------------------------------------------*/

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

void test00018() {

	//Creation of MiMMO container.
	MimmoGeometry * g1 = new MimmoGeometry();
	MimmoGeometry * g2 = new MimmoGeometry();
	
	g1->setRead(true);
	g1->setReadDir("geo_data");
	g1->setReadFilename("sphere");
	g1->setReadFileType(FileType::STL);

	g2->setRead(true);
	g2->setReadDir("geo_data");
	g2->setReadFilename("stanfordBunny_red");
	g2->setReadFileType(FileType::STL);
	
	StitchGeometry * stitcher = new StitchGeometry(1);
	stitcher->setPlotInExecution(true);
	stitcher->setOutputPlot(".");
	
	StitchGeometry * stitcher2 = new StitchGeometry(1);
	stitcher2->setPlotInExecution(true);
	stitcher2->setOutputPlot(".");
	
	SplitVectorField* sfv1 = new SplitVectorField(1);
	sfv1->setPlotInExecution(true);
	sfv1->setOutputPlot(".");
	dvecarr3E f1(300000, {{10,0,0}});
	sfv1->setField(f1);
	
	SplitVectorField* sfv2 = new SplitVectorField(1);
	sfv2->setPlotInExecution(true);
	sfv2->setOutputPlot(".");

	dvecarr3E f2(300000, {{-5.0,0,0}});
	sfv2->setField(f2);
	
	mimmo::OverlapVectorFields * ovr = new mimmo::OverlapVectorFields();
	
	MimmoGeometry * g3 = new MimmoGeometry();
	MimmoGeometry * g4 = new MimmoGeometry();
	
	g3->setWrite(true);
	g3->setWriteDir(".");
	g3->setWriteFilename("sphere_deformed");
	g3->setWriteFileType(FileType::STVTU);
	
	g4->setWrite(true);
	g4->setWriteDir(".");
	g4->setWriteFilename("stanfordBunny_red_deformed");
	g4->setWriteFileType(FileType::STVTU);
	
	MultiApply * applier = new MultiApply();

	
	g1->setClassCounter(0);
	g2->setClassCounter(1);
	g3->setClassCounter(2);
	g4->setClassCounter(3);
	stitcher->setClassCounter(4);
	stitcher2->setClassCounter(5);
	sfv1->setClassCounter(6);
	sfv2->setClassCounter(7);
	ovr->setClassCounter(8);
	applier->setClassCounter(9);
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(g1, g3, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(g2, g4, PortType::M_GEOM, PortType::M_GEOM) << endl;
	
	cout << "add pin info 3 : " << boolalpha << addPin(g1, stitcher, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 4 : " << boolalpha << addPin(g2, stitcher, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 5 : " << boolalpha << addPin(g1, stitcher2, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	
	cout << "add pin info 6 : " << boolalpha << addPin(stitcher, sfv1, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 7 : " << boolalpha << addPin(stitcher2, sfv2, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	
	cout << "add pin info 5 : " << boolalpha << addPin(stitcher, sfv1, PortType::M_VECGEOM, PortType::M_VECGEOM)<<endl;
	cout << "add pin info 6 : " << boolalpha << addPin(stitcher2, sfv2, PortType::M_VECGEOM, PortType::M_VECGEOM)<<endl;
	
	cout << "add pin info 7 : " << boolalpha << addPin(stitcher, sfv1, PortType::M_MAPDVERT, PortType::M_MAPDVERT)<<endl;
	cout << "add pin info 8 : " << boolalpha << addPin(stitcher2, sfv2, PortType::M_MAPDVERT, PortType::M_MAPDVERT)<<endl;
	
	cout << "add pin info 9 : " << boolalpha << addPin(sfv1, ovr, PortType::M_UMGEOVFD, PortType::M_UMGEOVFD) << endl;
	cout << "add pin info 10 : " << boolalpha << addPin(sfv2, ovr, PortType::M_UMGEOVFD, PortType::M_UMGEOVFD) << endl;
	
	cout << "add pin info 11 : " << boolalpha << addPin(ovr, applier, PortType::M_UMGEOVFD, PortType::M_UMGEOVFD) << endl;
	

	cout << "set pins done" << endl;

	//Create chain
	Chain ch0,ch1;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(g1);
	ch0.addObject(g2);
	ch0.addObject(stitcher);
	ch0.addObject(stitcher2);
	ch0.addObject(sfv1);
	ch0.addObject(sfv2);
	ch0.addObject(ovr);
	ch0.addObject(applier);

	ch1.addObject(g3);
	ch1.addObject(g4);

	duration<double> time_span;
	steady_clock::time_point t1,t2;
	//Execution of chain
	cout << "execution start" << endl;
	t1 = steady_clock::now();
		
//		CALLGRIND_START_INSTRUMENTATION;
		ch0.exec(true);
		ch1.exec(true);
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
	
	//Delete and nullify pointer
	delete g1;
	delete g2;
	delete g3;
	delete g4;
	delete stitcher;
	delete stitcher2;
	delete sfv1;
	delete sfv2;
	delete applier;
	delete ovr;
	
	g1		= NULL;
	g2 		= NULL;
	g3		= NULL;
	g4 		= NULL;
	stitcher = NULL;
	stitcher2 = NULL;
	sfv1 = NULL;
	sfv2 = NULL;
	ovr=NULL;
	applier = NULL;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test00018() ;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
