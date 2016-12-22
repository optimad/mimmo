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

void test00013() {

	//Creation of MiMMO container.
	MultipleMimmoGeometries * g1 = new MultipleMimmoGeometries(1,false);
	
	g1->setRead(true);
	g1->setAddReadFile("geo_data", "sphere", FileType::STL);
	g1->setAddReadFile("geo_data", "stanfordBunny_red", FileType::STL);
	
	MultipleMimmoGeometries * g2 = new MultipleMimmoGeometries(1,true);
	
	g2->setWrite(true);
	g2->setAddWriteFile(".", "sphere", FileType::STVTU);
	g2->setAddWriteFile(".", "stanfordBunny_red", FileType::STVTU);
	
	StitchGeometry * stitcher = new StitchGeometry(1);
	stitcher->setPlotInExecution(true);
	stitcher->setOutputPlot(".");
	
	SplitVectorField* sfv = new SplitVectorField(1);
	sfv->setPlotInExecution(true);
	sfv->setOutputPlot(".");
	dvecarr3E f1(300000, {{1,2,3}});
	sfv->setField(f1);
	
	SplitScalarField* sfs = new SplitScalarField(1);
	sfs->setPlotInExecution(true);
	sfs->setOutputPlot(".");

	dvector1D f2(700000, 3.5);
	sfs->setField(f2);
	
	
	
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(g1, g2, PortType::M_VECGEOM, PortType::M_VECGEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(g1, stitcher, PortType::M_VECGEOM, PortType::M_VECGEOM)<<endl;
	cout << "add pin info 3 : " << boolalpha << addPin(stitcher, sfv, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 4 : " << boolalpha << addPin(stitcher, sfs, PortType::M_GEOM, PortType::M_GEOM)<<endl;
	cout << "add pin info 5 : " << boolalpha << addPin(g1, sfv, PortType::M_VECGEOM, PortType::M_VECGEOM)<<endl;
	cout << "add pin info 6 : " << boolalpha << addPin(g1, sfs, PortType::M_VECGEOM, PortType::M_VECGEOM)<<endl;
	cout << "add pin info 7 : " << boolalpha << addPin(stitcher, sfv, PortType::M_MAPDVERT, PortType::M_MAPDVERT)<<endl;
	cout << "add pin info 8 : " << boolalpha << addPin(stitcher, sfs, PortType::M_MAPDCELL, PortType::M_MAPDCELL)<<endl;
	
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(g1);
	ch0.addObject(g2);
	ch0.addObject(stitcher);
	ch0.addObject(sfv);
	ch0.addObject(sfs);
	
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
	
	//Delete and nullify pointer
	delete g1;
	delete g2;
	delete stitcher;
	delete sfs;
	delete sfv;
	
	g1		= NULL;
	g2 		= NULL;
	stitcher = NULL;
	sfs = NULL;
	sfv = NULL;
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
