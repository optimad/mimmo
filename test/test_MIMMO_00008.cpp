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
#include "customOperators.hpp"
//#include "callgrind.h"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test0008() {


	//Instantiation of geometry geometry Object.
	MimmoGeometry* geometry = new MimmoGeometry();
	geometry->setRead(true);
	geometry->setReadFileType(STL);
	geometry->setReadDir("./geo_data");
	geometry->setReadFilename("sphere2");
	geometry->setWrite(false);

	geometry->exec();
	
	geometry->getGeometry()->buildBvTree(1);
 	geometry->getGeometry()->buildKdTree();
	
	Cylinder * cyl = new Cylinder();
	cyl->setOrigin({{0.5,0.5,0}});
	cyl->setSpan(0.6,2*M_PI,2);
	
	steady_clock::time_point t1,t2;
	duration<double> time_span;
	
	t1 = steady_clock::now();
	livector1D testSimplex1 = cyl->excludeGeometry(geometry->getGeometry());
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout<<"done simplex w bvtree "<<time_span.count() << " seconds."<<std::endl;
	
	t1 = steady_clock::now();
	livector1D testSimplex2 = cyl->excludeGeometry(geometry->getGeometry()->getPatch());
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout<<"done simplex w one-by-one search in  "<<time_span.count() << " seconds."<<std::endl;
	
	//check if they are the same and return 
	bool check = testSimplex1.size() == testSimplex2.size();
	if(check){
		for (auto & val : testSimplex1){
			check = check && (testSimplex2.end() != std::find(testSimplex2.begin(), testSimplex2.end(), val));
		}
	}

	std::cout<<testSimplex1.size()<<std::endl;
	std::cout<<testSimplex2.size()<<std::endl;
	std::cout<<"Tessellation matching "<<check<<std::endl;
	//print something;
	MimmoObject * obj1 = new MimmoObject(1);
	//fill obj 1 to plot;
	
	{ 
		livector1D idVert = geometry->getGeometry()->getVertexFromCellList(testSimplex1);
		
		for(auto &idV :idVert){
			darray3E coords = geometry->getGeometry()->getVertexCoords(idV);
			obj1->addVertex(coords, idV);
		}
		
		for(auto &idC :testSimplex1){
			bitpit::Cell & cc =  geometry->getGeometry()->getPatch()->getCell(idC);
			livector1D conn(cc.getVertexCount()) ;
			for(int i=0; i<cc.getVertexCount(); ++i){
				conn[i] = cc.getVertex(i);
			}
			obj1->addConnectedCell(conn, cc.getType(), idC);
		}
		
		obj1->getPatch()->write();
	}
	
	
	delete obj1;
	obj1 = NULL;

	
	t1 = steady_clock::now();
	//CALLGRIND_START_INSTRUMENTATION;
	livector1D testNode1 = cyl->includeCloudPoints(geometry->getGeometry());
	//CALLGRIND_STOP_INSTRUMENTATION;
	//CALLGRIND_DUMP_STATS;
	
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout<<"done vertex w kdtree search in  "<<time_span.count() << " seconds."<<std::endl;
	
	t1 = steady_clock::now();
	livector1D testNode2 = cyl->includeCloudPoints(geometry->getGeometry()->getPatch());
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout<<"done vertex w one-by-one search in  "<<time_span.count() << " seconds."<<std::endl;

	
	
	
	std::cout<<testNode1.size()<<std::endl;
	std::cout<<testNode2.size()<<std::endl;
	
	//check if they are the same and return 
	check = testNode1.size() == testNode2.size();
	
	if(check){
		for (auto & val : testNode1){
			check = check && (testNode2.end() != std::find(testNode2.begin(), testNode2.end(), val));
		}
	}
	
	std::cout<<"point cloud 1/2 matching "<<check<<std::endl;
	
	
	{
		dvecarr3E points(testNode1.size());
		ivector1D conn(testNode1.size());
		int counter = 0;
		for(auto & idV : testNode1){
			points[counter] = geometry->getGeometry()->getVertexCoords(idV);
			conn[counter] = counter;
			++counter;
		}
		
		bitpit::VTKUnstructuredGrid   output( "./", "Cloud1", bitpit::VTKElementType::VERTEX, points, conn);
		output.write() ;
		
	}
	
	//Delete and nullify pointer
	delete geometry;
	delete cyl;
	
	geometry 	= NULL;
	cyl = NULL;
	
	return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test0008() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
