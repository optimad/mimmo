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

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test0009() {


	std::string file = "bar"; 
	
	//Instantiation of geometry geometry Object.
	MimmoGeometry* geometry = new MimmoGeometry();
	geometry->setRead(true);
	geometry->setReadFileType(STL);
	geometry->setReadDir("./geo_data");
	geometry->setReadFilename(file);
	geometry->setWrite(false);

	geometry->exec();
	OBBox * obb = new OBBox();
		
	obb->setGeometry(geometry->getGeometry());
		
	obb->execute();
	obb->plot(".", "obboxTriangulation", 0, true);
	
// 	MimmoGeometry* geometry2 = new MimmoGeometry();
// 	geometry2->setRead(true);
// 	geometry2->setReadFileType(STL);
// 	geometry2->setReadDir("./geo_data");
// 	geometry2->setReadFilename("stanfordBunnyRot");
// 	geometry2->setWrite(false);
// 		
// 	
// 	geometry->exec();
// 	geometry2->exec();
// 	
// 
// 	OBBox * obb = new OBBox();
// 	
// 	obb->setGeometry(geometry->getGeometry());
// 	
// 	obb->execute();
// 	obb->plot(".", "obboxTriangulation", 0, true);
// 
// 
// 	OBBox * obb2 = new OBBox();
// 	
// 	obb2->setGeometry(geometry2->getGeometry());
// 	
// 	obb2->execute();
// 	obb2->plot(".", "obboxTriangulationRot", 0, true);
// 	
// 	
// // 	std::cout<<obb->getAxes()<<std::endl;
// // 	std::cout<<obb2->getAxes()<<std::endl;
// // 	
// // 	std::cout<<obb->getSpan()<<std::endl;
// // 	std::cout<<obb2->getSpan()<<std::endl;
// // 	
// // 	std::cout<<obb->getOrigin()<<std::endl;
// // 	std::cout<<obb2->getOrigin()<<std::endl;
// 	
// 	
// 	
// 	
// 	dvecarr3E vertex = geometry->getGeometry()->getVertexCoords();
// 	MimmoObject * obj= new MimmoObject(3, vertex);
// 	
// 	obb->setGeometry(obj);
// 	
// 	obb->execute();
// 	obb->plot(".", "obboxPointCloud", 1, true);
// 	
 	delete obb;
// 	delete obb2;
// 	delete obj;
 	delete geometry;
// 	delete geometry2;
// 	
 	obb = NULL;
// 	obb2 = NULL;
// 	obj = NULL;
 	geometry = NULL;
// 	geometry2 = NULL;
// 	
 	return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

			test0009() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
