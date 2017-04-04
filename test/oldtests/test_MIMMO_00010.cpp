/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
 *
 *  mimmo is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmo is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "mimmo.hpp"
#include "ClipGeometry.hpp"
#include "SpecularPoints.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test00010() {

	
	//std::string file = "helmet"; 
	std::string file = "helmet"; 
	//Instantiation of geometry geometry Object.
	MimmoGeometry* geometry = new MimmoGeometry();
	geometry->setRead(true);
	geometry->setReadFileType(FileType::STL);
	geometry->setReadDir("./geo_data");
	geometry->setReadFilename(file);
	geometry->setWrite(false);

	ClipGeometry * clip = new ClipGeometry();
	SpecularPoints * mirror = new SpecularPoints();
	
	darray4E plane = {{0,1.0,0.0,0.0}};
	clip->setClipPlane(plane);
	clip->setInsideOut(true);
	clip->setPlotInExecution(true);
	
	mirror->setPlane(plane);
	mirror->setInsideOut(true);
	mirror->setPlotInExecution(true);
	
	dvecarr3E points(4, {{0,0,0}});
	points[0] = {{0.2,-0.5,1.2}};
	points[1] = {{0.3,-0.3,0.99}};
	points[2] = {{0.5,-0.01,1.13}};
	points[3] = {{0.8,-0.9,0.7}};

	dvector1D data(4,1);
	for(int i=0; i<4; ++i) data[i] += double(i);
	
	mirror->setCoords(points);
	mirror->setScalarData(data);
	
	bool connbuilt = true;
	std::cout<<(connbuilt && addPin(geometry, clip, PortType::M_GEOM, PortType::M_GEOM))<<std::endl;
	std::cout<<(connbuilt && addPin(geometry, mirror, PortType::M_GEOM, PortType::M_GEOM))<<std::endl;
	
	Chain c0;
	c0.addObject(geometry);
	c0.addObject(clip);
	c0.addObject(mirror);
	
	c0.exec(true);
	
	
	delete geometry;
	delete mirror;
	delete clip;
	geometry=  NULL;
	mirror = NULL;
	clip = NULL;
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

			test00010() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
