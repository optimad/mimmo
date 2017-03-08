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
\*---------------------------------------------------------------------------*/

#include "MiMMO.hpp"
#include <functional>
#include "customOperators.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

void test0007() {


	//Instantiation of mimmo geometry Object.
	MimmoGeometry* geometry = new MimmoGeometry();
	geometry->setRead(true);
	geometry->setReadFileType(FileType::OFP);
	geometry->setReadDir("./geo_data");
	geometry->setReadFilename("points");
	geometry->setWrite(true);
	geometry->setWriteFileType(FileType::OFP);
	geometry->setWriteDir(".");
	string filename = "points";
	geometry->setWriteFilename(filename);

	//Set PINS
	cout << "set pins" << endl;

	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;

	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(geometry);

	//Execution of chain
	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	ch0.exec();
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	//Delete and nullify pointer
	delete geometry;

	geometry 	= NULL;

	//Print execution time
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	return;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		test0007() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
