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
#include "MiMMO.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


// =================================================================================== //


void testRC() {
	GenericInput * read = new GenericInput(true);
	read->setFilename("input/input_readClouds.txt");
	
	read->execute();
	
	dvecarr3E points = read->getResult<dvecarr3E>();

	GenericOutput * write = new GenericOutput();
	write->setFilename("output_readClouds.txt");
	write->setInput(points);
	write->execute();
	
	std::cout<<"Here my points"<<std::endl;
	for(auto & val: points){
		std::cout<<val<<std::endl;
	}
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

		testRC();
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

