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
#include "ControlDeformMaxDistance.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

double test00011() {

	//Creation of MiMMO container.
	MimmoGeometry * mimmoP = new MimmoGeometry();
	
	mimmoP->setRead(true);
	mimmoP->setReadDir("geo_data");
	mimmoP->setReadFileType(FileType::STL);
	mimmoP->setReadFilename("plane");

	dvecarr3E displ(1, darray3E{{0.0,0.0,0.0}}), nodes(1,darray3E{{0.0,0.0,0.0}});

	displ[0][2] = 1.0;	
	
	//createRBF
	MRBF* mrbf = new MRBF();
	mrbf->setMode(MRBFSol::NONE);
	mrbf->setSupportRadius(-1.0);
	mrbf->setDisplacements(displ);
	mrbf->setNode(nodes);

	ControlDeformMaxDistance * volconstr = new ControlDeformMaxDistance();
	volconstr->setLimitDistance(0.8);
	
	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(mimmoP, mrbf, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(mimmoP, volconstr, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 3 : " << boolalpha << addPin(mrbf, volconstr, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(mimmoP);
	ch0.addObject(mrbf);
	ch0.addObject(volconstr);
	
	//Execution of chain
	cout << "execution start" << endl;
	steady_clock::time_point t1 = steady_clock::now();
	ch0.exec(true);
	steady_clock::time_point t2 = steady_clock::now();
	cout << "execution done" << endl;

	double result = volconstr->getViolation();	
	
	//Delete and nullify pointer
	delete mimmoP;
	delete mrbf;
	delete volconstr;
	
	mrbf 		= NULL;
	mimmoP 		= NULL;
	volconstr 	= NULL;
	
	//Print execution time
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;

	return result;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		double checkViolation = test00011() ;

		std::cout<<"Constraint violation is : "<< checkViolation<<std::endl;
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
