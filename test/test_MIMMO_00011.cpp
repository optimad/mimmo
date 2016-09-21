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
#include "ControlDeformExtSurface.hpp"
#include "IOViolationFile.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;

// =================================================================================== //

dvector1D test00011() {

	//Creation of MiMMO container.
	MimmoGeometry * mimmoP = new MimmoGeometry();
	
	mimmoP->setRead(true);
	mimmoP->setReadDir("geo_data");
	mimmoP->setReadFileType(FileType::STL);
	mimmoP->setReadFilename("plane");
	mimmoP->setBuildBvTree(true);

	dvecarr3E displ(1, darray3E{{0.0,0.0,0.0}}), nodes(1,darray3E{{-0.23,0.0,0.0}});

	displ[0][2] = 0.3;	
	
	//createRBF
	MRBF* mrbf = new MRBF();
	mrbf->setMode(MRBFSol::NONE);
	mrbf->setSupportRadius(-1.0);
	mrbf->setDisplacements(displ);
	mrbf->setNode(nodes);

	ControlDeformMaxDistance * volconstr = new ControlDeformMaxDistance();
	volconstr->setLimitDistance(0.2);
	volconstr->setPlotInExecution(true);
	
	ControlDeformExtSurface * extconstr1 = new ControlDeformExtSurface();
	extconstr1->addFile(std::make_pair("geo_data/box4plane.stl", 0));
	extconstr1->setBackgroundDetails(20);
	extconstr1->setToleranceWithinViolation(1.0E-8);
	extconstr1->setClassCounter(1);
	extconstr1->setPlotInExecution(true);
	
	ControlDeformExtSurface * extconstr2 = new ControlDeformExtSurface();
	extconstr2->addFile(std::make_pair("geo_data/P4plane1.stl", 0));
	extconstr2->addFile(std::make_pair("geo_data/P4plane2.stl", 0));
	extconstr2->setBackgroundDetails(20);
	extconstr2->setPlotInExecution(true);
	extconstr2->setClassCounter(2);
	
	Apply * applier = new Apply();
	
	MimmoGeometry * mimmoW = new MimmoGeometry();
	
	mimmoW->setRead(false);
	mimmoW->setWrite(true);
	mimmoW->setWriteDir(".");
	mimmoW->setWriteFileType(FileType::STL);
	mimmoW->setWriteFilename("plane_deformed");

	IOViolationFile * fileViolation = new IOViolationFile();
	fileViolation->setDir(".");
	fileViolation->setFileName("violationReport");

	
	//Set PINS
	cout << "set pins" << endl;

	cout << "add pin info 1 : " << boolalpha << addPin(mimmoP, mrbf, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 2 : " << boolalpha << addPin(mimmoP, volconstr, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 3 : " << boolalpha << addPin(mrbf, volconstr, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;
	cout << "add pin info 4 : " << boolalpha << addPin(mimmoP, applier, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 5 : " << boolalpha << addPin(mrbf, applier, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;
	cout << "add pin info 6 : " << boolalpha << addPin(mimmoP, mimmoW, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 7 : " << boolalpha << addPin(mimmoP, extconstr1, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 8 : " << boolalpha << addPin(mrbf, extconstr1, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;
	cout << "add pin info 9 : " << boolalpha << addPin(mimmoP, extconstr2, PortType::M_GEOM, PortType::M_GEOM) << endl;
	cout << "add pin info 10 : " << boolalpha << addPin(mrbf, extconstr2, PortType::M_GDISPLS, PortType::M_GDISPLS) << endl;
	cout << "add pin info 11 : " << boolalpha << addPin(volconstr, fileViolation, PortType::M_VIOLATION, PortType::M_VIOLATION) << endl;
	cout << "add pin info 12 : " << boolalpha << addPin(extconstr1, fileViolation, PortType::M_VIOLATION, PortType::M_VIOLATION) << endl;
	cout << "add pin info 13 : " << boolalpha << addPin(extconstr2, fileViolation, PortType::M_VIOLATION, PortType::M_VIOLATION) << endl;
	cout << "set pins done" << endl;

	//Create chain
	Chain ch0,ch1;
	cout << "add inputs and objects to the chain" << endl;
	ch0.addObject(mimmoP);
	ch0.addObject(mrbf);
	ch0.addObject(volconstr);
	ch0.addObject(extconstr1);
	ch0.addObject(extconstr2);
	ch1.addObject(applier);
	ch1.addObject(mimmoW);
	ch1.addObject(fileViolation);
	
	duration<double> time_span;
	steady_clock::time_point t1,t2;
	//Execution of chain
	cout << "execution start" << endl;
	t1 = steady_clock::now();
		ch0.exec(true);
		ch1.exec(true);
	t2 = steady_clock::now();
	time_span = duration_cast<duration<double>>(t2 - t1);
	std::cout << "MiMMO execution took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
	
	cout << "execution done" << endl;

	dvector1D resultVol(3);
	resultVol[0] = volconstr->getViolation();	
	resultVol[1] = extconstr1->getViolation();
	resultVol[2] = extconstr2->getViolation();

	//Delete and nullify pointer
	delete mimmoP;
	delete mrbf;
	delete volconstr;
	delete extconstr1;
	delete extconstr2;
	delete mimmoW;
	delete applier;
	delete fileViolation;
	
	mrbf 		= NULL;
	mimmoP 		= NULL;
	volconstr 	= NULL;
	extconstr1  = NULL;
	extconstr2  = NULL; 
	applier = NULL;
	mimmoW = NULL;
	fileViolation = NULL;
	
	return resultVol;

}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling MiMMO Test routines*/

		dvector1D checkViolation = test00011() ;

// 		std::cout<<"++++Constraint violation report ++++"<<std::endl;
// 		std::cout<<""<<std::endl;
// 		std::cout<<"MaxDistance Control gets              : "<< checkViolation[0]<<std::endl;
// 		std::cout<<"ExtSurface  Control w/ box gets       : "<< checkViolation[1]<<std::endl;
// 		std::cout<<"ExtSurface  Control w/ planes gets    : "<< checkViolation[2]<<std::endl;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
}
