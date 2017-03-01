/*---------------------------------------------------------------------------*\
 *
 *  mimmino
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  mimmino is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  mimmino is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with mimmino. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include "bitpit.hpp"
#include "MiMMO.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;


// #include <chrono>
// 
// using namespace std::chrono;
// using namespace std::placeholders;


// =================================================================================== //

// void read_Dictionary( bitpit::Config::Section & slot, std::unordered_map<std::string, std::unique_ptr<BaseManipulation > >  & mapInst, std::unordered_map<std::string, BaseManipulation * >  & mapConn, bool debug) {
// 	
// 	std::unordered_map<std::string, int> candidateMap;
// 	
// 	candidateMap.insert(std::make_pair("MiMMO.Geometry", 0));
// 	candidateMap.insert(std::make_pair("MiMMO.MultiApply", 1));
// 	candidateMap.insert(std::make_pair("MiMMO.MRBF", 2));
// 	candidateMap.insert(std::make_pair("MiMMiNO.CreateSeedsOnSurface", 3));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SelectionByBox", 4));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SelectionByCylinder", 5));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SelectionBySphere", 6));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SelectionByMapping", 7));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SelectionByPID", 8));
// 	candidateMap.insert(std::make_pair("MiMMiNO.RCPoints", 9));
// 	candidateMap.insert(std::make_pair("MiMMiNO.ReconstructScalar", 10));
// 	candidateMap.insert(std::make_pair("MiMMiNO.ReconstructVector", 11));
// 	candidateMap.insert(std::make_pair("MiMMiNO.SurfaceConstraints", 12));
// 	candidateMap.insert(std::make_pair("MiMMiNO.IOPointListFF", 13));
// 	candidateMap.insert(std::make_pair("MiMMO.ClipGeometry", 14));
// 	candidateMap.insert(std::make_pair("MiMMO.SpecularPoints", 15));
// 	candidateMap.insert(std::make_pair("MiMMO.ControlDeformMaxDistance",16));
// 	candidateMap.insert(std::make_pair("MiMMO.ControlDeformExtSurface",17));
// 	candidateMap.insert(std::make_pair("MiMMiNO.IOOptimInputWriter",18));
// 	candidateMap.insert(std::make_pair("MiMMO.SplitVectorField",19));
// 	candidateMap.insert(std::make_pair("MiMMiNO.OverlapVectorFields",20));
// 	candidateMap.insert(std::make_pair("MiMMO.StitchGeometry",21));
// 	candidateMap.insert(std::make_pair("MiMMO.IOViolationFile",22));
// 	candidateMap.insert(std::make_pair("ClassNONE", 23));
// 	
// 	for(auto & sect : slot.getSections()){
//         
// 		int fallback_counter = 0;
// 		std::string fallback_name = "ClassNONE";
// 		
// 		int classCounter = sect.second->get<int>("ClassID", fallback_counter);
// 		std::string className = sect.second->get("ClassName", fallback_name);
// 		className = bitpit::utils::trim(className);
// 		
// 		bool iomode;
// 		int  nSect;
// 		int topo;
// 		
// 		int itype = 23;
// 		if(candidateMap.count(className))	itype = candidateMap[className];
// 		//instantiate your objects and get their pointer to a unique pointer map
// 		switch(itype){
// 			
// 			case 0 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::MimmoGeometry());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 1 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::MultiApply());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 2 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::MRBF());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 				
// 			case 3 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::CreateSeedsOnSurface());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 				
// 			case 4 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SelectionByBox());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 				
// 			break;
// 				
// 			case 5 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SelectionByCylinder());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 6 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SelectionBySphere());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 7 :
// 			{
// 				topo = sect.second->get<int>("Topology", fallback_counter);
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SelectionByMapping(topo));
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 				
// 			break;
// 				
// 			case 8 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SelectionByPID());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 9 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::RCPoints());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 10 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::ReconstructScalar());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 11 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::ReconstructVector());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 12 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::SurfaceConstraints());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 			
// 			case 13:
// 
// 			{
// 				int valmode = sect.second->get<int>("IOMode");
// 				iomode = (bool)valmode;
// 				nSect = sect.second->get<int>("NSect", fallback_counter);
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::IOPointListFF(iomode, nSect));
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 				
// 			break;
// 
// 			case 14 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::ClipGeometry());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			
// 			break;
// 
// 			case 15 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::SpecularPoints());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 16 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::ControlDeformMaxDistance());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 17 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::ControlDeformExtSurface());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 
// 			case 18 :
// 			{
// 				nSect = sect.second->get<int>("NSect", fallback_counter);
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::IOOptimInputWriter(nSect));
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 19 :
// 			{
// 				topo = sect.second->get<int>("Topology", fallback_counter);
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::SplitVectorField(topo));
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 20 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmino::OverlapVectorFields());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 21 :
// 			{
// 				topo = sect.second->get<int>("Topology", fallback_counter);
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::StitchGeometry(topo));
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 			
// 			case 22 :
// 			{
// 				std::unique_ptr<BaseManipulation> temp(new mimmo::IOViolationFile());
// 				mapInst[sect.first] =  std::move(temp);
// 			}
// 			break;
// 
// 			default:
// 					//do nothing;
// 				break;
// 				
// 		}
// 			
// 			if(mapInst.count(sect.first) > 0){
// 				mapInst[sect.first]->setClassCounter(classCounter);
// 				mapInst[sect.first]->absorbSectionXML(*(sect.second.get()), sect.first);
// 			}
// 		
// 	}
// 	
// 	if(debug)	std::cout<<"Instantiated objects : "<<mapInst.size()<<std::endl;
// 		
// 	for(auto & iM : mapInst){
// 		
// 		mapConn[iM.first] = iM.second.get();
// 		
// 		if(iM.second->getName() == "MiMMiNO.IOPointListFF"){
// 			std::vector<IOPointSection * > subsection = (dynamic_cast<IOPointListFF * >(iM.second.get()))->getSections();
// 			int counter = 0;
// 			for(auto & val : subsection){
// 				std::string newname = iM.first+".Section"+std::to_string(counter);
// 				mapConn[newname] = val;
// 				++counter;
// 			}
// 		}
// 		
// 		if(iM.second->getName() == "MiMMiNO.IOOptimInputWriter"){
// 			std::vector<DeformationInfo * > subsection = (dynamic_cast<IOOptimInputWriter * >(iM.second.get()))->getMyDefInfoList();
// 			int counter = 0;
// 			for(auto & val : subsection){
// 				std::string newname = iM.first+".Section"+std::to_string(counter);
// 				mapConn[newname] = val;
// 				++counter;
// 			}
// 		}
// 		
// 	}
// 	
// 	if(debug) 	std::cout<<"Connectable objects : "<<mapConn.size()<<std::endl;
// 	
// 	
// 	//absorb connections from file if any
// 	IOConnections_CAMILO * conns = new IOConnections_CAMILO(mapConn);
// 	
// 	
// 	if(config::root.hasSection("Connections")){
// 		
// 		bitpit::Config::Section & connXML = config::root.getSection("Connections");
// 		conns->absorbConnections(connXML, true);
// 	}	
// 	
// 	delete conns; 
// 	conns=NULL;
// 	
// }

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		//reading file dictionary name from argv;
		std::string dictName;
		bool debug = false;
		
		//reading arguments
		if(argc <2)	{
			std::cout<<"Error. No source XML dictionary found"<<std::endl;
			exit(1);
		}
		
		std::vector<std::string> input(argc-1);
		for(int i=1; i<argc; ++i){
				input[i-1]= argv[i];
				input[i-1] = bitpit::utils::trim(input[i-1]);
		}
		
		if(input[0] == "--help"){
			std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
			std::cout<<""<<std::endl;
			std::cout<<"	Brief MiMMOXml_exec helper"<<std::endl;
			std::cout<<""<<std::endl;
			std::cout<<""<<std::endl;
			std::cout<<"	This is the executable command of XML Dictionaries of SOUTH Geometric Module based on         "<<std::endl;
			std::cout<<"	MiMMO/CAMiLO API."<<std::endl;
			std::cout<<"	Need to specify a path to your XML dictionary in first position after command, or "<<std::endl;
			std::cout<<"	alternatively, you can launch this helper with --help option."<<std::endl;
			std::cout<<"	A verbose execution log can be activated postponing --debug to XML dictionary path."<<std::endl;
			std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
			exit(1);
		}else{
			dictName = input[0];
		}
		
		for(int i=1; i<argc-1; ++i){
			if(input[i] == "--debug"){
				debug = true;
			}
		}
		
		//get into the mood.
		
		bitpit::config::reset("MiMMOXML");
		bitpit::config::read(dictName);
		auto globMap = getManipulatorList();

		std::unordered_map<std::string, std::unique_ptr<BaseManipulation > > mapInst;
		
		for(auto & sect: config::root.getSections()){
			
			std::string fallback_name = "ClassNONE";
			 		
			std::string className = sect.second->get("ClassName", fallback_name);
			className = bitpit::utils::trim(className);
			
			if(globMap.count(className) > 0){
				auto fptr = globMap[className];
				mapInst[sect.first] = std::move(fptr(*(sect.second.get())));
			}
			
			std::cout<<"name "<<mapInst[sect.first]->getName()<<std::endl;
			std::cout<<"counter "<<mapInst[sect.first]->getClassCounter()<<std::endl;
			std::cout<<"refresh "<<(dynamic_cast<Apply *>(mapInst[sect.first].get()))->getRefreshGeometryTrees()<<std::endl;
		}
		

// 		std::unordered_map<std::string, BaseManipulation * > mapConn;
// 		std::unordered_map<std::string, BaseManipulation * > mapInstPtr;
// 		
// 		if(debug)	std::cout<<"Reading Dictionary, instantiate objects and their connections... "<<std::endl;;
// 		
// 		read_Dictionary(config::root, mapInst, mapConn, debug);
// 		
// 		if(debug)	std::cout<<" DONE."<<std::endl;	
// 		
// 		//create execution chains;
// 		mimmo::Chain chainRead,chainRead_nPointLayout, chainExec, chainApply, chainWrite, chainReRead_nPointLayout, chainRecollectInfo, chainWriteOptimInit;
// 		
// 		if(debug)	std::cout<<"Creating Execution chains... ";
// 
// 		//activate plotInExecution on objects if debug mode.
// 		if (debug){
// 			for(auto val : mapConn){
// 				val.second->setPlotInExecution(true);
// 				val.second->setOutputPlot(".");
// 			}
// 		}
// 	
// 		for(auto &val : mapInst){
// 			mapInstPtr[val.first] = val.second.get();
// 		}
// 		
// 		for(auto val : mapInstPtr){
// 			
// 			if(!config::root.hasSection(val.first))	continue;
// 			
// 			bitpit::Config::Section & XMLclasses = config::root.getSection(val.first);
// 			std::string nclass = XMLclasses.get("ClassName");
// 			nclass = bitpit::utils::trim(nclass);
// 			
// 			std::string readMode,nsections,addinfo;
// 			int nsect;
// 			bool checkAddInfo;
// 			
// 			int itype = 0;
// 			if(nclass == "MiMMiNO.IOPointListFF") 		itype = 1;
// 			if(nclass == "MiMMO.Geometry")		 		itype = 2;
// 			if(nclass == "MiMMO.StitchGeometry")		itype = 3;
// 			if(nclass == "MiMMO.MultiApply")			itype = 4;
// 			if(nclass == "MiMMiNO.IOOptimInputWriter")	itype = 5;
// 			
// 			
// 			switch (itype){
// 				
// 				case 1 :
// 					
// 					nsections = XMLclasses.get("NSect");
// 					nsections = bitpit::utils::trim(nsections);
// 					nsect = std::atoi(nsections.c_str());
// 					
// 					if(XMLclasses.get<int>("IOMode") == 1){
// 						chainWrite.addObject(val.second);
// 						
// 						for(int i=0; i<nsect; ++i){
// 							std::string temp = val.first + ".Section"+ std::to_string(i);
// 							if(mapConn.count(temp) > 0){
// 								chainExec.addObject(mapConn[temp]);
// 							}
// 						}
// 					}else{
// 						checkAddInfo = XMLclasses.hasOption("AdditionalInfo");
// 						if(checkAddInfo){
// 							addinfo = XMLclasses.get("AdditionalInfo");
// 							addinfo = bitpit::utils::trim(addinfo);
// 						}
// 						
// 						if(checkAddInfo && addinfo == "ReReadTemplate"){
// 								chainReRead_nPointLayout.addObject(val.second);
// 								
// 								for(int i=0; i<nsect; ++i){
// 									std::string temp = val.first + ".Section"+ std::to_string(i);
// 									if(mapConn.count(temp) > 0){
// 										chainRecollectInfo.addObject(mapConn[temp]);
// 									}
// 								}
// 								
// 						}else{
// 								chainRead_nPointLayout.addObject(val.second);
// 								
// 								for(int i=0; i<nsect; ++i){
// 									std::string temp = val.first + ".Section"+ std::to_string(i);
// 									if(mapConn.count(temp) > 0){
// 										chainExec.addObject(mapConn[temp]);
// 									}
// 								}
// 							
// 						}	
// 					}	
// 
// 				break;
// 				
// 				case 2 :
// 					readMode = XMLclasses.get("ReadFlag");
// 					if( readMode == "1"){
// 						chainRead.addObject(val.second);
// 					}else{	
// 						chainWrite.addObject(val.second);
// 					}
// 				break;	
// 				
// 				case 3 :
// 						chainRead.addObject(val.second);
// 				break;	
// 				
// 				case 4:
// 					  chainApply.addObject(val.second);
// 				break;	
// 					
// 				case 5:
// 						nsections = XMLclasses.get("NSect");
// 						nsections = bitpit::utils::trim(nsections);
// 						nsect = std::atoi(nsections.c_str());
// 						
// 						chainWriteOptimInit.addObject(val.second);
// 						
// 						for(int i=0; i<nsect; ++i){
// 							std::string temp = val.first + ".Section"+ std::to_string(i);
// 							if(mapConn.count(temp) > 0){
// 								chainRecollectInfo.addObject(mapConn[temp]);
// 							};
// 						}
// 				break;	
// 				
// 				default://for every other package
// 						chainExec.addObject(val.second);
// 				break;	
// 				
// 			}
// 		} //end for loop
// 		
// 		if(debug)	std::cout<<" DONE."<<std::endl;		
// 		
// 		//Execute
// 		if(debug)	std::cout<<"Executing your workflow... "<<std::endl;
// 
// 			
// 		if (chainRead.getNObjects() > 0){
// 			if(debug)	std::cout<<"Chain of prep geos"<<std::endl;
// 			chainRead.exec(debug);
// 		}	
// 		if (chainRead_nPointLayout.getNObjects() > 0)	{
// 			if(debug)	std::cout<<"Chain of reading RBF points layout"<<std::endl;
// 			chainRead_nPointLayout.exec(debug);
// 		}	
// 		if (chainExec.getNObjects() > 0)	{
// 			if(debug)	std::cout<<"Chain of workflow execution"<<std::endl;
// 			chainExec.exec(debug);
// 		}	
// 		if (chainApply.getNObjects()> 0)	{
// 			if(debug)	std::cout<<"Chain to apply modifications on geos"<<std::endl;
// 			chainApply.exec(debug);
// 		}	
// 		if (chainWrite.getNObjects()> 0)	{
// 			if(debug)	std::cout<<"Chain write results"<<std::endl;
// 			chainWrite.exec(debug);
// 		}	
// 		if (chainReRead_nPointLayout.getNObjects()> 0)	{
// 			if(debug)	std::cout<<"Chain to read and check template of RBF points layout"<<std::endl;
// 			chainReRead_nPointLayout.exec(debug);
// 		}	
// 		if (chainRecollectInfo.getNObjects()> 0)	{
// 			if(debug)	std::cout<<"Chain to postprocess relevant parameterization data"<<std::endl;
// 			chainRecollectInfo.exec(debug);
// 		}	
// 		if (chainWriteOptimInit.getNObjects()> 0)	{
// 			if(debug)	std::cout<<"Chain to write resume file optimInit.py"<<std::endl;
// 			chainWriteOptimInit.exec(debug);
// 		}	
// 		if(debug)	std::cout<<" DONE."<<std::endl;		
// 		
		//Done, now exiting;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

