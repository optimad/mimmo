/*---------------------------------------------------------------------------*\
 *
 *  mimmino
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of mimmo.
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


//=================================================================================== //

void read_Dictionary( bitpit::Config::Section & slot, std::unordered_map<std::string, std::unique_ptr<BaseManipulation > >  & mapInst, std::unordered_map<std::string, BaseManipulation * >  & mapConn, Factory<BaseManipulation> & rootFactory,  bool debug) {
	
	if(debug) std::cout<< "Currently reading XML dictionary"<<std::endl;
	
	for(auto & sect : slot.getSections()){
        
		std::string fallback_name = "ClassNONE";
		
		std::string className = sect.second->get("ClassName", fallback_name);
		className = bitpit::utils::trim(className);
		std::string idstring = sect.first;
		idstring = bitpit::utils::trim(idstring);
		
		
		if(rootFactory.containsCreator(className)){
			std::unique_ptr<BaseManipulation >temp (rootFactory.create(className, *(sect.second.get())));
			mapInst[idstring] = std::move(temp);
	
			if(debug) std::cout<<"...Instantiated MiMMO block: "<<sect.first<<" of type "<<className<<std::endl;
		}else if(idstring != "Connections") {
			if(debug) std::cout<<"...Failed instantiation of "<<sect.first<<". MiMMO block of type "<<className<<" not registrated in the API"<<std::endl;
		}
	}
	
	if(debug)	std::cout<<" "<<std::endl;
	if(debug)	std::cout<<"Instantiated objects : "<<mapInst.size()<<std::endl;
		
	for(auto & iM : mapInst){
		mapConn[iM.first] = iM.second.get();
		//need to find a way to define objects inside another object
	}

	if(debug)	std::cout<<" "<<std::endl;
	if(debug) 	std::cout<<"Connectable objects : "<<mapConn.size()<<std::endl;
	

	//absorb connections from file if any
	IOConnections_MIMMO * conns = new IOConnections_MIMMO (mapConn);
	
	if(config::root.hasSection("Connections")){
		bitpit::Config::Section & connXML = config::root.getSection("Connections");
		conns->absorbConnections(connXML, debug);
	}	
	
	delete conns; 
	conns=NULL;
	
	if(debug) std::cout<< "Finished reading XML dictionary"<<std::endl;
		
}

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
			std::cout<<"	This is the executable command for running MiMMO instructions from XML Control Dictionaries"<<std::endl;
			std::cout<<" "<<std::endl;
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
		
		std::unordered_map<std::string, std::unique_ptr<BaseManipulation > > mapInst;
		std::unordered_map<std::string, BaseManipulation * > mapConn;
		std::unordered_map<std::string, BaseManipulation * > mapInstPtr;
		
		auto &factory = Factory<BaseManipulation>::instance();
		
		read_Dictionary(config::root, mapInst, mapConn, factory, debug);

		if(debug)	std::cout<<"Creating Execution chains... ";
		
		//create map of chains vs priorities.
		std::map<uint,mimmo::Chain> chainMap;
		
		for(auto &val : mapConn){
			uint priority = val.second->getPriority();
			chainMap[priority].addObject(val.second);
	
			//activate plotInExecution on objects if debug mode.
			if (debug){
				val.second->setPlotInExecution(true);
				val.second->setOutputPlot(".");
			}
		}
		if(debug)	std::cout<<" DONE."<<std::endl;		
		
		//Execute
		if(debug)	std::cout<<"Executing your workflow... "<<std::endl;

		
		for(auto &val : chainMap){
			if (val.second.getNObjects() > 0){
				if(debug)	std::cout<<"...executing Chain w/ priority "<<val.first<<std::endl;
				val.second.exec(debug);
			}
		}
		
		if(debug)	std::cout<<"Workflow DONE."<<std::endl;		
		
		//Done, now exiting;
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

