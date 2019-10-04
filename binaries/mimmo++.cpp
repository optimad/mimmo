/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
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

#include "bitpit.hpp"
#include "mimmo.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;

/*!
 * Global logger of the process
 */
static bitpit::Logger * mimmo_log;

/*!
 * Global xml dictionary parser of the process
 */
static std::unique_ptr<bitpit::ConfigParser> mimmo_parser;


/*!
 * \enum Verbose
 * \brief set verbosity of message returned by mimmo++ execution
 */
enum class Verbose{
    QUIET=0 /**< no info returned*/,
    NORMAL=1 /**< medium info returned*/,
    FULL=2 /**< full info returned*/
};

//=================================================================================== //
/*!
 * \struct InfoMimmoPP
 * \brief database of essential information absorbed from mimmo++ custom arguments
 *
 * The class store from argv string of mimmo++ main the following information:
 *  - ditcname: (string) name of the target xml dictionary
 *  - vconsole: (enum VERBOSE) type of message verbosity returned by mimmo++ on console.
 *  - vlog: (enum VERBOSE) type of message verbosity returned by mimmo++ on log file.
 *  - optres: (bool) if true, return partial results of mimmo++ execution, i.e. all optional results of every block involved in the execution
 *  - optres_path: (string) specify path to save optional results of execution. Meaningful only if optres is active
 */
struct InfoMimmoPP{

    std::string dictName;       /**< target name of xml dictionary */
    Verbose vconsole;           /**< type of console verbosity.*/
    Verbose vlog;               /**< type of log file verbosity */
    bool optres;                /**< boolean to activate writing of execution optional results */
    bool expert;                /**< boolean to override mandatory ports checking */
    std::string optres_path;    /**< path to store optional results */

    /*! Base constructor*/
    InfoMimmoPP(){
        dictName    = "mimmo.xml";
        vlog        = Verbose::FULL;
        vconsole    = Verbose::NORMAL;
        optres      = false;
        optres_path = ".";
        expert      = false;
    }
    /*! Destructor */
    ~InfoMimmoPP(){};

    /*! Copy constructor */
    InfoMimmoPP(const InfoMimmoPP & other){
        *this = other;
    }

    /*!Assignement operator */
    InfoMimmoPP & operator=(const InfoMimmoPP & other){
        dictName = other.dictName;
        vlog = other.vlog;
        vconsole = other.vconsole;
        optres = other.optres;
        optres_path = other.optres_path;
        expert = other.expert;
        return *this;
    }
};
//=================================================================================== //
/*!
 * Read arguments argv of mimmo++ main
 * \return InfoMimmoPP structure filled;
 */
InfoMimmoPP readArguments(std::vector<std::string> & sanitized_args ){
    //reading arguments
    std::size_t argsize = sanitized_args.size();
    if(argsize < 2) {
        std::cout<<"Error. Not enough arguments found launching mimmo++"<<std::endl;
        std::cout<<"Please run mimmo++ --help or -h for a brief guide on how to use it"<<std::endl;
        exit(1);
    }

    //basic trimming.
    std::unordered_set<std::string> input;
    for(std::size_t i=1; i<argsize; ++i){
        input.insert(bitpit::utils::string::trim(sanitized_args[i]));
    }

    if(input.count("--help") || input.count("-h")){
        std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<"    Brief mimmo++ helper, version 1.0"<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<""<<std::endl;
        std::cout<<"    This is the executable command for running mimmo instructions from XML Control Dictionaries"<<std::endl;
        std::cout<<"    The command needs mandatorily a XML dictionary to run. It can return execution info on   "<<std::endl;
        std::cout<<"    both console(screen) and external log file. As further debug info, it can plot optional    "<<std::endl;
        std::cout<<"    results of its execution.   "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    The full list of options for running mimmo++ are: "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --help,-h                                       : print this helper "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --dictionary,-d=<dictionary name>               : full path to the target xml dictionary. Mandatory. "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --console-verbosity,-cv=<quiet/normal/full>     : print info on mimmo++ execution on console(screen).        "<<std::endl;
        std::cout<<"                                                    full is meant for full debug message verbosity, normal for "<<std::endl;
        std::cout<<"                                                    a medium verbosity, quiet shut off messaging on console.   "<<std::endl;
        std::cout<<"                                                    Default verbosity is medium.                               "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --log-verbosity,-lv=<quiet/normal/full>         : print info on mimmo++ execution external file mimmo.log.   "<<std::endl;
        std::cout<<"                                                    full is meant for full debug message verbosity, normal for "<<std::endl;
        std::cout<<"                                                    a medium verbosity, quiet shut off messaging on file.      "<<std::endl;
        std::cout<<"                                                    Default verbosity is full.                                 "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --optional-results,-or=yes                      : activate plot of optional/debug results of execution       "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --optional-results-path,-orp=<path>             : specify directory to store optional results.               "<<std::endl;
        std::cout<<"                                                    Default directory is the current one ./                        "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    --expert,-e=yes                                 : override mandatory ports connection checking.              "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<" "<<std::endl;
        std::cout<<"    For any problem, bug and malfunction please contact mimmo developers.                       "<<std::endl;
        std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
        exit(1);
    }

    std::unordered_map<int, std::string> keymap;
    int nkeys = 6;
    keymap[0] = "--dictionary=";
    keymap[1] = "--log-verbosity=";
    keymap[2] = "--console-verbosity=";
    keymap[3] = "--optional-results=";
    keymap[4] = "--optional-results-path=";
    keymap[5] = "--expert=";

    keymap[nkeys] = "-d=";
    keymap[nkeys+1] = "-lv=";
    keymap[nkeys+2] = "-cv=";
    keymap[nkeys+3] = "-or=";
    keymap[nkeys+4] = "-orp=";
    keymap[nkeys+5] = "-e=";

    keymap[2*nkeys] = "dict=";
    keymap[2*nkeys+1] = "vlog=";
    keymap[2*nkeys+2] = "vconsole=";
    keymap[2*nkeys+3] = "opt-res=";
    keymap[2*nkeys+4] = "opt-res-path=";
    keymap[2*nkeys+5] = "expert=";

    std::map<int, std::string> final_map;
    //visit input list and search for each key string  in key map. If an input string positively match a key,
    // clean the entry from key, give her the key map marker and store it in the final map.

    for(auto val: input){
        std::size_t pos = std::string::npos;
        int counter=-1;
        while (pos == std::string::npos && counter <3*nkeys){
        	++counter;
        	pos = val.find(keymap[counter]);
        }

        if(pos != std::string::npos) final_map[counter%nkeys] =val.substr(pos+keymap[counter].size());
    }

    if(final_map.count(0) < 1) {
        std::cout<<"Error. Not valid xml dictionary found launching mimmo++"<<std::endl;
        std::cout<<"Please run mimmo++ --help for a brief guide on how to use it"<<std::endl;
        exit(1);
    }


    //now assign arguments to InfoMimmoPP class
    InfoMimmoPP result;
    result.dictName = final_map[0];
    if(final_map.count(4)) result.optres_path = final_map[4];
    if(final_map.count(3)) result.optres = (final_map[3]=="yes");
    if(final_map.count(5)) result.expert = (final_map[5]=="yes");

    if(final_map.count(1)){
        int check = -1 + int(final_map[1]=="quiet") + 2*int(final_map[1]=="normal") + 3*int(final_map[1]=="full");
        if(check == -1) check = 2;
        result.vlog = static_cast<Verbose>(check);
    }

    if(final_map.count(2)){
        int check = -1 + int(final_map[2]=="quiet") + 2*int(final_map[2]=="normal") + 3*int(final_map[2]=="full");
        if(check == -1) check = 1;
        result.vconsole = static_cast<Verbose>(check);
    }

    return result;
}

//=================================================================================== //
/*!
 * read a mimmoXML dictionary, made as :
 * ##xml header string
 * ##dictionary mimmoXML header with version
 * ## <B>Blocks</B> ...declaration of all executable objects <B>/Blocks</B>
 * ## <B>Connections</B> ...declaration of all links between executable objects <B>/Connections</B>
 *
 * \param[out] mapInst map of all declared and instantiated blocks
 * \param[out] mapConn map of all blocks that need to be linked/connected
 * \param[in] rootFactory reference to internal register of mimmo API
 */
void read_Dictionary(std::map<std::string, std::unique_ptr<BaseManipulation > >  & mapInst, std::unordered_map<std::string, BaseManipulation * >  & mapConn, Factory<BaseManipulation> & rootFactory) {

    mimmo_log->setPriority(bitpit::log::NORMAL);
    (*mimmo_log)<< "Currently reading XML dictionary"<<std::endl;

    if(mimmo_parser->hasSection("Blocks")){
        bitpit::Config::Section & blockXML = mimmo_parser->getSection("Blocks");
        for(auto & sect : blockXML.getSections()){

            std::string fallback_name = "ClassNONE";

            std::string className = sect.second->get("ClassName", fallback_name);
            className = bitpit::utils::string::trim(className);
            std::string idstring = sect.first;
            idstring = bitpit::utils::string::trim(idstring);


            if(rootFactory.containsCreator(className)){
                std::unique_ptr<BaseManipulation >temp (rootFactory.create(className, *(sect.second.get())));
                //class is instantiated. I can push my custom name now.
                temp->setName(sect.first);
                mapInst[idstring] = std::move(temp);

                (*mimmo_log)<< "...Instantiated mimmo block: "<<sect.first<<" of type "<<className<<std::endl;
            }else {
                (*mimmo_log)<<"...Failed instantiation of "<<sect.first<<". mimmo block of type "<<className<<" not registered in the API"<<std::endl;
            }
        }

        (*mimmo_log)<<" "<<std::endl;
        (*mimmo_log)<<"Instantiated objects : "<<mapInst.size()<<std::endl;

        for(auto & iM : mapInst){
            mapConn[iM.first] = iM.second.get();
            //need to find a way to define objects inside another object
        }

        (*mimmo_log)<<" "<<std::endl;
        (*mimmo_log)<<"Connectable objects : "<<mapConn.size()<<std::endl;
        mimmo_log->setPriority(bitpit::log::DEBUG);
    }else{
        (*mimmo_log)<<"No Blocks section available in the XML dictionary"<<std::endl;
        mimmo_log->setPriority(bitpit::log::DEBUG);
    }

	//absorb connections from file if any
	std::unique_ptr<IOConnections_MIMMO> conns (new IOConnections_MIMMO (mapConn));

	if(mimmo_parser->hasSection("Connections")){
		bitpit::Config::Section & connXML = mimmo_parser->getSection("Connections");
		conns->absorbConnections(connXML, false);
	}else{
        mimmo_log->setPriority(bitpit::log::NORMAL);
        (*mimmo_log)<<"No Connections section available in the XML dictionary"<<std::endl;
        mimmo_log->setPriority(bitpit::log::DEBUG);
    }

    mimmo_log->setPriority(bitpit::log::NORMAL);
    (*mimmo_log)<< "Finished reading XML dictionary"<<std::endl;
    mimmo_log->setPriority(bitpit::log::DEBUG);


}

// =================================================================================== //
//core of xml handler

void mimmocore(const InfoMimmoPP & info) {
        //set the logger and the verbosity of the output messages in execution
        std::string log = "mimmo";
        mimmo::setLogger(log);
        auto mimmo_log = &bitpit::log::cout(log);
        switch(int(info.vconsole)){
            case 1 :
                bitpit::log::setConsoleVerbosity((*mimmo_log), bitpit::log::Verbosity::NORMAL);
                break;
            case 2 :
                bitpit::log::setConsoleVerbosity((*mimmo_log), bitpit::log::Verbosity::DEBUG);
                break;
            case 0 :
                bitpit::log::setConsoleVerbosity((*mimmo_log), bitpit::log::Verbosity::QUIET);
                break;
            default: //never been reached
                break;
        }

        switch(int(info.vlog)){
            case 1 :
                bitpit::log::setFileVerbosity((*mimmo_log), bitpit::log::Verbosity::NORMAL);
                break;
            case 2 :
                bitpit::log::setFileVerbosity((*mimmo_log), bitpit::log::Verbosity::DEBUG);
                break;
            case 0 :
                bitpit::log::setFileVerbosity((*mimmo_log), bitpit::log::Verbosity::QUIET);
                break;
            default: //never been reached
                break;

        }

        mimmo::setExpertMode(info.expert);

        //print resume args info.
        mimmo_log->setPriority(bitpit::log::NORMAL);
        {
            std::vector<std::string> verb(3, "quiet");
            verb[1] = "normal";
            verb[2] = "full";
            std::vector<std::string> yesno(2, "no");
            yesno[1] = "yes";

            (*mimmo_log)<< "Resume of arguments in mimmo++:"<<std::endl;
            (*mimmo_log)<< "dictionary:         "<<info.dictName<<std::endl;
            (*mimmo_log)<< "console verbosity:  "<<verb[static_cast<int>(info.vconsole)]<<std::endl;
            (*mimmo_log)<< "log file verbosity: "<<verb[static_cast<int>(info.vlog)]<<std::endl;
            (*mimmo_log)<< "debug results:      "<<yesno[int(info.optres)]<<std::endl;
            (*mimmo_log)<< "debug results path: "<<info.optres_path<<std::endl;
            (*mimmo_log)<< "expert mode:        "<<yesno[int(info.expert)]<<std::endl;
            (*mimmo_log)<< " "<<std::endl;
            (*mimmo_log)<< " "<<std::endl;
        }
        mimmo_log->setPriority(bitpit::log::DEBUG);

    	//Instantiate of global parser
      	mimmo_parser = std::unique_ptr<bitpit::ConfigParser>(new bitpit::ConfigParser("mimmoXML", 1, true));
      	mimmo_parser->read(info.dictName);

		std::map<std::string, std::unique_ptr<BaseManipulation > > mapInst;
		std::unordered_map<std::string, BaseManipulation * > mapConn;
		std::unordered_map<std::string, BaseManipulation * > mapInstPtr;

		auto &factory = Factory<BaseManipulation>::instance();
		read_Dictionary(mapInst, mapConn, factory);

        mimmo_log->setPriority(bitpit::log::NORMAL);
		(*mimmo_log)<<"Creating Execution chains... ";

		//create map of chains vs priorities.
		std::map<uint,mimmo::Chain> chainMap;

		for(auto &val : mapConn){
			uint priority = val.second->getPriority();
			chainMap[priority].addObject(val.second);
		}
		(*mimmo_log)<<" DONE."<<std::endl;

		//Execute
        (*mimmo_log)<<"Executing your workflow... "<<std::endl;
        mimmo_log->setPriority(bitpit::log::DEBUG);


		for(auto &val : chainMap){
			if (val.second.getNObjects() > 0){
                mimmo_log->setPriority(bitpit::log::NORMAL);
                (*mimmo_log)<<"...executing Chain w/ priority "<<val.first<<std::endl;
                mimmo_log->setPriority(bitpit::log::DEBUG);
                val.second.setPlotDebugResults(info.optres);
                val.second.setOutputDebugResults(info.optres_path);
				val.second.exec(true);
			}
		}

		mimmo_log->setPriority(bitpit::log::NORMAL);
		(*mimmo_log)<<"Workflow DONE."<<std::endl;
		//Done, now exiting;
        mimmo_log->setPriority(bitpit::log::DEBUG);
}


//main program
int main( int argc, char *argv[] ) {

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);

    {
        #endif
    	//Instantiate of global logger
    	mimmo_log = &bitpit::log::cout("mimmo");

        //sanitize argv
        std::vector<std::string> s_argv;
        s_argv.reserve(argc);
        for(int i=0; i<argc; ++i){
            char * sstr = strdup(argv[i]);
            if(sstr != 0){
                s_argv.push_back(std::string(sstr));
            }
            free(sstr);
        }


        try{
            InfoMimmoPP info = readArguments(s_argv);
            mimmocore(info);
        }
        catch(std::exception & e){
        	(*mimmo_log)<<"mimmo++ exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
    }

    MPI_Finalize();
    #endif

    return 0;
}
