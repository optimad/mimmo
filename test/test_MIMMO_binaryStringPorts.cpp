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
#include "IOPointListFF.hpp"
#include "RCPoints.hpp"
#include "MimmoGeometry.hpp"
#include "Chain.hpp"

using namespace std;
using namespace bitpit;
using namespace mimmo;
using namespace mimmino;


#include <chrono>

using namespace std::chrono;
using namespace std::placeholders;


// =================================================================================== //

class Checker: public BaseManipulation {
public:
	int m_c;
	std::string * m_dir = new std::string();
	std::string *m_file = new std::string();
	
	Checker(int c){
		m_c = c;

// 		std::cout<<"class "<<m_c<<" m_name "<<&m_name<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_counter "<<&m_counter<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_geometry "<<&m_geometry<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_parent "<<&m_parent<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_child "<<&m_child<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_portsType "<<&m_portsType<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_portIn"<<&m_portIn<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_portOut "<<&m_portOut<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_arePortsBuilt "<<&m_arePortsBuilt<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_active "<<&m_active<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_execPlot "<<&m_execPlot<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_outputPlot "<<&m_outputPlot<<std::endl;
// 		std::cout<<"class "<<m_c<<" sm_baseManipulationCounter "<<&sm_baseManipulationCounter<<std::endl;
// 		
// 		std::cout<<"class "<<m_c<<" m_dir "<<&m_dir<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_file "<<&m_file<<std::endl;
// 		std::cout<<"class "<<m_c<<" m_c "<<&m_c<<std::endl;
	};
	virtual ~Checker(){
		
		delete m_dir;
		delete m_file;
		m_dir = NULL;
		m_file = NULL;
		std::cout<<"called destructor of checker"<<m_counter<<std::endl;
	};
	
	void setDir( std::string  val){
		*m_dir = val;
	};
	
	void setFile(std::string  val){
		*m_file = val;
		
	};
	
	void setDir( std::string * val){
		*m_dir = *val;
	};
	
	void setFile(std::string * val){
		*m_file = *val;
		
	};
	
	std::string * getDir(){
		
		return m_dir;
	};
	std::string * getFile(){
		return m_file;
	};
	
	void buildPorts(){
		bool built = true;
		
		//input
		built = (built && createPortIn<std::string*, Checker>(this, &Checker::setFile, PortType::M_FILENAMEPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
		built = (built && createPortIn<std::string*, Checker>(this, &Checker::setDir, PortType::M_DIRPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));	
		//output
		built = (built && createPortOut<std::string*,Checker>(this, &Checker::getFile, PortType::M_FILENAMEPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));
		built = (built && createPortOut<std::string*,Checker>(this, &Checker::getDir, PortType::M_DIRPTR, mimmo::pin::containerTAG::SCALAR, mimmo::pin::dataTAG::STRING_));	
		m_arePortsBuilt = built;
	};
	void execute(){};
};



void testBinStrings() {

	
 	Checker * reader = new Checker(0);
   	Checker * writer = new Checker(1);

	//settings;
	reader->setDir("./geo_data");
	reader->setFile("provaPoints2");
	
	writer->setDir("./franco");
	writer->setFile("carni");
	

	//Set PINS
	cout << "set pins" << endl;
 	
 	std::cout << " add pin " <<  boolalpha << addPin(reader, writer, PortType::M_FILENAMEPTR, PortType::M_FILENAMEPTR) << std::endl;
 	std::cout << " add pin " <<  boolalpha << addPin(reader, writer, PortType::M_DIRPTR, PortType::M_DIRPTR) << std::endl;

	cout << "set pins done" << endl;
	
	// 	//Create chain
	mimmo::Chain ch0;

 	ch0.addObject(reader);
  	//ch0.addObject(writer);
	

	//Execution of chain
// 	cout << "execution start" << endl;
// 	steady_clock::time_point t1 = steady_clock::now();
// 	ch0.exec();
// 	steady_clock::time_point t2 = steady_clock::now();
// 	cout << "execution done " << endl;
// 	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
// 	std::cout << "execution took me " << time_span.count() << " seconds."<<std::endl;	

	reader->exec();
	
	std::cout<<*(writer->m_dir)<<'\t'<<writer->m_dir->empty()<<'\t'<<&(writer->m_dir)<<'\t'<<std::endl;	
	std::cout<<*(writer->m_file)<<'\t'<<writer->m_file->empty()<<'\t'<<&(writer->m_file)<<'\t'<<std::endl;	
	
  	delete reader;
  	reader = NULL;

	std::cout<<&(writer->m_dir)<<'\t'<<writer->m_dir->empty()<<'\t'<<*(writer->m_dir)<<'\t'<<std::endl;
	std::cout<<&(writer->m_file)<<'\t'<<writer->m_file->empty()<<'\t'<<*(writer->m_file)<<'\t'<<std::endl;	
	
// 	std::cout<<writer->m_dir<<'\t'<<writer->m_dir.empty()<<'\t'<<&(writer->m_dir)<<'\t'<<std::endl;
// 	std::cout<<writer->m_file<<'\t'<<writer->m_file.empty()<<'\t'<<&(writer->m_file)<<'\t'<<std::endl;	

	
// 	delete reader;
// 	reader = NULL;
	delete writer;
	writer = NULL;
	
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

#if ENABLE_MPI==1
	MPI::Init(argc, argv);

	{
#endif
		/**<Calling mimmino Test routines*/

		testBinStrings();
		
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return 0;
}

