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
 \ *---------------------------------------------------------------------------*/

#include "mimmo_iogeneric.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Reading a generic input from file with GenericInput.
 */
int test2_1() {
	
    GenericInput * ginput = new GenericInput();
    ginput->setReadFromFile(true);
    ginput->setCSV(true);
    ginput->setReadDir("input");
    ginput->setFilename("generic_input_00001.csv");
    ginput->exec();
    
    dvecarr3E result = ginput->getResult<dvecarr3E>();
    
    bool check = result.size() == 7;

    std::cout<<"test passed :"<<check<<std::endl;
    
    delete ginput;
    return int(!check);
}

/*! 
 * Read a write a MimmoPiercedVector with classes GenericInputMPVData and GenericOutputMPVData.
 */
int test2_2() {
    
    //create a fake geometry();
    dvecarr3E points(4, {{0,0,0}});
    livector2D  conn(2, livector1D(3,0));
    points[1][0] = 1.0;
    points[2][1] = 1.0;
    points[3][0] = 1.0;
    points[3][1] = 1.0;
    conn[0][1] = 1; conn[0][2]=2;
    conn[1][0] = 1; conn[1][1]=3; conn[1][2] = 2;
    
    MimmoObject * geo = new MimmoObject(1);
    for(int i=0; i<(int)points.size(); ++i) geo->addVertex(points[i], i);
    geo->addConnectedCell(conn[0], bitpit::ElementInfo::TRIANGLE, short(0), long(0));
    geo->addConnectedCell(conn[1], bitpit::ElementInfo::TRIANGLE, short(0), long(1));
    
    //create a MPV vector of doubles and darray3E;
    MimmoPiercedVector<double> scalar(geo, mimmo::MPVLocation::CELL);
    MimmoPiercedVector<darray3E> vector(geo, mimmo::MPVLocation::POINT);
    scalar.insert(0, 12.12345);
    scalar.insert(1, -3.456);
    vector.insert(1,{{-1.0, 0, 2.0}});
    vector.insert(2,{{-0.976, -0.976, -0.976}});
    vector.insert(0,{{1.2, 1.3, 1.4}});
    vector.insert(3,{{0, 0, 12.0}});
    
    // write on file : scalar in csv and vector in raw ascii.
    GenericOutputMPVData * write_scalar = new GenericOutputMPVData();
    write_scalar->setWriteDir(".");
    write_scalar->setFilename("scalarCSV.csv");
    write_scalar->setCSV(true);
    write_scalar->setInput(scalar);

    GenericOutputMPVData * write_vector = new GenericOutputMPVData();
    write_vector->setWriteDir(".");
    write_vector->setFilename("vectorRAW.dat");
    write_vector->setCSV(false);
    write_vector->setBinary(false);
    write_vector->setInput(vector);
    
    write_scalar->exec();
    write_vector->exec();
    
    
    //re read files and absorb structures.
    GenericInputMPVData * read_scalar = new GenericInputMPVData();
    read_scalar->setReadDir(".");
    read_scalar->setFilename("scalarCSV.csv");
    read_scalar->setCSV(true);
    read_scalar->setGeometry(geo);
    
    GenericInputMPVData * read_vector = new GenericInputMPVData();
    read_vector->setReadDir(".");
    read_vector->setFilename("vectorRAW.dat");
    read_vector->setCSV(false);
    read_vector->setBinary(false);
    read_vector->setGeometry(geo);

    read_scalar->exec();
    read_vector->exec();
    
    
    //check re read structure are the same.
    auto rscalar = read_scalar->getResult<double>();
    auto rvector = read_vector->getResult<darray3E>();
    
    bool check = (rscalar.getGeometry() == scalar.getGeometry()); 
    check = check && (rvector.getGeometry() == vector.getGeometry()); 
    check = check && (rscalar.getDataLocation() == scalar.getDataLocation());
    check = check && (rvector.getDataLocation() == vector.getDataLocation());
    check = check && (rscalar.size() == scalar.size());
    check = check && (rvector.size() == vector.size());
    
    delete write_scalar;
    delete write_vector;
    delete read_scalar;
    delete read_vector;
    delete geo;

    std::cout<<"test passed :"<<check<<std::endl;

    return int(!check);
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

        int val = test2_1() ;
        std::max(val, test2_2());
#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
