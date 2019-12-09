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

#include "mimmo_core.hpp"

/*
 * Test 00000
 * Testing Binary streams of a pool of 15 more or less complex structures.
 */

// =================================================================================== //

int test1() {

    mimmo::OBinaryStream outbuf;

    //1. testing a vector of simple type double.
    {
        std::vector<double> input(5, 1.12);

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<double>output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<5; ++i){
            check = check && input[i] == output[i];
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<double>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<double> buffering succeeded"<<std::endl;
    	}
    }

    //2. testing a vector of double pointers.
    {
        std::vector<double> input(5, 1.12);
        std::vector<double*> input_ptr(5, nullptr);
        for(int i=0; i<5; ++i){
            input_ptr[i] = &input[i];
        }

        outbuf.seekg(0); //clean it
        outbuf << input_ptr;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<double*>output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<5; ++i){
            check = check && input_ptr[i] == output[i];
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<double*>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<double*> buffering succeeded"<<std::endl;
    	}
    }

    //3. testing a vector of arrays
    {
        std::vector<std::array<long,4>> input(5, {{1,2,3,4}});

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<std::array<long,4> >output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<5; ++i){
            for (int j = 0; j<4; ++j){
                check = check && input[i][j] == output[i][j];
            }
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<array<.,.>>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<array<.,.>> buffering succeeded"<<std::endl;
    	}
    }

    //4. testing a vector of MimmoPiercedVector pointers
    {
        std::vector<mimmo::MimmoPiercedVector<std::string> * > input;
        mimmo::MimmoPiercedVector<std::string> m1;
        mimmo::MimmoPiercedVector<std::string> m2;

        input.push_back(&m1);
        input.push_back(nullptr);
        input.push_back(nullptr);
        input.push_back(&m2);


        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<mimmo::MimmoPiercedVector<std::string> *>output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<4; ++i){
            check = check && input[i] == output[i];
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<MimmoPiercedVector * >"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<MimmoPiercedVector * > buffering succeeded"<<std::endl;
    	}
    }

    //5. testing a vector of pairs
    {
        std::vector<std::pair<int*,std::string>> input(4);
        std::vector<int> integers(4, 126);

        input[0] = std::make_pair(&integers[0], "Led");
        input[1] = std::make_pair(&integers[1], "Zeppelin");
        input[2] = std::make_pair(&integers[2], "Still");
        input[3] = std::make_pair(&integers[3], "Rocks");

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<std::pair<int*,std::string>>output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<4; ++i){
            check = check && input[i].first == output[i].first;
            check = check && input[i].second == output[i].second;
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<pair<.,.>>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<pair<.,.>> buffering succeeded"<<std::endl;
    	}
    }

    //6. testing a vector of vector
    {
        std::vector<std::vector<double>> input(4,{{6.1,5.0,4.0,3.3,2.912}});

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::vector<std::vector<double> >output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        for(int i=0; i<4; ++i){
            for(int j=0; j<5; ++j){
                check = check && input[i][j] == output[i][j];
            }
        }
        if(!check){
    		std::cout<<"FAILED to buffer vector<vector<.>>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"vector<vector<.>> buffering succeeded"<<std::endl;
    	}
    }

    //7. testing a pair of pointers
    {
        double a = 2.3456;
        long b = 124567;
        std::pair<double*, long*> input = std::make_pair(&a,&b);

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::pair<double*, long* >output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        check = check && input.first == output.first;
        check = check && input.second == output.second;

        if(!check){
    		std::cout<<"FAILED to buffer pair<double*, long*>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"pair<double*, long*> buffering succeeded"<<std::endl;
    	}
    }

    //8. testing a pair of pointer, data
    {
        double a = 2.3456;
        long b = 124567;
        std::pair<long, double*> input = std::make_pair(b,&a);

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::pair<long, double* >output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= true;
        check = check && input.first == output.first;
        check = check && input.second == output.second;

        if(!check){
    		std::cout<<"FAILED to buffer pair<long, double*>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"pair<long, double*> buffering succeeded"<<std::endl;
    	}
    }

    //9. testing an unordered_map of pointer-key, pointer-argument
    {
        double a(2.3456), b(-3.0), c(12.12);
        long d(-1),e(-2345678),f(123123);
        std::unordered_map<double*, long*> input;
        input[&a] = &d;
        input[&b] = &e;
        input[&c] = &f;

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::unordered_map<double*, long*> output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= input.size() == output.size();
        std::unordered_map<double*, long*>::iterator inIt;
        inIt =  input.begin();
        while(check && inIt != input.end()){
            check = check && output.count(inIt->first)>0;
            check = check && inIt->second == output[inIt->first];
            ++inIt;
        }
        if(!check){
    		std::cout<<"FAILED to buffer unordered_map<double*, long*>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"unordered_map<double*, long*> buffering succeeded"<<std::endl;
    	}
    }

    //10. testing an unordered_map of data-key, data-argumet
    {
        double a(2.3456), b(-3.0), c(12.12);
        long d(-1),e(-2345678),f(123123);
        std::unordered_map<double, long> input;
        input[a] = d;
        input[b] = e;
        input[c] = f;

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::unordered_map<double, long> output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= input.size() == output.size();
        std::unordered_map<double, long>::iterator inIt;
        inIt =  input.begin();
        while(check && inIt != input.end()){
            check = check && output.count(inIt->first)>0;
            check = check && inIt->second == output[inIt->first];
            ++inIt;
        }
        if(!check){
    		std::cout<<"FAILED to buffer unordered_map<double, long>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"unordered_map<double, long> buffering succeeded"<<std::endl;
    	}
    }

    //11. testing an unordered_map of string-key, pair-argumet
    {
        std::unordered_map<std::string, std::pair<long, int>> input;
        input["bron"] = std::make_pair(1,2);
        input["y-aur"] = std::make_pair(3,4);
        input["stomp"] = std::make_pair(5,6);

        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        std::unordered_map<std::string, std::pair<long, int>> output;
        inbuf>>output;

        //check if vectors are equal;
        bool check= input.size() == output.size();
        std::unordered_map<std::string, std::pair<long, int>>::iterator inIt;
        inIt =  input.begin();
        while(check && inIt != input.end()){
            check = check && output.count(inIt->first)>0;
            check = check && inIt->second.first == output[inIt->first].first;
            check = check && inIt->second.second == output[inIt->first].second;
            ++inIt;
        }
        if(!check){
    		std::cout<<"FAILED to buffer unordered_map<string, pair<long,int> >"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"unordered_map<string, pair<long,int> > buffering succeeded"<<std::endl;
    	}
    }

    //12. testing an MimmoPiercedVector of strings
    {
        mimmo::MimmoObject geo;
        std::string name = "stringfield";
        mimmo::MimmoPiercedVector<std::string> input(&geo, mimmo::MPVLocation::POINT);
        input.setName(name);
        input.insert(0, "just");
        input.insert(1, "something");
        input.insert(19, "tofill");
        input.insert(22, "this");
        input.insert(99, "field");


        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        mimmo::MimmoPiercedVector<std::string> output;
        inbuf>>output;

        //check if they are equal;
        bool check= input.size() == output.size();
        check = check && input.getGeometry() == output.getGeometry();
        check = check && input.getName() == output.getName();
        check = check && static_cast<int>(input.getDataLocation()) == static_cast<int>(input.getDataLocation());
        mimmo::MimmoPiercedVector<std::string>::iterator inIt =  input.begin();
        while(check && inIt != input.end()){
            check = check && output.exists(inIt.getId());
            check = check && *inIt== output[inIt.getId()];
            ++inIt;
        }
        if(!check){
    		std::cout<<"FAILED to buffer MimmoPiercedVector<string>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"MimmoPiercedVector<string> buffering succeeded"<<std::endl;
    	}
    }

    //13. testing an MimmoPiercedVector of array<double,3>
    {
        mimmo::MimmoObject geo;
        std::string name = "vectorfield";
        mimmo::MimmoPiercedVector<std::array<double,3> > input(&geo, mimmo::MPVLocation::CELL);
        input.setName(name);
        input.insert(0, {{1.,2.,3.}});
        input.insert(1, {{4.,5.,6.}});
        input.insert(19,{{-7.,8.,-9.}});
        input.insert(22, {{1111.1,2.21,-3.3456}});
        input.insert(99, {{-8.765,0.,0.}});


        outbuf.seekg(0); //clean it
        outbuf << input;
        mimmo::IBinaryStream inbuf(outbuf.data(), outbuf.getSize());
        mimmo::MimmoPiercedVector<std::array<double,3> > output;
        inbuf>>output;

        //check if they are equal;
        bool check= input.size() == output.size();
        check = check && input.getGeometry() == output.getGeometry();
        check = check && input.getName() == output.getName();
        check = check && static_cast<int>(input.getDataLocation()) == static_cast<int>(input.getDataLocation());
        mimmo::MimmoPiercedVector<std::array<double,3> >::iterator inIt =  input.begin();
        while(check && inIt != input.end()){
            check = check && output.exists(inIt.getId());
            for(int i=0; i<3; ++i){
                check = check && inIt->at(i)== output[inIt.getId()][i];
            }
            ++inIt;
        }
        if(!check){
    		std::cout<<"FAILED to buffer MimmoPiercedVector<array<double,3>>"<<std::endl;
    		return 1;
    	}else{
    		std::cout<<"MimmoPiercedVector<array<double,3>> buffering succeeded"<<std::endl;
    	}
    }


    return 0;
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif

    int val = 1;
	/**<Calling mimmo Test routines*/
    try{
        val = test1() ;
    }
    catch(std::exception & e){
        std::cout<<"test_core_00001 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
