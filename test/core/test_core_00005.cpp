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
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;

/*
 * Test 00005
 * Testing MimmoPiercedVector Copy/Assignment/Swap/squeezeOutExcept
 */

// =================================================================================== //

int test5() {

	MimmoPiercedVector<double> mpv1;
	MimmoObject * obj = new MimmoObject();

	mpv1.insert(13, 1.234);
	mpv1.insert(18, -2.0);
    mpv1.setGeometry(obj);
	MimmoPiercedVector<double> mpv2(mpv1);
	MimmoPiercedVector<double> mpv3 = mpv1;

	bool check = mpv2.exists(13) && mpv2.exists(18);
	check = check && (mpv3.exists(13) && mpv3.exists(18));
    check = check && (mpv2.getGeometry()== obj);
    check = check && (mpv3.getGeometry()== obj);

	if(!check){
		std::cout<<"Copy or assignment of MimmoPiercedVector failed -kernel verification"<<std::endl;
        delete obj;
        return 1;
	}

	check = check && ((mpv2[13] == mpv1[13]) && (mpv2[18] == mpv1[18]));
	check = check && ((mpv3[13] == mpv1[13]) && (mpv3[18] == mpv1[18]));

	if(!check){
		std::cout<<"Copy or assignment of MimmoPiercedVector failed -storage verification"<<std::endl;
		delete obj;
        return 1;
	}else{
		std::cout<<"Copy and assignment of MimmoPiercedVector succeded"<<std::endl;
	}

    //check the squeezeOutExcept method with order inside it.
    MimmoPiercedVector<long> mpv4;
    mpv4.setGeometry(obj);

    mpv4.insert(0, 13);
    mpv4.insert(1, 16);
    mpv4.insert(2, 19);
    mpv4.insert(3, 8);
    mpv4.insert(5, 1);
    mpv4.insert(6, 0);

    mpv4.squeezeOutExcept(std::vector<long>({{2,5}}), true);

    check = (mpv4.size() == 2);
    if(!check){
        std::cout<<"Squeeze out of MimmoPiercedVector failed"<<std::endl;
        delete obj;
        return 1;
    }else{
        std::cout<<"Squeeze out of MimmoPiercedVector succeded"<<std::endl;
        for(auto it=mpv4.begin(); it!=mpv4.end(); ++it){
            std::cout<<"Retained id "<<it.getId()<<" with value "<<*it<<std::endl;
        }
    }

    delete obj;
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
        val = test5() ;
    }
    catch(std::exception & e){
        std::cout<<"test_core_00005 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
