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
using namespace std;
using namespace bitpit;
using namespace mimmo;

/*
 * Test 00005
 * Testing MimmoPiercedVector Copy/Assignment/Swap
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
	std::cout<<mpv1<<std::endl;
	delete obj;
    return 0;
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

        int val = test5() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
