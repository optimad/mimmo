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

#include "mimmo_utils.hpp"
#include <exception>
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Test: mirroring a point with data attached w.r.t a plane
 */
int test1() {

    dvecarr3E points(1,{{0,1,0}});
    dvecarr3E data(1, {{1,2,3}});

    SpecularPoints * sp = new SpecularPoints();
    sp->setCoords(points);
    sp->setVectorData(data);
    sp->setPlane({{0,0,0}},{{0,1,0}});
    sp->exec();

    auto coords = sp->getCloudResult();
    auto datas = sp->getCloudVectorData();

    bool check = ((coords.size() ==2 ) && (datas.size() ==2) );
    if(!check)  {
        delete sp;
        return 1;
    }
    check  = ( std::abs(coords[0][1] + coords[1][1])  < 1.E-18 );
    check  = check && ( std::abs(datas[1][0] - datas[0][0])  < 1.E-18 );
    check  = check && ( std::abs(datas[1][2] - datas[0][2])  < 1.E-18 );
    check  = check && ( std::abs(datas[1][1] + datas[0][1])  < 1.E-18 );
//     std::cout<<coords[0]<<'\t'<<datas[0]<<std::endl;
//     std::cout<<coords[1]<<'\t'<<datas[1]<<std::endl;
    std::cout<<"test passed : "<<check<<std::endl;

    delete sp;

    return int(!check);
}

// =================================================================================== //

int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
		/**<Calling mimmo Test routines*/
        int val = 1;
        try{
            val = test1();
        }

        catch(std::exception & e){
            std::cout<<"test_utils_00001 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
