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

#include "mimmo_geohandlers.hpp"
using namespace std;
using namespace bitpit;
using namespace mimmo;



// =================================================================================== //
/*!
 * Testing geohandlers module. Stitching geometries, Splitting and Overlapping fields defined on them 
 */
int test1() {
	
    //define 3 single triangle mesh
    MimmoObject * m1 = new MimmoObject(1);
    MimmoObject * m2 = new MimmoObject(1);
    MimmoObject * m3 = new MimmoObject(1);
    
    dvecarr3E points(5, {{0.0,0.0,0.0}});
    points[1] = {{1.0,0.0,0.0}};
    points[2] = {{2.0,0.0,0.0}};
    points[3] = {{0.5,1.0,0.0}};
    points[4] = {{1.5,1.0,0.0}};

    livector1D conn(3, 0);
    
    m1->addVertex(points[0],0);
    m1->addVertex(points[1],1);
    m1->addVertex(points[4],4);
    conn[0] = 0; conn[1] = 1; conn[2] = 4;
    m1->addConnectedCell(conn, bitpit::ElementInfo::Type::TRIANGLE, 0, 4);
    
    m2->addVertex(points[1],1);
    m2->addVertex(points[3],3);
    m2->addVertex(points[4],4);
    conn[0] = 1; conn[1] = 3; conn[2] = 4;
    m2->addConnectedCell(conn, bitpit::ElementInfo::Type::TRIANGLE, 1, 9);
    
    
    m3->addVertex(points[1],1);
    m3->addVertex(points[2],2);
    m3->addVertex(points[3],3);
    conn[0] = 1; conn[1] = 2; conn[2] = 3;
    m3->addConnectedCell(conn, bitpit::ElementInfo::Type::TRIANGLE, 2, 12);


    //stitch geometries
    StitchGeometry * stitch1 = new StitchGeometry(1);
    StitchGeometry * stitch2 = new StitchGeometry(1);
    
    stitch1->setAddGeometry(m1);
    stitch1->setAddGeometry(m2);
    
    stitch2->setAddGeometry(m2);
    stitch2->setAddGeometry(m3);
    
    stitch1->exec();
    stitch2->exec();
    
    if(stitch1->getGeometry()->getNCells() !=2 || stitch2->getGeometry()->getNCells() !=2) return 1;

    //define to fake fields on stitched geometries
    dvector1D field1(2, 1.0);
    dvector1D field2(2, 2.0);
    
    //Split fields
    SplitScalarField * split1 = new SplitScalarField();
    SplitScalarField * split2 = new SplitScalarField();
    
    split1->setGeometry(stitch1->getGeometry());
    split1->setField(field1);
    split1->setSplittedGeometries(stitch1->getOriginalGeometries());
    split1->setCellDivisionMap(stitch1->getCellDivisionMap());
    
    split2->setGeometry(stitch2->getGeometry());
    split2->setField(field2);
    split2->setSplittedGeometries(stitch2->getOriginalGeometries());
    split2->setCellDivisionMap(stitch2->getCellDivisionMap());
    
    split1->exec();
    split2->exec();
    
    //overlap Fields
    OverlapScalarFields * olap = new OverlapScalarFields();
    
    olap->setOverlapCriterium(4);
    olap->setDataFieldMap(split1->getSplittedData());
    olap->setDataFieldMap(split2->getSplittedData());
    
    olap->exec();
    
    auto datafield = olap->getDataFieldMap();
    
    bool check = true;
    for(auto &val: datafield ){
        
        if(val.first == m1) check = check && ( (*(val.second))[0] == 1.0);
        if(val.first == m2) check = check && ( (*(val.second))[0] == 3.0);
        if(val.first == m3) check = check && ( (*(val.second))[0] == 2.0);
        
    }
    std::cout<<"test passed :"<<check<<std::endl; 
    delete m1;
    delete m2;
    delete m3;
    delete stitch1;
    delete stitch2;
    delete split1;
    delete split2;
    delete olap;
    
    
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

        int val = test1() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
