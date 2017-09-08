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
 * Testing geohandlers module. Copy and Assignment of ExtractXXX abstract and derived class.
 * Copy and assignment of Selection classes.
 */
int test3_1() {

    ExtractScalarField *sca = new ExtractScalarField();
    MimmoObject * obj = new MimmoObject();
    dmpvector1D ciccio;
    ciccio.insert(23,-12.456767);
    
    sca->setGeometry(obj);
    sca->setId(34);
    sca->setMode(ExtractMode::MAPPING);
    sca->setTolerance(1.2345E-4);
    sca->setField(ciccio);
        
    ExtractScalarField *scaCC = new ExtractScalarField(*sca);
    ExtractScalarField *scaAO = new ExtractScalarField();
    *scaAO = *sca;
    
    bool check= true;
    //verify copy constructor result content
    check = check && (scaCC->getGeometry() == obj);
    check = check && (scaCC->getId() != sca->getId());
    check = check && (scaCC->getMode() == sca->getMode());
    check = check && (scaCC->getTolerance() == sca->getTolerance());
    auto fCC = scaCC->getOriginalField();
    check = check  && (fCC.exists(23));

    if(!check){
        delete sca;
        delete scaCC;
        delete scaAO;
        delete obj;
        std::cout<<"Failing copy construction of geohandlers ExtractFields"<<std::endl;
        return 1;
    }

    //verify copy constructor result content
    check = check && (scaAO->getGeometry() == obj);
    check = check && (scaAO->getId() != sca->getId());
    check = check && (scaAO->getMode() == sca->getMode());
    check = check && (scaAO->getTolerance() == sca->getTolerance());
    auto fAO = scaAO->getOriginalField();
    check = check  && (fAO.exists(23));
    
    if(!check){
        delete sca;
        delete scaCC;
        delete scaAO;
        delete obj;
        std::cout<<"Failing assignment of geohandlers ExtractFields"<<std::endl;
        return 1;
    }
    std::cout<<"test 1 passed :"<<check<<std::endl; 
    
    delete sca;
    delete scaCC;
    delete scaAO;
    delete obj;
    
    return 0;
}

int test3_2() {
    
    SelectionByBoxWithScalar * sca = new SelectionByBoxWithScalar();
    MimmoObject * obj = new MimmoObject();
    dmpvector1D ciccio;
    ciccio.insert(23,-12.456767);
    
    sca->setGeometry(obj);
    sca->setOrigin({{1,2,3}});
    sca->setDual(true);
    sca->setField(ciccio);
    
    SelectionByBoxWithScalar *scaCC = new SelectionByBoxWithScalar(*sca);
    SelectionByBoxWithScalar *scaAO = new SelectionByBoxWithScalar();
    *scaAO = *sca;
    
    bool check= true;
    
    auto fCC = scaCC->getField();
    
    //verify copy constructor result content
    check = check && (scaCC->getGeometry() == obj);
    check = check && (scaCC->getOrigin() == sca->getOrigin());
    check = check && (scaCC->isDual());
    
    check = check  && (fCC.exists(23));
    
    if(!check){
        delete sca;
        delete scaCC;
        delete scaAO;
        delete obj;
        std::cout<<"Failing copy construction of geohandlers SelectionByBoxWithScalar"<<std::endl;
        return 1;
    }
    

    auto fAO = scaAO->getField();
    
//     std::cout<<scaAO->getGeometry()<<'\t'<<obj<<std::endl;
//     std::cout<<scaAO->getOrigin()<<'\t'<<sca->getOrigin()<<std::endl;
//     std::cout<<scaAO->isDual()<<std::endl;
//     std::cout<<fAO.exists(23)<<std::endl;
    
    //verify copy constructor result content
    check = check && (scaAO->getGeometry() == obj);
    check = check && (scaAO->getOrigin() == sca->getOrigin());
    check = check && (scaAO->isDual());
    
    check = check  && (fAO.exists(23));
    
    if(!check){
        delete sca;
        delete scaCC;
        delete scaAO;
        delete obj;
        std::cout<<"Failing assignment of geohandlers SelectionByBoxWithScalar"<<std::endl;
        return 1;
    }
    std::cout<<"test 2 passed :"<<check<<std::endl; 
    
    delete sca;
    delete scaCC;
    delete scaAO;
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

        int val = test3_1();
        val = std::max(val, test3_2());

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
