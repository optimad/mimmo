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
#include <exception>
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
   dvecarr3E pp(3,{{0,0,0}});
   pp[0][0] =1;
   pp[1][1] =1;
   pp[2][2] =1;
   livector1D conn(3,0);
   conn[1] = 1;
   conn[2] = 2;

   obj->addVertex(pp[0],0);
   obj->addVertex(pp[1],1);
   obj->addVertex(pp[2],2);
   obj->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, 0, 0);

   dmpvector1D ciccio;
   ciccio.insert(23,-12.456767);

   sca->setGeometry(obj);
   sca->setId(34);
   sca->setMode(ExtractMode::MAPPING);
   sca->setTolerance(1.2345E-4);
   sca->setField(&ciccio);

   ExtractScalarField *scaCC = new ExtractScalarField(*sca);
   ExtractScalarField *scaAO = new ExtractScalarField();
    *scaAO = *sca;

   MimmoPiercedVector<double> fCC = scaCC->getOriginalField();

   // std::cout<<scaCC->getGeometry()<<'\t'<<obj<<std::endl;
   // std::cout<<scaCC->getId()<<'\t'<<sca->getId()<<std::endl;
   // std::cout<<(int)scaCC->getMode()<<'\t'<<(int)sca->getMode()<<std::endl;
   // std::cout<<scaCC->getTolerance()<<'\t'<<sca->getTolerance()<<std::endl;
   // std::cout<<fCC.exists(23)<<std::endl;

   bool check= true;
   //verify copy constructor result content
   // std::cout << "check scaCC" << std::endl;
   check = check && (scaCC->getGeometry() == obj);
   check = check && (scaCC->getId() != sca->getId());
   check = check && (scaCC->getMode() == sca->getMode());
   check = check && (scaCC->getTolerance() == sca->getTolerance());
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
 //  std::cout << "check scaAO" << std::endl;
   check = check && (scaAO->getGeometry() == obj);
   // std::cout << "check get Id" << std::endl;
   check = check && (scaAO->getId() != sca->getId());
   // std::cout << "check get Mode" << std::endl;
   check = check && (scaAO->getMode() == sca->getMode());
   // std::cout << "check get tol" << std::endl;
   check = check && (scaAO->getTolerance() == sca->getTolerance());
   // std::cout << "check get field" << std::endl;
   MimmoPiercedVector<double> fAO;
   fAO = scaAO->getOriginalField();
   // std::cout << "check exists" << std::endl;
   check = check  && (fAO.exists(23));

   // std::cout << "if !check" << std::endl;
   if(!check){
       // std::cout << "delete" << std::endl;
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
   dvecarr3E pp(3,{{0,0,0}});
   pp[0][0] =1;
   pp[1][1] =1;
   pp[2][2] =1;
   livector1D conn(3,0);
   conn[1] = 1;
   conn[2] = 2;

   obj->addVertex(pp[0],0);
   obj->addVertex(pp[1],1);
   obj->addVertex(pp[2],2);
   obj->addConnectedCell(conn, bitpit::ElementType::TRIANGLE, 0, 0);

   dmpvector1D ciccio;
   ciccio.insert(23,-12.456767);

   sca->setGeometry(obj);
   sca->setOrigin({{1,2,3}});
   sca->setDual(true);
   sca->setField(&ciccio);

   SelectionByBoxWithScalar *scaCC = new SelectionByBoxWithScalar(*sca);
   SelectionByBoxWithScalar *scaAO = new SelectionByBoxWithScalar();
   *scaAO = *sca;

   bool check= true;

   MimmoPiercedVector<double>* fCC = scaCC->getField();

       // std::cout<<scaCC->getGeometry()<<'\t'<<obj<<std::endl;
       // std::cout<<scaCC->getOrigin()<<'\t'<<sca->getOrigin()<<std::endl;
       // std::cout<<scaCC->isDual()<<std::endl;
       // std::cout<<fCC.exists(23)<<std::endl;

   //verify copy constructor result content
   check = check && (scaCC->getGeometry() == obj);
   check = check && (scaCC->getOrigin() == sca->getOrigin());
   check = check && (scaCC->isDual());

   check = check  && (fCC->exists(23));

   if(!check){
       delete sca;
       delete scaCC;
       delete scaAO;
       delete obj;
       std::cout<<"Failing copy construction of geohandlers SelectionByBoxWithScalar"<<std::endl;
       return 1;
   }


   MimmoPiercedVector<double>*fAO = scaAO->getField();

//     std::cout<<scaAO->getGeometry()<<'\t'<<obj<<std::endl;
//     std::cout<<scaAO->getOrigin()<<'\t'<<sca->getOrigin()<<std::endl;
//     std::cout<<scaAO->isDual()<<std::endl;
//     std::cout<<fAO.exists(23)<<std::endl;

   //verify copy constructor result content
   check = check && (scaAO->getGeometry() == obj);
   check = check && (scaAO->getOrigin() == sca->getOrigin());
   check = check && (scaAO->isDual());

   check = check  && (fAO->exists(23));

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

#if MIMMO_ENABLE_MPI
	MPI_Init(&argc, &argv);
#endif
        int val = 1;

		/**<Calling mimmo Test routines*/
        try{
            val = test3_1();
            val = std::max(val, test3_2());
        }
        catch(std::exception & e){
            std::cout<<"test_geohandlers_00003 exited with an error of type : "<<e.what()<<std::endl;
            return 1;
        }
#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
