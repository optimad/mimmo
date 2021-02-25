/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
 * Test 00001
 * Testing Executable blocks infrastructure: BaseManipulation - Pins - Chain execution
 */

class ManipA: public mimmo::BaseManipulation{
public:
	int m_member;
	dvecarr3E m_coords;
	dvector1D m_field;

	ManipA(){m_member = 1;};
	virtual ~ManipA(){};
	dvecarr3E getCoords(){ return m_coords;};
	dvector1D getField(){ return m_field; };
	void buildPorts(){

        //Register ports
		mimmo::PortManager::instance().addPort(M_COORDS, MC_VECARR3, MD_FLOAT,"test_core_00001.cpp");
		mimmo::PortManager::instance().addPort("M_MYPERSONALFIELD", MC_VECTOR, MD_FLOAT,"test_core_00001.cpp");
        bool built = true;
		built = built && createPortOut<dvecarr3E, ManipA>(this, &ManipA::getCoords, M_COORDS);
		built = built && createPortOut<dvector1D, ManipA>(this, &ManipA::getField, "M_MYPERSONALFIELD");
		m_arePortsBuilt = built;
	};
	void execute(){m_member = 4;};

};

class ManipB: public mimmo::BaseManipulation{
public:
	int m_member;
	dvecarr3E m_coords;
	dvector1D m_field;

	ManipB(){m_member = 1;};
	virtual ~ManipB(){};
	void setCoords(dvecarr3E points){ m_coords = points;};
	void setField(dvector1D field){m_field = field;};
	void buildPorts(){

        //Register ports
		mimmo::PortManager::instance().addPort(M_COORDS, MC_VECARR3, MD_FLOAT,"test_core_00001.cpp");
		mimmo::PortManager::instance().addPort("M_MYSCALARFIELD", MC_VECTOR, MD_FLOAT,"test_core_00001.cpp");

		bool built = true;
		built = built && createPortIn<dvecarr3E, ManipB>(this, &ManipB::setCoords, M_COORDS);
		built = built && createPortIn<dvector1D, ManipB>(this, &ManipB::setField, "M_MYSCALARFIELD");
		m_arePortsBuilt = built;
	};
	void execute(){m_member = 4;};


};

// =================================================================================== //

int test1() {
	//testing instantiation
	ManipA * objA = new ManipA();
	ManipB * objB = new ManipB();

	if(objA->m_member != 1 || objB->m_member != 1){
		std::cout<<"Failed instantiation"<<std::endl;
        delete objA;
        delete objB;
		return 1;
	}else{
		std::cout<<"Correct instantiation"<<std::endl;
	}
	//testing connections

	bool checkPin = true;
	dvecarr3E points(3,{{1.5,1.5,1.5}});
	dvector1D field(3,-1.2);

	objA->m_coords = points;
	objA->m_field = field;

	checkPin = checkPin && mimmo::pin::addPin(objA, objB, M_COORDS, M_COORDS);
    checkPin = checkPin && mimmo::pin::addPin(objA, objB, "M_MYPERSONALFIELD", "M_MYSCALARFIELD");

	if(!checkPin){
		std::cout<<"Failed getting connections"<<std::endl;
		delete objA;
		delete objB;
		return 1;
	}else{
		std::cout<<"Connections created"<<std::endl;
	}

	//testing chain & execution

	mimmo::Chain * c0 = new mimmo::Chain();
	c0->addObject(objB);
	c0->addObject(objA);

	c0->exec(false);

	int sF = field.size();
	int sC = points.size();

	bool checkExec = (objA->m_member == 4);
	checkExec = checkExec && (objA->m_member == 4);
	checkExec = checkExec && ((int)objB->m_field.size() == sF);
	checkExec = checkExec && ((int)objB->m_coords.size() == sC);

    {
        std::unique_ptr<mimmo::Chain> chain2 = c0->clone();
        checkExec = checkExec && (chain2->getNObjects() == c0->getNObjects());
        checkExec = checkExec && (c0->getID() != chain2->getNObjects());
    }

    if(!checkExec){
		std::cout<<"Failed execution"<<std::endl;
		delete objA;
		delete objB;
		delete c0;
		return 1;
	}else{
		std::cout<<"Successfull execution"<<std::endl;
	}

	delete objA;
	delete objB;
	delete c0;

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
