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
 * Test 00001
 * Testing Executable blocks infrastructure: BaseManipulation - Pins - Chain execution 
 */



class ManipA: public BaseManipulation{
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
        PortManager::instance().addPort(M_COORDS, MC_VECARR3, MD_FLOAT);
        PortManager::instance().addPort("M_MYPERSONALFIELD", MC_VECTOR, MD_FLOAT);
        bool built = true;
		built = built && createPortOut<dvecarr3E, ManipA>(this, &ManipA::getCoords, M_COORDS);
		built = built && createPortOut<dvector1D, ManipA>(this, &ManipA::getField, "M_MYPERSONALFIELD");
		m_arePortsBuilt = built;
	};
	void execute(){m_member = 4;};
	
};

class ManipB: public BaseManipulation{
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
        PortManager::instance().addPort(M_COORDS, MC_VECARR3, MD_FLOAT);
        PortManager::instance().addPort(M_SCALARFIELD, MC_VECTOR, MD_FLOAT);
        
		bool built = true;
		built = built && createPortIn<dvecarr3E, ManipB>(this, &ManipB::setCoords, M_COORDS);
		built = built && createPortIn<dvector1D, ManipB>(this, &ManipB::setField, M_SCALARFIELD);
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
	
	checkPin = checkPin && addPin(objA, objB, M_COORDS, M_COORDS);
    checkPin = checkPin && addPin(objA, objB, "M_MYPERSONALFIELD", M_SCALARFIELD);
	
	if(!checkPin){ 
		std::cout<<"Failed getting connections"<<std::endl;
		delete objA;
		delete objB;
		return 1;
	}else{
		std::cout<<"Connections created"<<std::endl;
	}	
	
	//testing chain & execution
	
	Chain * c0 = new Chain();
	c0->addObject(objB);
	c0->addObject(objA);
	
	c0->exec(false);
	
	int sF = field.size();
	int sC = points.size();
	
	bool checkExec = (objA->m_member == 4);
	checkExec = checkExec && (objA->m_member == 4);
	checkExec = checkExec && ((int)objB->m_field.size() == sF);
	checkExec = checkExec && ((int)objB->m_coords.size() == sC);
	
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
