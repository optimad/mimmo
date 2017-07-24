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



/*!
 * Creating surface triangular mesh and return it in a MimmoObject.
 * \param[in,out] mesh pointer to a MimmoObject mesh to fill.
 * \return true if successfully created mesh
 */
bool createMimmoMesh(MimmoObject * mesh){
    
    double dx = 0.25, dy = 0.25;
    int nV, nC;
    //create vertexlist
    dvecarr3E vertex(35,{{0.0,0.0,0.0}});
    livector2D conn(48, livector1D(3));
    
    for(int i=0; i<7; ++i){
        for(int j=0; j<5; j++){
            nV = 5*i + j;
            vertex[nV][0] = i*dx;
            vertex[nV][1] = j*dy;
        }
    }
    
    for(int j=0; j<4; ++j){
        for(int i=0; i<3; ++i){
            nC = 8*i + 2*j;
            
            conn[nC][0] = 5*i + j; 
            conn[nC][1] = 5*(i+1) + j;
            conn[nC][2] = 5*i + j+1;
            
            conn[nC+1][0] = 5*(i+1) + j; 
            conn[nC+1][1] = 5*(i+1) + j+1;
            conn[nC+1][2] = 5*i + j+1;
        }
    }
    
    for(int j=0; j<4; ++j){
        for(int i=3; i<6; ++i){
            nC = 8*i + 2*j;
            
            conn[nC][0] = 5*i + j; 
            conn[nC][1] = 5*(i+1) + j;
            conn[nC][2] = 5*(i+1) + j+1;
            
            conn[nC+1][0] = 5*i + j;  
            conn[nC+1][1] = 5*(i+1) + j+1;
            conn[nC+1][2] = 5*i + j+1;
        }
    }
    
    mesh->getVertices().reserve(35);
    mesh->getCells().reserve(48);
    
    //fill the mimmoObject;
    long cV=0;
    for(auto & val: vertex){
        mesh->addVertex(val, cV);
        cV++;
    }
    
    long cC=0;
    bitpit::ElementInfo::Type eltype = bitpit::ElementInfo::TRIANGLE;
    for(auto & val: conn){
        mesh->addConnectedCell(val, eltype, cC);
        cC++;
    }
    
    bool check = (mesh->getNCells() == 48) && (mesh->getNVertex() == 35);
    
    mesh->buildAdjacencies();
    return check;
}
// =================================================================================== //
/*!
 * Testing geohandlers module. Clipping, Selecting and Reconstructing Fields 
 */
int test2() {
	
    //define 3 single triangle mesh
    MimmoObject * m1 = new MimmoObject(1);
    if(!createMimmoMesh(m1)){
        delete m1;
        return 1;
    }    

    ClipGeometry * clip = new ClipGeometry();
    clip->setClipPlane({{0.75,0.5,0.0}}, {{1.0,0.0,0.0}});
    clip->setGeometry(m1);
    clip->setPlotInExecution(true);
    clip->exec();
    
    ClipGeometry * clip2 = new ClipGeometry();
    clip2->setInsideOut(true);
    clip2->setGeometry(m1);
    clip2->setClipPlane({{0.75,0.5,0.0}}, {{1.0,0.0,0.0}});
    clip2->setPlotInExecution(false);
    clip2->exec();
    
    
    SelectionBySphere * sel1 = new SelectionBySphere();
    sel1->setOrigin({{1.5,0.5,0.0}});
    sel1->setSpan({{0.6,2.0*M_PI, M_PI}});
    sel1->setDual(false);
    sel1->setGeometry(clip->getClippedPatch());
    sel1->setPlotInExecution(false);
    sel1->exec();


    SelectionBySphere * sel2 = new SelectionBySphere((*sel1));
    sel2->setDual(true);
    sel2->setPlotInExecution(false);
    sel2->exec();
    
 
    dmpvector1D field1(clip2->getGeometry());
    dmpvector1D field2(sel1->getPatch());
    dmpvector1D field3(sel2->getPatch());
    
    for (auto vertex : clip2->getGeometry()->getVertices()){
        long int ID = vertex.getId();
        field1.data().insert(ID, 1.0);
    }

    for (auto vertex : sel1->getGeometry()->getVertices()){
        long int ID = vertex.getId();
        field2.data().insert(ID, 1.0);
    }

    for (auto vertex : sel2->getGeometry()->getVertices()){
        long int ID = vertex.getId();
        field3.data().insert(ID, 1.0);
    }

	
    ReconstructScalar * recon = new ReconstructScalar();
    
    recon->setGeometry(m1);
    recon->setOverlapCriterium(3);
    recon->addData(field1);
    recon->addData(field2);
    recon->addData(field3);
    
    recon->exec();
    

    auto finalfield = recon->getResultField();
    
    bool check= true;
    for(auto & val : finalfield.data()){
        check = check && (val == 1.0);
    }


    SwitchScalarField * swtch = new SwitchScalarField();

    swtch->setFields(recon->getResultFields());
    swtch->setGeometry(sel1->getPatch());
    swtch->setPlotInExecution(true);
    swtch->exec();

    
    for(auto & val : swtch->getSwitchedField().getDataAsVector()){
       	check = check && (val == 1.0);
    }

    std::cout<<"plotting check"<<std::endl;	
    delete m1;
    delete clip;
    delete clip2;
    delete sel1;
    delete sel2;
    delete recon;
    delete swtch;
    
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

        int val = test2() ;

#if ENABLE_MPI==1
	}

	MPI::Finalize();
#endif
	
	return val;
}
