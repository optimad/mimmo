/* test 002 -> Basic Elementary Shapes of libRA_BasicEleShapes
 * Testing all methods of all classes defining an Elementary shape.
 */

 using namespace std;
#include <iostream>
#include <fstream> 
#include "libRA_BasicEleShapes.hpp"

int main(int argc, char*argv[]){
	
std::ofstream out;
out.open("testingElementaryShapes.log");
if(out.is_open()){
		
	//testing BASE methods.
	{
		out<<"Testing BASE_ElementalShapes methods..."<<endl;
		darray3E origin_, origin2_;
		bool check;
		origin_.fill(-4); origin_[1]=2;
		origin2_.fill(1); 
		
		ElementalCube * cube1 = new ElementalCube;
		ElementalCube * cube2 = new ElementalCube(origin_,4,4,4);
		
		check = (cube1->getClassType()=="Cube");
		if(check){
			out<<"...test getClassType method...OK"<<endl;
		}else{
			out<<"...test getClassType method...FAILED"<<endl;		
		}
		
		check = (cube2->getShapeOrigin()==origin_);
		if(check){
			out<<"...test getShapeOrigin method...OK"<<endl;
		}else{
			out<<"...test getShapeOrigin method...FAILED"<<endl;		
		}
		
		cube1->setShapeOrigin(origin2_[0], origin2_[1],origin2_[2]);
		check = (cube1->getShapeOrigin()==origin2_);
		if(check){
			out<<"...test setShapeOrigin method...OK"<<endl;
		}else{
			out<<"...test setShapeOrigin method...FAILED"<<endl;		
		}
		
		origin2_.fill(0.0);
		cube1->clearShape();
		
		check = (cube1->getShapeOrigin()==origin2_);
		if(check){
			out<<"...test clearShape method...OK"<<endl;
		}else{
			out<<"...test clearShape method...FAILED"<<endl;		
		}
		
		
		*(static_cast<BASE_ElementalShapes* > (cube1)) = *(static_cast<BASE_ElementalShapes* > (cube2));
		check = (cube1->getShapeOrigin()==cube2->getShapeOrigin());
		if(check){
			out<<"...test copyOperator method...OK"<<endl;
		}else{
			out<<"...test copyOperator method...FAILED"<<endl;		
		}
		delete cube1; cube1=NULL;
		delete cube2; cube2=NULL;
	}	
	
	//testing ElementalCube;
	{
		out<<"Testing ElementalCube methods..."<<endl;
		dvector1D span(3,1);
		span[0]=2.0; span[2]=3.0;
		darray3E or_;
		or_.fill(4.0);
		
		ElementalCube * cube1 = new ElementalCube;
		ElementalCube * cube2 = new ElementalCube(or_, span[0],span[1], span[2]);
		bool check;	
		
		check = (cube2->getShapeSpan()==span);
		if(check){
			out<<"...test getShapeSpan method...OK"<<endl;
		}else{
			out<<"...test getShapeSpan method...FAILED"<<endl;		
		}
		cube1->setShapeSpan(span[0], span[1],span[2]);
		check = (cube1->getShapeSpan() == span);
		if(check){
			out<<"...test setShapeSpan method...OK"<<endl;
		}else{
			out<<"...test setShapeSpan method...FAILED"<<endl;		
		}
		
		dvector1D span2(3,0);
		cube1->clearShape();
		
		check = (cube1->getShapeSpan()==span2);
		if(check){
			out<<"...test clearShape method...OK"<<endl;
		}else{
			out<<"...test clearShape method...FAILED"<<endl;		
		}
		
		*cube1 = *cube2;
		check = (cube1->getShapeSpan()== cube2->getShapeSpan());
		check = check && (cube1->getShapeOrigin()== cube2->getShapeOrigin());
		check = check && (cube1->getClassType()== cube2->getClassType());
		if(check){
			out<<"...test copyOperator method for "<<cube1->getClassType()<<" ...OK"<<endl;
		}else{
			out<<"...test copyOperator method for "<<cube1->getClassType()<<" ...FAILED"<<endl;		
		}
		delete cube1; cube1=NULL;
		delete cube2; cube2=NULL;
		
	}
	
	//testing ElementalCylinder;
	{
		out<<"Testing ElementalCylinder methods..."<<endl;
		dvector1D span(2,2);
		span[0]=1.0;
		darray3E or_;
		or_.fill(2.0);
		
		ElementalCylinder * cyl1 = new ElementalCylinder;
		ElementalCylinder * cyl2 = new ElementalCylinder(or_, span[0],span[1]);
		bool check;	
		
		check = (cyl2->getShapeSpan()==span);
		if(check){
			out<<"...test getShapeSpan method...OK"<<endl;
		}else{
			out<<"...test getShapeSpan method...FAILED"<<endl;		
		}
		
		cyl1->setShapeSpan(span[0], span[1]);
		check = (cyl1->getShapeSpan()==span);
		if(check){
			out<<"...test setShapeSpan method...OK"<<endl;
		}else{
			out<<"...test setShapeSpan method...FAILED"<<endl;		
		}
		
		dvector1D span2(2,0);
		cyl1->clearShape();
		
		check = (cyl1->getShapeSpan()==span2);
		if(check){
			out<<"...test clearShape method...OK"<<endl;
		}else{
			out<<"...test clearShape method...FAILED"<<endl;		
		}
		
		
		*cyl1 = *cyl2;
		check = (cyl1->getShapeSpan()== cyl2->getShapeSpan());
		check = check && (cyl1->getShapeOrigin()== cyl2->getShapeOrigin());
		check = check && (cyl1->getClassType()== cyl2->getClassType());
		if(check){
			out<<"...test copyOperator method for "<<cyl1->getClassType()<<" ...OK"<<endl;
		}else{
			out<<"...test copyOperator method for "<<cyl1->getClassType()<<" ...FAILED"<<endl;		
		}
		delete cyl1; cyl1=NULL;
		delete cyl2; cyl2=NULL;
		
	}

	//testing ElementalSphere;
	{
		out<<"Testing ElementalSphere methods..."<<endl;
		double span = 3;
		darray3E or_;
		or_.fill(0.0);
		
		ElementalSphere * s1 = new ElementalSphere;
		ElementalSphere * s2 = new ElementalSphere(or_, span);
		bool check;	
		
		check = (s2->getShapeSpan()==span);
		if(check){
			out<<"...test getShapeSpan method...OK"<<endl;
		}else{
			out<<"...test getShapeSpan method...FAILED"<<endl;		
		}
		
		s1->setShapeSpan(span);
		check = (s1->getShapeSpan()==span);
		if(check){
			out<<"...test setShapeSpan method...OK"<<endl;
		}else{
			out<<"...test setShapeSpan method...FAILED"<<endl;		
		}
		
		double span2=0.0;
		s1->clearShape();
		
		check = (s1->getShapeSpan()==span2);
		if(check){
			out<<"...test clearShape method...OK"<<endl;
		}else{
			out<<"...test clearShape method...FAILED"<<endl;		
		}
		
		
		*s1 = *s2;
		check = (s1->getShapeSpan()== s2->getShapeSpan());
		check = check && (s1->getShapeOrigin()== s2->getShapeOrigin());
		check = check && (s1->getClassType()== s2->getClassType());
		if(check){
			out<<"...test copyOperator method for "<<s1->getClassType()<<" ...OK"<<endl;
		}else{
			out<<"...test copyOperator method for "<<s1->getClassType()<<" ...FAILED"<<endl;		
		}
		delete s1; s1=NULL;
		delete s2; s2=NULL;
		
	}
	
	
	//testing ElementalHexaHedron;
	{
		out<<"Testing ElementalHexahedron methods..."<<endl;
		
		dvecarr3E vertList(8,darray3E{0,0,0});
		vertList[1][0] = vertList[5][0]=1.0 ; 
		vertList[2][0] = vertList[2][1]= vertList[6][0] = vertList[6][1]=1.0;
		vertList[3][1] = vertList[7][1]= 1.0;
		vertList[4][2] = vertList[5][2]=vertList[6][2]=vertList[7][2] =1.0;
		
		ivector2D face_conn(6, ivector1D(4) );
		ivector2D vert_conn(8, ivector1D(3) );
		
		face_conn[0][0] = 0;  face_conn[0][1] = 3;   face_conn[0][2] = 7; face_conn[0][3] = 4;  //west side 
		face_conn[1][0] = 0;  face_conn[1][1] = 4;   face_conn[1][2] = 5; face_conn[1][3] = 1;  //south side 
		face_conn[2][0] = 0;  face_conn[2][1] = 1;   face_conn[2][2] = 2; face_conn[2][3] = 3;  //front side             
		face_conn[3][0] = 1;  face_conn[3][1] = 5;   face_conn[3][2] = 6; face_conn[3][3] = 2;  //east side             
		face_conn[4][0] = 3;  face_conn[4][1] = 2;   face_conn[4][2] = 6; face_conn[4][3] = 7; //north side 
		face_conn[5][0] = 4;  face_conn[5][1] = 7;   face_conn[5][2] = 6; face_conn[5][3] = 5; // back side 
	
		vert_conn[0][0] = 0; vert_conn[0][1] = 1; vert_conn[0][2] =2; //W-S-F plane intersection  
		vert_conn[1][0] = 1; vert_conn[1][1] = 2; vert_conn[1][2] =3; //S-F-E plane intersection  
		vert_conn[2][0] = 2; vert_conn[2][1] = 3; vert_conn[2][2] =4; //F-E-N plane intersection
		vert_conn[3][0] = 0; vert_conn[3][1] = 2; vert_conn[3][2] =4; //W-F-N plane intersection
		vert_conn[4][0] = 0; vert_conn[4][1] = 1; vert_conn[4][2] =5; //W-S-B plane intersection
		vert_conn[5][0] = 1; vert_conn[5][1] = 3; vert_conn[5][2] =5; //S-E-B plane intersection
		vert_conn[6][0] = 3; vert_conn[6][1] = 4; vert_conn[6][2] =5; //E-N-B plane intersection
		vert_conn[7][0] = 0; vert_conn[7][1] = 4; vert_conn[7][2] =5;
		
		dvecarr4E planes(6,darray4E{0,0,0,0});
		planes[0][0] =1.0; 
		planes[1][1] =1.0; 
		planes[2][2] =1.0;
		planes[3][0] =-1.0; planes[3][3] = 1.0;
		planes[4][1] =-1.0; planes[4][3] = 1.0;
		planes[5][2] =-1.0; planes[5][3] = 1.0;
		
		//////////////////////////////////////////////////////////////////
		
		ElementalHexaHedron * s1 = new ElementalHexaHedron(vertList);
		
		bool check=true;	
		dvecarr3E vert2  = s1->getHexaVertices();
		
		for(int i =0; i<8; i++){
			check = check && (vertList[i] == vert2[i]);
		}
		out<<"...test getHexaVertex/Vertices methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}
		
		check=true;	
		dvecarr4E pla  = s1->getHexaPlanes();
		for(int i =0; i<6; i++){
			check = check && (planes[i] == pla[i]);
		}
		out<<"...test getPlane/hexaPlanes methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}

		check=true;	
		ivector2D vconn  = s1->getVertConnectivity();
		for(int i =0; i<8; i++){
			check = check && (vert_conn[i] == vconn[i]);
		}
		out<<"...test getVertConn/Connectivity methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}
		
		check=true;	
		ivector2D fconn  = s1->getFaceConnectivity();
		for(int i =0; i<6; i++){
			check = check && (face_conn[i] == fconn[i]);
		}
		out<<"...test getFaceConn/Connectivity methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}

		
		darray3E P; P.fill(0.5); P[0]=0.0;
		check = s1->checkBelongPolygon(P,0);
		out<<"...test checkBelongPolygon methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}
		
		s1->clearShape();
		dvecarr3E vertDum = s1->getHexaVertices();
		check = (vertDum.size() == 0);
		out<<"...test clearShape/resetHexaHedron methods...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}
		
		s1->setBaseVertex(&vertList);
		dvecarr3E vert3 = s1->getHexaVertices();
		check = true;
		for(int i =0; i<8; i++){
			check = check && (vert3[i] == vert2[i]);
		}
		out<<"...test setBaseVertex method...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}

		
		P[0] = -0.5; P[1] = 1;P[2] = 0;

		s1->changeVertex(3, P, false);
		
		check = (s1->getHexaVertex(3) == P);
		out<<"...test changeVertex method...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}

		s1->setShapeOrigin(-0.5,0,0);
		darray3E prova1; prova1.fill(0.0);
		prova1[0] = -0.5;
		
		check = (s1->getShapeOrigin() == prova1);
		out<<"...test setShapeOrigin method...";
		if(check){
			out<<"OK"<<endl;
		}else{
			out<<"FAILED"<<endl;		
		}
		
		
		ElementalHexaHedron * s2 = new ElementalHexaHedron(vertList);
		*s1 = *s2;
		
		darray3E dum1 = s1->getHexaVertex(6);
		darray3E dum2 = s2->getHexaVertex(6);
		
		check = (dum1 == dum2);
		check = check && (s1->getShapeOrigin()== s2->getShapeOrigin());
		check = check && (s1->getClassType()== s2->getClassType());
		if(check){
			out<<"...test copyOperator method for "<<s1->getClassType()<<" ...OK"<<endl;
		}else{
			out<<"...test copyOperator method for "<<s1->getClassType()<<" ...FAILED"<<endl;		
		}
		delete s1; s1=NULL;
		delete s2; s2=NULL;
		
		exit(1);
	}
	
	
}else{
	cout<<"Not able to open log file"<<endl;
}

out.close();	
return(0);

} //end program 
