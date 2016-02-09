/* test 000 -> custom Operators testing
 * Testing all general low level functions, stored in customOperators
 */
 using namespace std;
#include <iostream>
#include <fstream> 
#include "customOperators.hpp"
 
int main(int argc, char*argv[]){

/* General test all utilities stored in customOperators */

std::ofstream out;
out.open("testingCustomOperators.log");
if(out.is_open()){
	
	//freeContainer
	{
		std::vector< dvector2D > test(10, dvector2D(12, dvector1D(4,-1)));
		freeContainer(test);
		if(test.size() != 0){
			out<<"Test for freeContainer function...FAILED"<<std::endl;
		}else{
			out<<"Test for freeContainer function...OK"<<std::endl;
		}
	}
	
	//summing vector
	{
		ivector1D test(10, 1);
		bvector1D test2(10, true);
		for(int i=1; i<10; ++i){test[i] = test[i-1]+1;}
		int result = summing_Vector(test);
		bool result2 = false;
		result2 = summing_Vector(test2);
		if(result != 55 || !result2){
			std::cout<<"Test for summing_Vector function...FAILED"<<std::endl;
		}else{
			out<<"Test for summing_Vector function...OK"<<std::endl;
		}
	}
	
	// logical_summing
	{
		{
			bool check = true;
			//boolean overloading;
			bool r1,t1;
			r1= t1 = true;
			bool r2,t2;
			r2= t2 = false;
			bool * policy = new bool[3];
			
			//first policy -> returning target
			policy[0] = true; policy[1]=false; policy[2]=false;
			check = check && (t1 == logical_summing(r1,t1,policy));
			check = check && (t2 == logical_summing(r1,t2,policy));
			
			//second policy -> returning sum
			policy[0] = false; policy[1]=true; policy[2]=false;
			check = check && (true == logical_summing(r1,t1,policy));
			check = check && (true == logical_summing(r2,t1,policy));
			check = check && (true == logical_summing(r1,t2,policy));
			check = check && (false == logical_summing(r2,t2,policy));
			
			//third policy -> returning difference
			policy[0] = false; policy[1]=false; policy[2]=true;
			check = check && (false == logical_summing(r1,t1,policy));
			check = check && (false == logical_summing(r2,t1,policy));
			check = check && (true == logical_summing(r1,t2,policy));
			check = check && (false == logical_summing(r2,t2,policy));
			
			delete[] policy; policy=NULL;
			
			if(!check){
				out<<"Test for boolean logical_summing function...FAILED"<<std::endl;
			}else{
				out<<"Test for boolean logical_summing function...OK"<<std::endl;
			}
		}
		
		{
			bool check = true;
			//boolean overloading;
			int r1,t1;
			r1 = t1 = 1;
			int r2,t2;
			r2 = t2 = 0;
			
			//first policy -> returning sum
			check = check && (1 == logical_summing(r1,t1,true));
			check = check && (1 == logical_summing(r2,t1,true));
			check = check && (1 == logical_summing(r1,t2,true));
			check = check && (0 == logical_summing(r2,t2,true));
			
			//second policy -> returning difference
			check = check && (0 == logical_summing(r1,t1,false));
			check = check && (0 == logical_summing(r2,t1,false));
			check = check && (1 == logical_summing(r1,t2,false));
			check = check && (0 == logical_summing(r2,t2,false));
			
			if(!check){
				out<<"Test for int logical_summing function...FAILED"<<std::endl;
			}else{
				out<<"Test for int logical_summing function...OK"<<std::endl;
			}
		}
		
		
	}
	
	//getSign
	{
		double test = -1.23456;
		bool check = (-1.0 == getSign(test));
		if(!check){
			out<<"Test for getSign function...FAILED"<<std::endl;
		}else{
			out<<"Test for getSign function...OK"<<std::endl;
		}
		
	}
	
	// posVectorFind
	{
		bool check=true;
		ivector1D test1(4,2);
		std::array<double,3> test2;
		test2.fill(-1.234);
		int mark1 = 2;
		double mark2 = -1.234;
		double mark3 = -2;
		int mark4 = 1;
		
		check = check && (posVectorFind(test1,mark1) == 0);
		check = check && (posVectorFind(test1,mark4) == -1);
		check = check && (posVectorFind(test2,mark2) == 0);
		check = check && (posVectorFind(test2,mark3) == -1);
		
		if(!check){
			out<<"Test for posVectorFind function...FAILED"<<std::endl;
		}else{
			out<<"Test for posVectorFind function...OK"<<std::endl;
		}
		
	}
	
	//checkVectorFind
	{
		bool check=true;
		ivector1D test1(4,2);
		std::array<double,8> test2;
		test2.fill(-1.234);
		int mark1 = 2;
		double mark2 = -1.234;
		double mark3 = -2;
		int mark4 = 1;
		
		check = check && (checkVectorFind(test1,mark1) == true);
		check = check && (checkVectorFind(test1,mark4) == false);
		check = check && (checkVectorFind(test2,mark2) == true);
		check = check && (checkVectorFind(test2,mark3) == false);
		
		if(!check){
			out<<"Test for checkVectorFind function...FAILED"<<std::endl;
		}else{
			out<<"Test for checkVectorFind function...OK"<<std::endl;
		}
		
	}
	
	//conVect
	{
		bool check = true;
		dvecarr3E testA(1);
		testA[0].fill(2);
		dvector2D testV = conVect(testA);
		for(int i=0; i<3; ++i){
			check = check && (testV[0][i] == testA[0][i]);
		}
		
		if(!check){
			out<<"Test for conVect function...FAILED"<<std::endl;
		}else{
			out<<"Test for conVect function...OK"<<std::endl;
		}
	}
	
	//conArray
	{
		bool check = true;
		dvector2D testV(1, dvector1D(3, 2));
		dvecarr3E testA = conArray<double,3>(testV);
		for(int i=0; i<3; ++i){
			check = check && (testV[0][i] == testA[0][i]);
		}
		
		if(!check){
			out<<"Test for conArray function...FAILED"<<std::endl;
		}else{
			out<<"Test for conArray function...OK"<<std::endl;
		}
	}
	
	//findPosition
	{
		bool check = true;
		dvector1D test(6,-1);
		double mark = 12;
		test[4] =mark;
		
		check = check && (4 == findPosition(mark, test));
		if(!check){
			out<<"Test for findPosition function...FAILED"<<std::endl;
		}else{
			out<<"Test for findPosition function...OK"<<std::endl;
		}
	}
	
	//getVectorSubset
	{
		bool check = true;
		dvector1D test(6,-1);
		dvector1D  result(2), result1(4,-1);
		result[0] = 8; result[1] = -2;
		
		test[3] = 8;
		test[4] = -2;
		
		//first test
		dvector1D prova = getVectorSubset(3,5,test);
		for(int i=0; i<prova.size(); ++i){
			check = check && (prova[i]==result[i]);
		}
		
		//second test
		dvector1D prova1 = getVectorSubset(5,3,test);
		for(int i=0; i<prova1.size(); ++i){
			check = check && (prova1[i]==result1[i]);
		}
		//third test
		dvector1D prova2 = getVectorSubset(3,3,test);
		for(int i=0; i<prova2.size(); ++i){
			int loc = (i+3)%prova2.size();
			check = check && (prova2[i]==test[loc]);
		}
		
		//fourth test	
		dvector1D prova3 = getVectorSubset(7,2,test);
		check = check && (prova3.size() == 0);
		
		//fifth test
		dvector1D prova4 = getVectorSubset(3,-1,test);
		check = check && (prova4.size() == 0);
	
		if(!check){
			out<<"Test for getVectorSubset function...FAILED"<<std::endl;
		}else{
			out<<"Test for getVectorSubset function...OK"<<std::endl;
		}
	}	

	//fillVectorSubset
	{
		bool check = true;
		dvector1D test(6,-1), dum(6,-1);
		dvector1D result(6,-1), result1(6,0), result2(6,0);
		result[3] = 0; result[4] = 0;
		result1[3] = -1; result1[4] = -1;
		
		//first test
		fillVectorSubset(3,5,test,0.0);
		for(int i=0; i<test.size(); ++i){
			check = check && (test[i]==result[i]);
		}
		
		test = dum;
		//second test
		fillVectorSubset(5,3,test,0.0);
		for(int i=0; i<test.size(); ++i){
			check = check && (test[i]==result1[i]);
		}
		
		test = dum;
		//second test
		fillVectorSubset(3,3,test,0.0);
		for(int i=0; i<test.size(); ++i){
			check = check && (test[i]==result2[i]);
		}
		
		test = dum;
		fillVectorSubset(7,2,test,0.0);
		for(int i=0; i<test.size(); ++i){
			check = check && (test[i]==dum[i]);
		}
		
		test = dum;
		fillVectorSubset(3,-1,test,0.0);
		for(int i=0; i<test.size(); ++i){
			check = check && (test[i]==dum[i]);
		}
		
		if(!check){
			out<<"Test for fillVectorSubset function...FAILED"<<std::endl;
		}else{
			out<<"Test for fillVectorSubset function...OK"<<std::endl;
		}
	}	
	
	

}else{
	std::cout<<"cannot open the log file for custom Operators"<<endl;
	
}

out.close();

return(0);
} //end test 
