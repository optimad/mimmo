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

#include <mimmo_common.hpp>
#include "MimmoSharedPointer.hpp"
#include <unordered_set>

/*
 * Test 6_1
 * Testing MimmoSharedPointer base construction, dereferencing
 */
// =================================================================================== //
int test6_1() {

    std::vector<int> * v1 = new std::vector<int>(3, -4);
    mimmo::MimmoSharedPointer<std::vector<int>> sp_v1(v1);

    bool check = sp_v1.getCounter() == 1;
    check = check && ( sp_v1->size() == (*sp_v1).size() );
    check = check && ( sp_v1->size() == 3);
    std::cout<<"mimmo shared pointer base construction and dereferencing test ";
    if(check)
        std::cout<<"...PASSED"<<std::endl;
    else
        std::cout<<"...FAILED"<<std::endl;
    return !check;
}
// =================================================================================== //

/*
 * Test 6_2
 * Testing MimmoSharedPointer swap and reset method
 */
// =================================================================================== //
int test6_2() {

    std::vector<int> * v1 = new std::vector<int>(6, -4);
    std::vector<int> * v2 = new std::vector<int>(3, 12);
    std::vector<int> * v3 = new std::vector<int>(10, 632);

    mimmo::MimmoSharedPointer<std::vector<int>> sp_v1(v1);
    mimmo::MimmoSharedPointer<std::vector<int>> sp_v2(v2);

    sp_v1.swap(sp_v2);

    bool check = (sp_v1->size()) == 3 && (sp_v2->size()) == 6 ;

    //augmenting reference of sp_v1
    mimmo::MimmoSharedPointer<std::vector<int>> copyspv1(sp_v1);
    //reset sp_v1 to v3. This will decrement reference of old contents (see copyspv1 must be to 1)
    // and reset sp_v1 to new v3 content with counter 1.
    sp_v1.reset(v3);

    check = check && ( sp_v1->size() == 10 && sp_v1.getCounter() == 1);
    check = check && ( copyspv1.getCounter() == 1);
    std::cout<<"mimmo shared pointer swap and reset test ";
    if(check)
        std::cout<<"...PASSED"<<std::endl;
    else
        std::cout<<"...FAILED"<<std::endl;
    return !check;
}

// =================================================================================== //
/*
 * Test 6_3
 * Testing MimmoSharedPointer Copy and Move Assignments.
 */
// =================================================================================== //
int test6_3() {

    std::vector<int> * v1 = new std::vector<int>(6, -4);
    std::vector<int> * v2 = new std::vector<int>(3, 12);

    //create 2 shared pointer with different multiplicity instance.
    //sp_v1 with counter 4.
    mimmo::MimmoSharedPointer<std::vector<int>> sp_v1(v1);
    mimmo::MimmoSharedPointer<std::vector<int>> c1_sp_v1(sp_v1);
    mimmo::MimmoSharedPointer<std::vector<int>> c2_sp_v1(sp_v1);
    mimmo::MimmoSharedPointer<std::vector<int>> c3_sp_v1(sp_v1);

    //sp_v2 with counter 2.
    mimmo::MimmoSharedPointer<std::vector<int>> sp_v2(v2);
    mimmo::MimmoSharedPointer<std::vector<int>> c1_sp_v2(sp_v2);

    bool check=true;
    //Copy Assignment testing
    {
        mimmo::MimmoSharedPointer<std::vector<int>> test;
        test = sp_v1; // this should have counter = 5;
        check = check && ( test.getCounter() == 5);

        test = sp_v2; // this should have counter = 3;
        check = check && ( test.getCounter() == 3);
    }
    check = check && ( sp_v1.getCounter() == 4);
    check = check && ( sp_v2.getCounter() == 2);

    //Move Assignment testing
    //test1 substitute sp_v1 without incrementing counter. sp_v1 will become empty
    mimmo::MimmoSharedPointer<std::vector<int>> test1 = std::move(sp_v1);
    check = check && ( test1.getCounter() == 4);
    check = check && ( sp_v1.get() == nullptr);

    //test1 substituting sp_v2 without incrementing counter. sp_v2 will become empty.
    //old content owned still shared with c1_sp_v1, c2_sp_v1,c3_sp_v1 will be decreased by 1.
    test1 = std::move(sp_v2);
    check = check && ( test1.getCounter() == 2);
    check = check && ( sp_v2.get() == nullptr);
    check = check && ( c1_sp_v1.getCounter() == 3);

    std::cout<<"mimmo shared pointer copy and move assignment test ";
    if(check)
        std::cout<<"...PASSED"<<std::endl;
    else
        std::cout<<"...FAILED"<<std::endl;
    return !check;
}

// =================================================================================== //

/*
    class for testing in 6_4;
*/
class CustomClass{
    public:
        CustomClass(){}
        void setInstance(mimmo::MimmoSharedPointer<std::vector<int>> instance)
        {
            m_instance = instance;
        }
        mimmo::MimmoSharedPointer<std::vector<int>> getInstance()
        {
            return m_instance;
        }
        std::size_t getInternalInstanceCounter()
        {
            return m_instance.getCounter();
        }
    private:
        mimmo::MimmoSharedPointer<std::vector<int>> m_instance;
};

/*
 * Test 6_4
 * Testing MimmoSharedPointer exchange through Class-set/get Methods and Binary streams.
 */
// =================================================================================== //
int test6_4() {
    mimmo::MimmoSharedPointer<std::vector<int>> targetmsp(new std::vector<int>(4, -1));
    bool check = true;
    std::size_t beforeCounter, afterCounter;

    //binary stream passing just the target shared pointer. The expected result is
    // the targetmsp target counter incremented by 1.
    beforeCounter = targetmsp.getCounter();
    {
        mimmo::OBinaryStream obuffer;
        obuffer<<targetmsp;
        //creating a input binary stream from the previous output buffer
        mimmo::IBinaryStream ibuffer(obuffer.data(), obuffer.getSize());
        // use the input binary stream to fill another shared pointer.
        mimmo::MimmoSharedPointer<std::vector<int>> fromStream;
        ibuffer>>fromStream;
        afterCounter = fromStream.getCounter();
    }
    check = check && ( (afterCounter - beforeCounter) == 1);

    //binary stream passing target shared pointer through get method of a class. The expected result is
    // the targetmsp target counter incremented by 1.
    CustomClass p1, p2;
    //using a set bu msp copy //this will increment the counter of targetmsp of 1
    p1.setInstance(targetmsp);
    beforeCounter = targetmsp.getCounter();

    //Try through the class
    {
        // write msp data into a buffer
        mimmo::OBinaryStream obuffer;
        {
            auto temp = p1.getInstance();
            obuffer<<temp;
        }
        //creating a input binary stream from the previous output buffer
        mimmo::IBinaryStream ibuffer(obuffer.data(), obuffer.getSize());
        // use the input binary stream to fill another shared pointer (works as a copy construction, incrementing of 1).
        mimmo::MimmoSharedPointer<std::vector<int>> fromStream;
        ibuffer>>fromStream;
        p2.setInstance(fromStream);
    }
    afterCounter = targetmsp.getCounter();
    check = check && ( (afterCounter - beforeCounter) == 1);

    //I have now tre instances shared in targetmsp, p1,p2.
    check = check && (targetmsp.getCounter()== 3);

    std::cout<<"mimmo shared pointer binary streams exchange test ";
    if(check)
        std::cout<<"...PASSED"<<std::endl;
    else
        std::cout<<"...FAILED"<<std::endl;
    return !check;
}
// =================================================================================== //

/*
 * Test 6_5
 * Testing MimmoSharedPointer compare operators and hashing.
 */
// =================================================================================== //
int test6_5() {

    bool check = true;
    std::unordered_map<mimmo::MimmoSharedPointer<std::string>, int> hashmap;

    {
        std::unordered_map<int, mimmo::MimmoSharedPointer<std::string>> umap;
        umap[0]=std::move(mimmo::MimmoSharedPointer<std::string>(new std::string("Stannis")));
        umap[1]=std::move(mimmo::MimmoSharedPointer<std::string>(new std::string("Boris")));
        umap[2]=std::move(mimmo::MimmoSharedPointer<std::string>(new std::string("Rene")));
        umap[3]=std::move(mimmo::MimmoSharedPointer<std::string>(umap[1]));
        umap[4]=std::move(mimmo::MimmoSharedPointer<std::string>(umap[2]));

        //testing boolean op
        bool check = umap[0];

        //testing compare op
        check = check && (umap[1] == umap[3]);
        check = check && (umap[1] != umap[2]);

        //testing hashing and compare through unordered map.
        for(auto & tuple: umap){
            hashmap[tuple.second] = tuple.first;
        }
    }
    //the inverse hashmap now will have 3 elements, and all of them will count 1
    //because umap is destroyed.
    check = check && (hashmap.size() == 3);

    for(auto & tuple: hashmap){
        check = check && (tuple.first.getCounter() == 1);
    }

    std::cout<<"mimmo shared pointer compare ops and hashing test ";
    if(check)
        std::cout<<"...PASSED"<<std::endl;
    else
        std::cout<<"...FAILED"<<std::endl;
    return !check;
}


int main( int argc, char *argv[] ) {

	BITPIT_UNUSED(argc);
	BITPIT_UNUSED(argv);

#if MIMMO_ENABLE_MPI
    MPI_Init(&argc, &argv);
#endif

	int val = 1;
    /**<Calling mimmo Test routines*/
    try{
        val = test6_1() ;
        val = std::max(val, test6_2());
        val = std::max(val, test6_3());
        val = std::max(val, test6_4());
        val = std::max(val, test6_5());
    }
    catch(std::exception & e){
        std::cout<<"test_core_00006 exited with an error of type : "<<e.what()<<std::endl;
        return 1;
    }

#if MIMMO_ENABLE_MPI
	MPI_Finalize();
#endif

	return val;
}
