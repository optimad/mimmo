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
\*---------------------------------------------------------------------------*/
#include "InOut.hpp"
#include "BaseManipulation.hpp"

using namespace std;


/*!
    Output stream operator for dvector1D 
    \param[in] buffer is the output stream
    \param[in] var is the element to be streamed
    \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const dvector1D &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        buffer << var[i];
    }
    return buffer;
}


/*!
    Input stream operator for dvector1D
    \param[in] buffer is the input stream
    \param[in] var is the element to be streamed
    \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, dvector1D &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> var[i];
    }
    return buffer;
}


/*!
* Output stream operator for livector1D
* \param[in] buffer is the output stream
* \param[in] var is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const livector1D &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        buffer << var[i];
    }
    return buffer;
}


/*!
* Input stream operator for livector1D
* \param[in] buffer is the input stream
* \param[in] var is the element to be streamed
* \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, livector1D &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> var[i];
    }
    return buffer;
}

/*!
* Output stream operator for livector2D
* \param[in] buffer is the output stream
* \param[in] var is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const livector2D &var)
{
    int nQ;
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        nQ = var[i].size();
        buffer << nQ;
        for (int j = 0; j < nQ; ++j) {
            buffer << var[i][j];
        }
    }
    return buffer;
}


/*!
* Input stream operator for livector2D
* \param[in] buffer is the input stream
* \param[in] var is the element to be streamed
* \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, livector2D &var)
{
    int nP, nQ;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> nQ;
        var[i].resize(nQ);
        for (int j = 0; j < nQ; ++j) {
            buffer >> var[i][j];
        }
    }
    return buffer;
}

/*!
* Output stream operator for shivector1D
* \param[in] buffer is the output stream
* \param[in] var is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const shivector1D &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        buffer << var[i];
    }
    return buffer;
}


/*!
* Input stream operator for shivector1D
* \param[in] buffer is the input stream
* \param[in] var is the element to be streamed
* \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, shivector1D &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> var[i];
    }
    return buffer;
}

/*!
    Output stream operator for dvecarr3E
    \param[in] buffer is the output stream
    \param[in] var is the element to be streamed
    \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const dvecarr3E &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 3; ++j) {
            buffer << var[i][j];
        }
    }

    return buffer;
}


/*!
    Input stream operator for dvecarr3E
    \param[in] buffer is the input stream
    \param[in] var is the element to be streamed
    \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, dvecarr3E &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 3; ++j) {
            buffer >> var[i][j];
        }
    }

    return buffer;
}

/*!
*	Output stream operator for mimmo::ShapeType enum
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const mimmo::ShapeType &var)
{
    buffer << static_cast<int> (var);
    return buffer;
}


/*!
*	Input stream operator for mimmo::ShapeType enum
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, mimmo::ShapeType &var)
{
    int val;
    buffer >> val;
    var = static_cast<mimmo::ShapeType>	(val);
    return buffer;
}

/*!
*	Output stream operator for std::array\<mimmo::CoordType,3\> enum
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const std::array<mimmo::CoordType,3> &var)
{
    std::array<int,3> dum;
    for(int i=0; i<3; ++i) dum[i] = static_cast<int> (var[i]);
    buffer << dum;
    return buffer;
}


/*!
*	Input stream operator for std::array\<mimmo::CoordType,3\> enum
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::array<mimmo::CoordType,3> &var)
{
    std::array<int,3> val;
    buffer >> val;
    for(int i=0; i<3; ++i)	var[i] = static_cast<mimmo::CoordType>(val[i]);
    return buffer;
}

/*!
*	Input stream operator for std::pair\<MimmoObject*, dvecarr3E *\> 
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::pair<mimmo::MimmoObject*, dvecarr3E *>& element){
    
    mimmo::MimmoObject * geo;
    dvecarr3E * data;
    buffer >> geo >> data ;
    element = std::make_pair(geo, data);
    return buffer;
};

/*!
*	Output stream operator for std::pair\<MimmoObject*, dvecarr3E *\>
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::pair<mimmo::MimmoObject*, dvecarr3E *>& element){
        buffer<<element.first<<element.second;
        return buffer;
};

/*!
*	Input stream operator for std::pair\<MimmoObject*, dvector1D *\> 
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::pair<mimmo::MimmoObject*, dvector1D *>& element){
    mimmo::MimmoObject * geo;
    dvector1D * data;
    buffer >> geo >> data ;
    element = std::make_pair(geo, data);
    return buffer;
};

/*!
*	Input stream operator for std::pair\<MimmoObject*, dvector1D *\> 
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::pair<mimmo::MimmoObject*, dvector1D *>& element){
    buffer<<element.first<<element.second;
    return buffer;
};


/*!
* Input stream operator for std::vector\<mimmo::TrackingPointer * \>
* \param[in] buffer is the input stream
* \param[in] var is the element to be streamed
* \result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector<mimmo::TrackingPointer * > & var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> var[i];
    }
    return buffer;
}

/*!
* Output stream operator for std::vector\<mimmo::TrackingPointer * \>
* \param[in] buffer is the output stream
* \param[in] var is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const std::vector<mimmo::TrackingPointer * > &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        buffer << var[i];
    }
    return buffer;
}

/*!
*	Input stream operator for std::pair\<BaseManipulation*, double\> 
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::pair<mimmo::BaseManipulation*, double>& element){
    mimmo::BaseManipulation * obj;
    double data;
    buffer >> obj >> data ;
    element = std::make_pair(obj, data);
    return buffer;
};

/*!
*	Input stream operator for std::pair\<BaseManipulation*, double\>
*	\param[in] buffer is the input stream
*	\param[in] element is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::pair<mimmo::BaseManipulation*, double>& element){
    buffer<<element.first<<element.second;
    return buffer;
};


/*!
*	Output stream operator for dvecarr2E
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const dvecarr2E &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 2; ++j) {
            buffer << var[i][j];
        }
    }
    
    return buffer;
}


/*!
*	Input stream operator for dvecarr2E
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, dvecarr2E &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 2; ++j) {
            buffer >> var[i][j];
        }
    }
    
    return buffer;
}

/*!
*	Output stream operator for ivecarr2E
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const ivecarr2E &var)
{
    int nP = var.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 2; ++j) {
            buffer << var[i][j];
        }
    }
    
    return buffer;
}


/*!
*	Input stream operator for ivecarr2E
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, ivecarr2E &var)
{
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        for (int j = 0; j < 2; ++j) {
            buffer >> var[i][j];
        }
    }
    return buffer;
}

/*!
*	Output stream operator for std::vector\<dvecarr2E\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream  &buffer, const std::vector<dvecarr2E> &var)
{
    int nP = var.size();
    int nP2;
    buffer << nP;
    for (int i = 0; i < nP; ++i) {
        nP2 = var[i].size();
        buffer << nP2;
        for (int j = 0; j < nP2; ++j) {
            for(int k=0; k<2; ++k){
                buffer << var[i][j][k];
            }
        }
    }
    
    return buffer;
}


/*!
*	Input stream operator for std::vector\<dvecarr2E\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector<dvecarr2E> &var)
{
    int nP, nP2;
    buffer >> nP;
    var.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> nP2;
        var[i].resize(nP2);
        for (int j=0; j<nP2; ++j){
            for (int k = 0; k < 2; ++k) {
                buffer >> var[i][j][k];
            }
        }	
    }
    return buffer;
}


/*!
*	Input stream operator for std::unordered_map\<std::string,std::pair\<int, MimmoObject*\> \>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::unordered_map<std::string, std::pair<int, mimmo::MimmoObject*> >&  var){
    
    std::string key;
    mimmo::MimmoObject * value;
    int ftype;
    int nP;
    buffer >> nP;
    for (int i = 0; i < nP; ++i) {
        buffer >> key>>ftype>>value;
        var[key] = std::make_pair(ftype,value);
    }
    return buffer;
};

/*!
*	Output stream operator for std::unordered_map\<std::string, std::pair\<int, MimmoObject*\>\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::unordered_map<std::string, std::pair<int, mimmo::MimmoObject*> >& var){
    
    int nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee.first<<ee.second.first<<ee.second.second;
    }
    return buffer;
};

/*!
*	Input stream operator for std::unordered_map\<long,std::pair\<int,long\> \>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::unordered_map<long,std::pair<int,long> >&  var){
    
    long key;
    int val1;
    long val2;
    int nP;
    buffer >> nP;
    for (int i = 0; i < nP; ++i) {
        buffer >> key>>val1>>val2;
        var[key] = std::make_pair(val1,val2);
    }
    return buffer;
};

/*!
*	Output stream operator for std::unordered_map\<long,std::pair\<int,long\> \>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::unordered_map<long,std::pair<int,long> >& var){
    
    int nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee.first<<ee.second.first<<ee.second.second;
    }
    return buffer;
};



/*!
*	Input stream operator for std::vector\<MimmoObject *\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::vector<mimmo::MimmoObject *>&  var){
    
    int nP;
    buffer >> nP;
    var.resize(nP);
    for (auto & ee : var) {
        buffer >> ee;
    }
    return buffer;
};

/*!
*	Output stream operator for std::vector\<MimmoObject *\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::vector<mimmo::MimmoObject *>& var){
    
    int nP = var.size();
    buffer << nP;
    for (auto & ee : var) {
        buffer << ee;
    }
    return buffer;
};



/*!
*	Input stream operator for mimmo::FileDataInfo
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, mimmo::FileDataInfo&  var){

    buffer >> var.ftype;
    buffer >> var.fdir;
    buffer >> var.fname;
    
    return buffer;
};

/*!
*	Output stream operator for mimmo::FileDataInfo
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const mimmo::FileDataInfo& var){
    
    buffer<<var.ftype<<var.fdir<<var.fname;
    return buffer;
    
    
};

/*!
*	Input stream operator for std::vector\<mimmo::FileDataInfo\>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector<mimmo::FileDataInfo>&  var){
    
    int nP;
    buffer>>nP;
    var.resize(nP);
    for(int i=0; i<nP; ++i){
        buffer>>var[i];
    }
    return buffer;
};

/*!
*	Output stream operator for std::vector\<mimmo::FileDataInfo\>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::vector<mimmo::FileDataInfo>& var){
    int nP = var.size();
    buffer<<nP;
    for(int i=0; i<nP; ++i){
        buffer<<var[i];
    }
    return buffer;
};

/*!
*	Input stream operator for std::unordered_map\< MimmoObject*, dvector1D* \>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::unordered_map< mimmo::MimmoObject*, dvector1D* >&  var){
    
    mimmo::MimmoObject * key;
    dvector1D * val;
    int size;
    buffer>>size;
    for (int i=0; i<size; ++i) {
        buffer >> key >> val;
        var[key] = val;
    }
    return buffer;
};

/*!
*	Output stream operator for std::unordered_map\< MimmoObject*, dvector1D* \>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::unordered_map< mimmo::MimmoObject*, dvector1D* >& var){
    int size = var.size();
    buffer<<size;
    for (auto & ee : var) {
        buffer << ee.first<<ee.second;
    }
    return buffer;
};

/*!
*	Input stream operator for std::unordered_map\< MimmoObject*, dvecarr3E* \>
*	\param[in] buffer is the input stream
*	\param[in] var is the element to be streamed
*	\result Returns the same input stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer,std::unordered_map< mimmo::MimmoObject*, dvecarr3E* >&  var){
    
    mimmo::MimmoObject * key;
    dvecarr3E * val;
    int size;
    buffer>>size;
    for (int i=0; i<size; i++) {
        buffer >> key >> val;
        var[key] = val;
    }
    return buffer;
};

/*!
*	Output stream operator for std::unordered_map\< MimmoObject*, dvecarr3E* \>
*	\param[in] buffer is the output stream
*	\param[in] var is the element to be streamed
*	\result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::unordered_map< mimmo::MimmoObject*, dvecarr3E* >& var){
    int size = var.size();
    buffer<<size;
    for (auto & ee : var) {
        buffer << ee.first<<ee.second;
    }
    return buffer;
};

/*!
* Input stream operator for std::vector\< std::pair\<mimmo::MimmoObject*, dvector1D *\> \>
* \param[in] buffer is the input stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector< std::pair<mimmo::MimmoObject*, dvector1D *> >& element){
    
    int nP;
    buffer >> nP;
    element.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> element[i];
    }
    return buffer;
};

/*!
* Output stream operator for std::vector\< std::pair\<mimmo::MimmoObject*, dvector1D *\> \>
* \param[in] buffer is the output stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::vector< std::pair<mimmo::MimmoObject*, dvector1D *> >& element){
    
    int nP = element.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) buffer << element[i];
    return buffer;
};

/*!
* Input stream operator for std::vector\< std::pair\<mimmo::MimmoObject*, dvecarr3E *\> \>
* \param[in] buffer is the input stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &buffer, std::vector< std::pair<mimmo::MimmoObject*, dvecarr3E *> >& element){
    
    int nP;
    buffer >> nP;
    element.resize(nP);
    for (int i = 0; i < nP; ++i) {
        buffer >> element[i];
    }
    return buffer;
};

/*!
* Output stream operator for std::vector\< std::pair\<mimmo::MimmoObject*, dvecarr3E *\> \>
* \param[in] buffer is the output stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &buffer, const std::vector< std::pair<mimmo::MimmoObject*, dvecarr3E *> >& element){
    
    int nP = element.size();
    buffer << nP;
    for (int i = 0; i < nP; ++i) buffer << element[i];
    return buffer;
};

//==============================================================//
// DATA TYPE  CLASS	IMPLEMENTATION								//
//==============================================================//
namespace mimmo{

/*!
 * Default constructor of DataType
 */
DataType::DataType(){
    m_conType = containerTAG::SCALAR;
    m_dataType = dataTAG::INT;
};

/*!
 * Custom constructor of DataType.
* \param[in] conType TAG of type of container.
* \param[in] dataType TAG of type of data.
*/
DataType::DataType(containerTAG conType, dataTAG dataType){
    m_conType 	= conType;
    m_dataType	= dataType;
    return;
};

/*!
 * Default destructor of DataType
 */
DataType::~DataType(){};

/*!
 * Copy constructor of DataType.
 */
DataType::DataType(const DataType & other){
    m_conType 	= other.m_conType;
    m_dataType	= other.m_dataType;
    return;
};

/*!
 * Compare operator of DataType.
 */
bool DataType::operator==(const DataType & other){
    bool check = true;
    check = check && (m_conType == other.m_conType);
    check = check && (m_dataType == other.m_dataType);
    return(check);
};


//==============================================================//
// BASE INOUT CLASS	IMPLEMENTATION                              //
//==============================================================//

/*!
 * Default constructor of PortOut
 */
PortOut::PortOut(){
    m_objLink.clear();
};

/*!
 * Default destructor of PortOut
 */
PortOut::~PortOut(){};

/*!
 * Copy constructor of PortOut.
 */
PortOut::PortOut(const PortOut & other){
    m_objLink 	= other.m_objLink;
    m_obuffer	= other.m_obuffer;
    m_portLink	= other.m_portLink;
    m_datatype	= other.m_datatype;
    return;
};

/*!
 * Compare operator of PortOut.
 */
bool PortOut::operator==(const PortOut & other){
    bool check = true;
    check = check && (m_portLink == other.m_portLink);
    check = check && (m_objLink == other.m_objLink);
    check = check && (m_datatype == other.m_datatype);
    return(check);
};

/*!
 * It gets the objects linked by this port.
* \return vector of pointer to linked objects.
*/
std::vector<mimmo::BaseManipulation*>
PortOut::getLink(){
    return(m_objLink);
}

/*!
 * It gets the input port ID of the objects linked by this port.
* \return Vector of PortID.
*/
std::vector<PortID>
PortOut::getPortLink(){
    return(m_portLink);
}

/*!
 * It gets the TAG of data type communicated by this port.
* \return TAG of datat type communicated.
*/
DataType
PortOut::getDataType(){
    return(m_datatype);
}

/*!
 * It empties the output buffer.
 */
void
mimmo::PortOut::cleanBuffer(){
    m_obuffer.seekg(0);
}

/*!
 * It clears the links to objects and the related ports.
 */
void
mimmo::PortOut::clear(){
    m_objLink.clear();
    m_portLink.clear();
}

/*!
 * It removes the link to an object and the related port ID.
* \param[in] j Index of the linked object in the links vector of this port.
*/
void
mimmo::PortOut::clear(int j){
    if (j < (int)m_objLink.size() && j >= 0){
        m_objLink.erase(m_objLink.begin() + j);
        m_portLink.erase(m_portLink.begin() + j);
    }
}

/*! 
 * Execution of the PIN.
 * All the pins are called in execution of the sending owner after its own execution.
 * Reading stage of pin linked receivers is automatically performed within this execution.
 */
void
mimmo::PortOut::exec(){
    if (m_objLink.size() > 0){
        writeBuffer();
        bitpit::IBinaryStream input(m_obuffer.data(), m_obuffer.getSize());
        cleanBuffer();
        for (int j=0; j<(int)m_objLink.size(); j++){
            if (m_objLink[j] != NULL){
                m_objLink[j]->setBufferIn(m_portLink[j], input);
                m_objLink[j]->readBufferIn(m_portLink[j]);
                m_objLink[j]->cleanBufferIn(m_portLink[j]);
            }
        }
    }
};

/*!
 * Default constructor of PortIn
 */
PortIn::PortIn(){
    m_mandatory =false;
    m_familym = 0;
};

/*!
 * Default destructor of PortIn
 */
PortIn::~PortIn(){};

/*!
 * Copy constructor of PortIn.
*/
PortIn::PortIn(const PortIn & other){
    m_objLink 	= other.m_objLink;
    m_ibuffer	= other.m_ibuffer;
    m_datatype  = other.m_datatype;
    m_mandatory = other.m_mandatory;
    m_familym   = other.m_familym;
    return;
};

/*!
 * Compare operator of PortIn.
 */
bool PortIn::operator==(const PortIn & other){
    bool check = true;
    check = check && (m_objLink == other.m_objLink);
    check = check && (m_datatype == other.m_datatype);
    check = check && (m_mandatory = other.m_mandatory);
    check = check && (m_familym   = other.m_familym);
    return(check);
};

/*!
 * It gets the linked object by this port.
 * \return Pointer to linked object.
 */
std::vector<mimmo::BaseManipulation*>
PortIn::getLink(){
    return(m_objLink);
}

/*!
 * It gets the TAG of data type communicated by this port.
 * \return TAG of datat type communicated.
 */
DataType
PortIn::getDataType(){
    return(m_datatype);
}

/*!
 * It gets if this port has to be mandatorily linked.
 * \return mandatory?.
 */
bool
PortIn::isMandatory(){
    return(m_mandatory);
}

/*!
 * It gets the family of this port.
 * \return mandatory family.
 */
int
PortIn::getFamily(){
    return(m_familym);
}

/*!
 * It clears all linked objects into this port.
 */
void
mimmo::PortIn::clear(){
    m_objLink.clear();
}

/*!
 * It removes the link to an object and the related port ID.
 * \param[in] j Index of the linked object in the links vector of this port.
 */
void
mimmo::PortIn::clear(int j){
    if (j < (int)m_objLink.size() && j >= 0){
        m_objLink.erase(m_objLink.begin() + j);
    }
}


/*!
 * It releases the memory occupied by the input buffer.
 */
void
mimmo::PortIn::cleanBuffer(){
    m_ibuffer.seekg(0);
}

}
