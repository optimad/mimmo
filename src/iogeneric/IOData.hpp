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
#ifndef __IODATA_HPP__
#define __IODATA_HPP__

namespace mimmo{

/*!
 *
 * \class IOData
 * \ingroup iogeneric
 * \brief IOData is the base class of generic data stored as input or result in a Port.
 *
 * IOData represents a generic data used in manipulation objects. Get and set methods of
 * this class are templated methods on the type of attached data.
 *
 */
class IOData{

public:
    /*!Default constructor of IOData.
     */
    IOData(){};

    /*!
     * Destructor
     */
    virtual ~IOData(){};

    /*!It gets the data stored in the object.
     * Even if not declared as pure virtual its behavior is
     * like a pure virtual method, i.e. an analogous method
     * is implemented in the template derived class IODataT.
     * \return Pointer to data stored.
     */
    template<typename T>
    T* getData();

    /*!It sets the data stored in the object.
     * Even if not declared as pure virtual its behavior is
     * like a pure virtual method, i.e. an analogous method
     * is implemented in the template derived class IODataT.
     * \param[in] data Data to be stored.
     */
    template<typename T>
    void setData(T data);

};

/*!
 * \class IODataT
 * \ingroup iogeneric
 * \brief IODataT is the templated class of generic data derived from IOData base class.
 *
 * IODataT stores a generic data used in manipulation objects.
 */
template<typename T>
class IODataT: public IOData{
public:
    T                 m_data;    /**<Data contained in the object.*/

public:
    /*!Default constructor of IODataT.
     */
    IODataT(){};

    /*!Custom constructor of IODataT.
     * \param[in] data Data to be stored.
     */
    IODataT(T &data){
        m_data = data;
    };

    /*!Default destructor of IODataT.
     */
    ~IODataT(){};

    /*!Copy constructor of IODataT.
     */
    IODataT(const IODataT & other){
        this->m_data     = other.m_data;
    }

    /*!It sets the data stored in the object.
     * \param[in] data Data to be stored.
     */
    void setData(T &data){
        m_data = data;
    }

    /*!It gets the data stored in the object.
     * \return Pointer to data stored.
     */
    T* getData(){
        return(&m_data);
    }

};

};


#endif /* __IODATA_HPP__ */
