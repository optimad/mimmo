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
\*---------------------------------------------------------------------------*/
#ifndef __MIMMOSHAREDPOINTER_HPP__
#define __MIMMOSHAREDPOINTER_HPP__

#include <mimmo_binary_stream.hpp>

//forward declaration.
namespace mimmo{
    template<class O>
    class MimmoSharedPointer;
}

/*!
 * \ingroup binaryStream
 * \{
 */
// Providing binary stream for class MimmoSharedPointer
template< class T>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, mimmo::MimmoSharedPointer<T>& element);

template< class T >
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const mimmo::MimmoSharedPointer<T>& element);
/*!
 *\}
 */

namespace mimmo{

/*!
   @class MimmoSharedPointer
   \ingroup core
   @brief MimmoSharedPointer is a custom implementation of shared pointer

  MimmoSharedPointer is a custom implementation of shared pointer of mimmo library.
  It is basically a smart pointer who shares the ownership of a base pointer of a
  dynamically allocated class O (using c++ new) throughout its own multiple instances.
  Instances multiplicity is tracked incrementing/decrementing an internal counter.
  Each copy assignment/construction of a shared pointer increment the counter
  with respect to the copy source.
  A move assignment/construction does not increment counter, but move the contents
  of the source as they are to the new shared pointer; hence the source becomes a void object.
  When the last instance of the smart pointer is deleted (counter to 0), the
  internal object pointed is deleted too.
  Please beware, once a pointer to an allocated object is delegated to the shared
  pointer, the last one is responsible of the object destruction. If a forced external
  destruction occurs with a direct delete on the raw pointer object, it leads to
  memory violations inside the shared pointer.\n
  Example of correct/incorrect usage of mimmo shared pointer are reported below.
  Please note, despite the examples don't due to the aim to highlight the incorrectness
  of usage, that it is a best practice to allocate dynamically using c++ new directly in a
  MimmoSharedPointer constructor.\n


  ---------------------
  <CODE>
  <B>SAMPLE 1 CORRECT</B>\n
  mimmo::MimmoSharedPointer<std::string> hold_ptr; //create a new empty pointer.\n
  {\n
  &nbsp; std::string * data = new std::string("haveSomeStringDataHere"); //allocate dynamically some string data\n
  &nbsp; mimmo::MimmoSharedPointer<std::string> create_ptr(data); //use data pointer to create a mimmo shared pointer with counter 1.\n
  &nbsp; hold_ptr = create_ptr; //share data instance with another shared pointer. Now counter is set to 2;\n
  }\n
   (after exiting scope create_ptr is destroyed, counter return to 1.\n
    Now hold_ptr owns data instance with counter 1 and it will worry to destroy data,\n
   when it destroys himself.)\n
  </CODE>


  -----------------------
  <CODE>
  <B>SAMPLE 1 WRONG</B>\n
  mimmo::MimmoSharedPointer<std::string> hold_ptr; //create a new empty pointer.\n
  {\n
  &nbsp; std::string data("haveSomeStringDataHere"); //allocate statically some string data\n
  &nbsp; mimmo::MimmoSharedPointer<std::string> create_ptr(&data); // use data pointer to create a mimmo shared pointer with counter 1.\n
  &nbsp; hold_ptr = create_ptr; // share data instance with another shared pointer. Now counter is set to 2;\n
  }\n
   (after exiting scope create_ptr is destroyed, counter return to 1,\n
    but string data allocated statically are destroyed when out of scope!\n
    Now hold_ptr has counter 1 and holds a data pointer pointing to an already destroyed object\n
   and when it tries to destroy himself(and data instance) this will lead to a memory violation).\n
  </CODE>


  -------------------------
  <CODE>
  <B>SAMPLE 2 CORRECT</B>\n
  std::string * data = new std::string("haveSomeStringDataHere"); //allocate dynamically some string data\n
  mimmo::MimmoSharedPointer<std::string> ptr1(data);//pass it to a mimmo shared pointer using base constructor\n
  {\n
  &nbsp; mimmo::MimmoSharedPointer<std::string> ptr2(ptr1); //pass the same data to another pointer usign base constructor\n
  &nbsp; (create a safe copy of the pointer using always copy constructor or assignment,\n
  &nbsp; that increment counter to 2.)\n
  }\n
  (after exiting ptr2 is destroyed and counter is decremented by 1, but data are still safe.\n
  ptr1 is still alive with counter = 1 and a pointer pointing safely to data).\n
  </CODE>


  -------------------------
  <CODE>
  <B>SAMPLE 2 WRONG</B>\n
  std::string * data = new std::string("haveSomeStringDataHere"); //allocate dynamically some string data\n
  mimmo::MimmoSharedPointer<std::string> ptr1(data);//pass it to a mimmo shared pointer using base constructor\n
  {\n
  &nbsp; mimmo::MimmoSharedPointer<std::string> ptr2(data);//pass the same data to another pointer using base constructor\n
  &nbsp; (both ptr1 and ptr2 owns the same data pointer with counter 1)\n
  }\n
  (after exiting ptr2 is destroyed and with him also data (ptr2 counter reached 0):\n
   but ptr1 is still alive with counter 1 and a pointer pointing to a data already destroyed.\n
   Destroying ptr1 will lead to a memory violation).\n
  </CODE>


   -------------------------
 */
template<class O>
class MimmoSharedPointer
{

template< class T>
friend mimmo::IBinaryStream& (::operator>>) (mimmo::IBinaryStream &buffer, mimmo::MimmoSharedPointer<T>& element);

template< class T >
friend mimmo::OBinaryStream& (::operator<<) (mimmo::OBinaryStream &buffer, const mimmo::MimmoSharedPointer<T>& element);

public:

    MimmoSharedPointer(O* object = nullptr);
    ~MimmoSharedPointer();
    MimmoSharedPointer(const MimmoSharedPointer<O> & other);
    MimmoSharedPointer(MimmoSharedPointer<O> && other);
    MimmoSharedPointer & operator=(const MimmoSharedPointer<O> & other);
    MimmoSharedPointer & operator=(MimmoSharedPointer<O> && other);

    //methods
    void  swap(MimmoSharedPointer<O>& other) noexcept;
    O* get() const;
    std::size_t getCounter() const;

    void reset(O* object = nullptr);

    //derefencing operators
    O& operator*() const;
    O* operator->() const;

    //boolean operators
    bool operator!() const;
    operator bool() const;

    //inline compare operators
    /*!
     * Equal to operator of MimmoSharedPointer. It compares only pointer to object.
     * \param[in] rhs Tested MimmoSharedPointer
     */
    inline bool operator==(const mimmo::MimmoSharedPointer<O>& rhs) const
    {
        return (get() == rhs.get());
    };

    /*!
     * Not equal to operator of MimmoSharedPointer. It compares only pointer to object.
     * \param[in] rhs Tested MimmoSharedPointer
     */
    inline bool operator!=(const mimmo::MimmoSharedPointer<O>& rhs) const
    {
        return !(*this == rhs);
    }

    /*!
     * Equal operator of MimmoSharedPointer to pure pointer to object type.
     * \param[in] _rhs Tested pointer to object
     */
    inline bool operator==(const O* _rhs) const
    {
        // Check only pointer to object
        return (m_object == _rhs);
    };

    /*!
     * Not equal operator of MimmoSharedPointer to pure pointer to object type.
     * \param[in] _rhs Constant tested pointer to object
     */
    inline bool operator!=(const O* _rhs) const
    {
        return !(*this == _rhs);
    }

    /*!
     * Equal operator of MimmoSharedPointer to void pointer.
     * \param[in] _rhs Tested void pointer
     */
    inline bool operator==(const void* _rhs) const
    {
        // Check only pointer to object
        return (m_object == _rhs);
    };

    /*!
     * Not equal operator of MimmoSharedPointer to void pointer.
     * \param[in] _rhs Tested void pointer
     */
    inline bool operator!=(const void* _rhs) const
    {
        return !(*this == _rhs);
    }

    /*!
     * Equal operator of MimmoSharedPointer to null pointer.
     * \param[in] _rhs Tested null pointer
     */
    inline bool operator==(const std::nullptr_t _rhs) const
    {
        // Check only pointer to object
        return (m_object == _rhs);
    };

    /*!
     * Not equal operator of MimmoSharedPointer to null pointer.
     * \param[in] _rhs Tested null pointer
     */
    inline bool operator!=(const std::nullptr_t _rhs) const
    {
        return !(*this == _rhs);
    }

    /*!
     * Equal operator of MimmoSharedPointer to long int.
     * \param[in] lrhs Tested long int
     */
    inline bool operator==(const long int lrhs) const
    {
        // Check only pointer to object
        if (lrhs == 0){
            return (!m_object);
        }
        else{
            return bool(m_object);
        }
    };

    /*!
     * Not equal operator of MimmoSharedPointer to long int.
     * \param[in] lrhs Tested long int
     */
    inline bool operator!=(const long int lrhs) const
    {
        return !(*this == lrhs);
    }


protected:

    void _init(O* object = nullptr, std::size_t* counter = nullptr);
    void _reset(O* object = nullptr, std::size_t* counter = nullptr);
    void _decrement();

    O*              m_object;   /**<Pointer to target objet. */
    std::size_t*    m_counter;  /**<Pointer to counter of MimmoSharedPointer objects
                                    sharing the pointer to target.*/

};

} // end namespace mimmo


namespace std{

/*!
 * Hash structure for MimmoSharedPointer class.
 */
template<class O>
struct hash<mimmo::MimmoSharedPointer<O> >
{
    size_t
    operator()(const mimmo::MimmoSharedPointer<O> & obj) const;

    size_t
    operator()(mimmo::MimmoSharedPointer<O> & obj);
};

} // end namespace std

#include "MimmoSharedPointer.tpp"

#endif /* __MIMMOSHAREDPOINTER_HPP__ */
