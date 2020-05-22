/*---------------------------------------------------------------------------*\
 *
 *  mimmo
 *
 *  Copyright (C) 2015-2020 OPTIMAD engineering Srl
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


/*!
* Input stream operator for bitpit::MimmoSharedPointer\< T \>
* \param[in] buffer is the input stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
template< class T>
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buffer, mimmo::MimmoSharedPointer<T>& element){
    T* object;
    std::size_t* counter;
    buffer >> object;
    buffer >> counter;
    element._reset(object, counter);
    return buffer;
};

/*!
* Output stream operator for bitpit::MimmoSharedPointer\< T \>
* \param[in] buffer is the output stream
* \param[in] element is the element to be streamed
* \result Returns the same output stream received in input.
*/
template< class T >
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buffer, const mimmo::MimmoSharedPointer<T>& element)
{
    buffer << element.m_object;
    buffer << element.m_counter;
    return buffer;
};


namespace mimmo{

/*!
   Default constructor of MimmoSharedPointer.
   The class must be initialized with a pointer to a O type class, dynamically
   allocated. If a valid pointer is passed, the class will set its counter
   to 1.
   \param[in] object raw pointer to a target object
 */
template<class O>
MimmoSharedPointer<O>::MimmoSharedPointer(O* object)
{
    _init(object, nullptr);
}

/*!
   Destructor. When a MimmoSharedPointer is destroyed the internal counter
   is decreased by one. If the counter reaches zero,
   the pointed object is explicitly deleted.
 */
template<class O>
MimmoSharedPointer<O>::~MimmoSharedPointer(){
    _decrement();
}

/*!
  Copy constructor of MimmoSharedPointer. Counter is increased by one
  with respect the copy source.
  \param[in] other Copied MimmoSharedPointer
 */
template<class O>
MimmoSharedPointer<O>::MimmoSharedPointer(const MimmoSharedPointer<O> & other)
{
    m_object = other.m_object;
    m_counter = other.m_counter;
    if (m_counter){
        (*m_counter)++;
    }
}

/*!
  Move constructor of MimmoSharedPointer. Contents are moved from the target source
  to the current class, without incrementing internal counter. Target source becomes
  an empty shell, ready to be destroyed, without affecting its original contents.
  \param[in] other Source MimmoSharedPointer to be moved
 */
template<class O>
MimmoSharedPointer<O>::MimmoSharedPointer(MimmoSharedPointer<O> && other)
{
    //attach the new contents to the class
    m_object = other.m_object;
    m_counter = other.m_counter;

    //empty the target source contents. They are owned by this class now
    other.m_object = nullptr;
    other.m_counter = nullptr;
}

/*!
   Copy Assignment operator of MimmoSharedPointer. First, original contents of the class
   are released decrementing the relative counter of 1 (destruction of contents occurs
   only if such counter reaches 0). Then new contents from target source are copied into
   the current class, incrementing the counter by one.
   \param[in] other Copied MimmoSharedPointer source
*/
template<class O>
MimmoSharedPointer<O> & MimmoSharedPointer<O>::operator=(const MimmoSharedPointer<O> & other)
{
    //check to avoid useless operations in case of auto-copy assignment.
    if(&other == this)
        return *this;

    //decrement the old contents in the class(if any).
    _decrement();

    //assign new contents from source
    m_object = other.m_object;
    m_counter = other.m_counter;

    //incrementing my counter.
    if (m_counter){
        (*m_counter)++;
    }
    return *this;
}

/*!
   Move Assignment operator of MimmoSharedPointer. First, original contents of the class
   are released decrementing the relative counter of 1 (destruction of contents occurs
   only if such counter reaches 0). Then, new contents are moved from the target source
   to the current class, without incrementing their internal counter. Finally, target source becomes
   an empty shell, ready to be destroyed, without affecting the contents.
   \param[in] other Source MimmoSharedPointer to be moved
*/
template<class O>
MimmoSharedPointer<O> & MimmoSharedPointer<O>::operator=(MimmoSharedPointer<O> && other)
{
    //check to avoid auto-move assignment.
    if(&other == this)
        return *this;

    //decrement the old contents in the class(if any).
    _decrement();

    //attach the new contents to the class
    m_object = other.m_object;
    m_counter = other.m_counter;

    //empty the target source contents. They are owned by this class now
    other.m_object = nullptr;
    other.m_counter = nullptr;

    return *this;
}

/*!
 * Swap MimmoSharedPointer method. Swap contents between current class and
   another one of the same type, keeping both reference counting untouched.
 * \param[in] other Target shared pointer to swap
 */
 template<class O>
 void MimmoSharedPointer<O>::swap(MimmoSharedPointer<O>& other) noexcept
{
    std::swap(m_object, other.m_object);
    std::swap(m_counter, other.m_counter);
}


/*!
   \return raw pointer to internal object controlled by the class
 */
template<class O>
O* MimmoSharedPointer<O>::get() const
{
    return m_object;
}

/*!
   Return the number of the multiple instances sharing the current pointed objects.
   In case of the raw internal pointer is nullptr, return 0.
   \return internal counter status
 */
template<class O>
std::size_t MimmoSharedPointer<O>::getCounter() const
{
    if(m_counter)   return *(m_counter);
    else            return 0;
}


/*!
 * Reset MimmoSharedPointer to a new target pointer or to nullptr.
   Old contents of the class are released decrementing their counter by one
   (ed eventually destroyed if the counter reaches 0).
   New target pointer will be held, with counter 1. If the new pointer is nullptr
   the class will be set to empty.
 * \param[in] object Pointer to target object
 *
 */
template<class O>
void MimmoSharedPointer<O>::reset(O* object)
{
    _reset(object, nullptr);
}

/*!
   Dereference operator.
   \return Reference to internal pointed object.
 */
template<class O>
O& MimmoSharedPointer<O>::operator*() const
{
    return *(m_object);
}

/*!
   Dereference operator.
   \return Raw pointer to object.
*/
template<class O>
O* MimmoSharedPointer<O>::operator->() const
{
    return m_object;
};

/*!
   Logical NOT operator.
   \return True if raw internal pointer to object is nullptr.
*/
template<class O>
bool MimmoSharedPointer<O>::operator!() const
{
    return !m_object;
}

/*!
   Boolean conversion operator.
   \return Boolean conversion of raw pointer to object.
*/
template<class O>
MimmoSharedPointer<O>::operator bool() const
{
    return bool(m_object);
}

/*!
   Initialize the shared pointer and the counter.
   \param[in] object Pointer to target object
   \param[in] counter Pointer to counter of sharing pointers. If nullptr a new counter object is created. If object pointer is not nullptr the counter is increased by one.
*/
template<class O>
void MimmoSharedPointer<O>::_init(O* object, std::size_t* counter )
{
    m_object = object;
    m_counter = counter;
    if (m_object){
        if(!m_counter)
            m_counter = new std::size_t(0);
        (*m_counter)++;
    }
    else{
        m_counter = nullptr;
    }
};

/*!
   Reset the class to new contents (raw object pointer and counter pointer). First old contents
   are released using _decrement() function. New contents are then allocated using _init(object, counter).
   \param[in] object Pointer to target object
   \param[in] counter Pointer to counter of sharing pointers.
*/
template<class O>
void MimmoSharedPointer<O>::_reset(O* object, std::size_t* counter )
{
    _decrement(),
    _init(object, counter);
};

/*!
   Method to decrement the internal class counter of 1. If the counter reaches zero,
   the internal pointed object is explicitly deleted.
 */
template<class O>
void MimmoSharedPointer<O>::_decrement(){
    if (m_counter){
        ( *m_counter )-- ;
        if ( (*m_counter) == 0 ){
            delete m_object;
            delete m_counter;
        }
    }
    m_object = nullptr;
    m_counter = nullptr;
}

} //end namespace mimmo

namespace std{

/*!
   Create hash of MimmoSharedPointer using the internal raw object pointer
   \param[in] obj const Reference to MimmoSharedPointer class
 */
template<class O>
size_t hash<mimmo::MimmoSharedPointer<O>>::operator()(const mimmo::MimmoSharedPointer<O> & obj) const
{
    return hash<const O*>()(obj.get());
}

/*!
   Create hash of MimmoSharedPointer using the internal raw object pointer
   \param[in] obj Reference to MimmoSharedPointer class
*/
template<class O>
size_t hash<mimmo::MimmoSharedPointer<O>>::operator()(mimmo::MimmoSharedPointer<O> & obj)
{
    return hash<O*>()(obj.get());
}

} // end namespace std
