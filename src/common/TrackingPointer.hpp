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

#ifndef __TRACKINGPOINTER_HPP__
#define __TRACKINGPOINTER_HPP__

namespace mimmo{

    /*!
     * \ingroup common
     * \{
     */

/*!
 * \class TrackingPointer
 * \brief Basic virtual class to derive a generic object whose pointer can return
 * an identifying name through the method whichClass.
 */
class TrackingPointer{
public:
    TrackingPointer(){};
    virtual ~TrackingPointer(){};

    /*!
     * Pure virtual class to recover the identifying name
     */
    virtual std::string whichClass() = 0;
};

/*!
 * \}
 */

}

#endif /* __TRACKINGPOINTER_HPP__ */
