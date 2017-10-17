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

#include "customOperators.hpp"

/*!
 * Custom Operator - for bitpit::Vertex
 * \ingroup common_Utils
 * \param[in] v1 first vertex
 * \param[in] v2 second vertex
 * \return difference vertex coordinates V1-V2;
 */
std::array<double, 3> operator-(const bitpit::Vertex &v1, const bitpit::Vertex &v2){
    std::array<double, 3> coords2 = v1.getCoords() - v2.getCoords();
    return ( coords2 );
}

