/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "Info.hpp"

///*!Default constructor of Info
// */
Info::Info(){};
//
///*!Default destructor of Info
// */
Info::~Info(){};

///*!Copy constructor of Info.
// */
Info::Info(const Info & other){
	m_naxes 	= other.m_naxes;
	m_axes 		= other.m_axes;
	m_origin 	= other.m_origin;
	m_npoints 	= other.m_npoints;
	m_coords 	= other.m_coords;
};

/*!Assignement operator of Info.
 */
Info & Info::operator=(const Info & other){
	m_naxes 	= other.m_naxes;
	m_axes 		= other.m_axes;
	m_origin 	= other.m_origin;
	m_npoints 	= other.m_npoints;
	m_coords 	= other.m_coords;
	return (*this);
};
