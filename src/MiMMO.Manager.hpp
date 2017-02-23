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
 * \*---------------------------------------------------------------------------*/
#ifndef __MIMMOMANAGER_HPP__
#define __MIMMOMANAGER_HPP__

#include "bitpit_common.hpp"
#include <unordered_map>

#define REGISTER_MANIPULATOR(name, type) \
static int manipulator_type = mimmo::registerManipulator(name); 

// \

// BITPIT_UNUSED(manipulator_type);


namespace mimmo{

extern std::unordered_map<std::string, int> _manipulatorList;
	
int registerManipulator(const std::string & name);

};

#endif /*! __MIMMOMANAGER_HPP__*/