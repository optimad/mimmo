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

#include "MiMMO.Manager.hpp"

using namespace mimmo;

namespace mimmo {

// Initial instantiation of the manipulators' list
std::unique_ptr< std::unordered_map<std::string, int> >_manipulatorList(nullptr);

/*!
 * Registering a manipulator class with a target name 
 */
int registerManipulator(const std::string & name)
{
    if (!_manipulatorList) {
		std::unique_ptr< std::unordered_map<std::string, int> > temp(new std::unordered_map<std::string, int>());
        _manipulatorList = std::move(temp);
    }

    int id = (int)_manipulatorList->size();
    _manipulatorList->insert(std::make_pair(name, id));
    return id;
}

/*!
 * Return the map of available objects
 */
const std::unordered_map<std::string, int> & getManipulatorList()
{
    return *(_manipulatorList.get());
}

};
