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
std::unique_ptr< std::unordered_map<std::string, xmlBuilder> >_manipulatorList(nullptr);

/*!
 * Registering map with class name and their xmlBuilder
 */
int registerManipulator(const std::string & name, xmlBuilder funct)
{
    if (!_manipulatorList) {
		std::unique_ptr< std::unordered_map<std::string, xmlBuilder> > temp(new std::unordered_map<std::string, xmlBuilder>());
        _manipulatorList = std::move(temp);
    }

    int id = (int)_manipulatorList->size();
    _manipulatorList->insert(std::make_pair(name, funct));
    return id;
}

/*!
 * Return the name-xmlBuilder map
 */
const std::unordered_map<std::string, xmlBuilder> & getManipulatorList()
{
    return *(_manipulatorList.get());
}

/*!
 * Return a instantiated Manipulator as unique pointer. Require name of the class and its xml slot where read its parameters
 */
std::unique_ptr<BaseManipulation> factoryManipulator(const std::string & name, const bitpit::Config::Section & xmlslot)
{
	std::unique_ptr<BaseManipulation> temp(nullptr);
	
	if(_manipulatorList->count(name)>0){
		auto fptr = (*(_manipulatorList.get()))[name];
		temp = std::move(fptr(xmlslot));
	}
	return std::move(temp);
};





};
