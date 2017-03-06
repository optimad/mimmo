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
 \ *---------------------------------------------------------------------------*/

#ifndef BLOCK_FACTORY_HPP
#define BLOCK_FACTORY_HPP

#include <unordered_map>
#include <vector>
#include "configuration.hpp"
#include "bitpit_common.hpp"

/*!
 *	\date			02/March/2017
 *	\authors		Andrea Iob
 *	\authors		Rocco Arpa
 *
 *	\brief Factory base template singleton for automatic factorization of manipulator classes
 */
namespace mimmo{
	
template <class Base>
class Factory {

private:
	Factory() : defaultCreator(0) {}         // Private constructor
	Factory(const Factory&);                 // Prevent copy-construction
	Factory& operator=(const Factory&);      // Prevent assignment
	~Factory() { deleteCreators(); } 

public:
	class AbstractCreator {
		public:
			virtual ~AbstractCreator() {}
			virtual Base* create(const bitpit::Config::Section & xml_root) const = 0;
	};

	static Factory& instance()
	{
		static Factory factory;
		return factory;
	}

	static Base * create(const std::string name, const bitpit::Config::Section & xml_root)
	{
		Factory<Base>& factory = instance();
		if (factory.creators.count(name) == 0) {
				return 0;
		}
		
		return (factory.creators[name])->create(xml_root);
	}
	
	void removeCreator(const std::string name)
	{
		if (creators.count(name) > 0) {
			delete creators[name];
			creators.erase(name);
		}
		return;
	}

	int addCreator(const std::string name, const AbstractCreator* creator)
	{
		removeCreator(name);
		creators[name]= creator;
		return (int) creators.size() + 1;
	}

	bool containsCreator(const std::string name)
	{
		return creators.count(name) > 0;
	}

	int setDefaultCreator(const AbstractCreator* creator)
	{
		defaultCreator = creator;
		return 0;
	}

	std::vector<std::string> mapRegisteredBlocks(){
		std::vector<std::string> result(creators.size());
		int counter = 0;
		for(auto &val : creators){
			result[counter] = val.first;
			++counter;
		}
		return result;
	}
private:
	const AbstractCreator* defaultCreator;
	std::unordered_map<std::string, const AbstractCreator*> creators;
	
	void deleteCreators(){
		for(auto & val : creators){
		    delete val.second;
		}
		creators.clear();
	};
	
};

/*!
 *	\date			02/March/2017
 *	\authors		Andrea Iob
 *
 *	\brief Template class to create a Base * = new Derived constructor() class.
 */
template <class Base, class Derived>
class Creator : public Factory<Base>::AbstractCreator {

typedef Derived* (*CreateFn) (const bitpit::Config::Section & );

public:
	Creator(CreateFn fn = NULL) : createFn(fn) {}
	virtual Derived* create(const bitpit::Config::Section & xml_root) const {
		return (createFn ? createFn(xml_root) : new Derived(xml_root));
	};

	CreateFn createFn;
};

};
/*!
 *	\date			02/March/2017
 *	\authors		Andrea Iob
 *
 *	\brief Macro collection for automatic class registration
 */
// #define REGISTER_DEFAULT(Base) \
// static int factory_##Base = Factory<Base>::instance().setDefaultCreator(new Creator<Base, Base>());

// #define REGISTER_DEFAULT_CUSTOM(Base, customCreator) \
// static int factory_##Base = Factory<Base>::instance().setDefaultCreator(new Creator<Base, Base>(&customCreator));

#define REGISTER(Base, Derived, name) \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new mimmo::Creator<Base, Derived>());

#define REGISTER_MANIP(Base, Derived, name) \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new mimmo::Creator<Base, Derived>());

// #define REGISTER_NO_UNUSED(Base, Derived, name) \
// static int factory_##Base##_##Derived = Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>()); \
// BITPIT_UNUSED(factory_##Base##_##Derived);
// 
// #define REGISTER_CUSTOM(Base, Derived, name, customCreator) \
// static int factory_##Base##_##Derived = Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>(&customCreator));

#define IS_REGISTERED(Base, name) \
mimmo::Factory<Base>::instance().containsCreator(name);

#endif
