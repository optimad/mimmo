/*---------------------------------------------------------------------------*\
 * 
 *  mimmo
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
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

#ifndef BLOCK_FACTORY_HPP
#define BLOCK_FACTORY_HPP

#include <unordered_map>
#include <vector>
#include "configuration.hpp"
#include "bitpit_common.hpp"

namespace mimmo{

/*!
 * \class Factory
 * \brief Factory base template singleton for automatic factorization of manipulator classes
 * 
 * Factory register a list of Base classes with their own custom creators and instance them automatically
 */
template <class Base>
class Factory {

private:
	/*!Private constructor*/
	Factory() : defaultCreator(0) {};          
	/*!Prevent copy-construction*/
	Factory(const Factory&);                  
	/*!Prevent assignment*/
	Factory& operator=(const Factory&);
	/*! Destructor */
	~Factory() { deleteCreators(); }; 

public:
	/*!
	 * \class AbstractCreator
	 * \brief Abstract class embedded in Factory to link creators of type Base* <>(const bitpit::Config::Section & xml_root)
	 */
	class AbstractCreator {
		public:
			/*! Destructor */
			virtual ~AbstractCreator() {}
			/*! Pre virtual create method*/
			virtual Base* create(const bitpit::Config::Section & xml_root) const = 0;
	};

	/*! Instance the singleton */
	static Factory& instance()
	{
		static Factory factory;
		return factory;
	}
	/*! Create a Base Class with linked creator */ 
	static Base * create(const std::string name, const bitpit::Config::Section & xml_root)
	{
		Factory<Base>& factory = instance();
		if (factory.creators.count(name) == 0) {
				return 0;
		}
		
		return (factory.creators[name])->create(xml_root);
	}
	
	/*! Remove a creator from registered list */ 
	void removeCreator(const std::string name)
	{
		if (creators.count(name) > 0) {
			delete creators[name];
			creators.erase(name);
		}
		return;
	}
	
	/*! Add a creator of type AbstractCreator to registered list */ 
	int addCreator(const std::string name, const AbstractCreator* creator)
	{
		removeCreator(name);
		creators[name]= creator;
		return (int) creators.size() + 1;
	}

	/*! Check if a class is already registered */
	bool containsCreator(const std::string name)
	{
		return creators.count(name) > 0;
	}

	/*! set the default creator for all your classes */
	int setDefaultCreator(const AbstractCreator* creator)
	{
		defaultCreator = creator;
		return 0;
	}

	/*! Return the whole map of registered classes */
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
 * \class Creator
 * \brief Template class to create an object Base * = new Derived constructor() class,
 * where Derived is a generic derived class of Base.
 */

template <class Base, class Derived>
class Creator : public Factory<Base>::AbstractCreator {

/** Typedef for generic Creator function. All function must have as argument a reference to 
 bitpit::Config::Section object */	
typedef Derived* (*CreateFn) (const bitpit::Config::Section & );

public:
	/*! Constructor */
	Creator(CreateFn fn = NULL) : createFn(fn) {}
	/*! Function to create a new object of class Derived, with function fn linked to Creator */
	virtual Derived* create(const bitpit::Config::Section & xml_root) const {
		return (createFn ? createFn(xml_root) : new Derived(xml_root));
	};
	
	CreateFn createFn; /**< function linked as creator */
};

};

/*!
 *	\brief Macro collection for automatic class registration
 */
#define REGISTER(Base, Derived, name) \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new mimmo::Creator<Base, Derived>());

#define REGISTER_NO_UNUSED(Base, Derived, name) \
static int factory_##Base##_##Derived = Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>()); \
BITPIT_UNUSED(factory_##Base##_##Derived);

#define REGISTER_CUSTOM(Base, Derived, name, customCreator) \
static int factory_##Base##_##Derived = Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>(&customCreator));

#define IS_REGISTERED(Base, name) \
mimmo::Factory<Base>::instance().containsCreator(name);

#endif
