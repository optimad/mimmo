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

#ifndef BLOCK_FACTORY_HPP
#define BLOCK_FACTORY_HPP

#include <unordered_map>
#include <vector>
#include "configuration.hpp"
#include "bitpit_common.hpp"

namespace mimmo{

/*!
 *   \ingroup common
 *   \{
 */    
    
/*!
 * \class Factory
 * \brief Factory base template singleton for automatic 
 *  factorization of executable classes
 *
 * Factory register a list of Base classes with their own custom 
 * method which create and instance them automatically.
 */
template <class Base>
class Factory {

private:
    /*!Private constructor of the singleton*/
    Factory() : defaultCreator(0) {};
    /*!Prevent copy-construction for singleton*/
    Factory(const Factory&);
    /*!Prevent assignment for singleton*/
    Factory& operator=(const Factory&);
    /*! Destructor */
    ~Factory() { 
        deleteCreators(); 
    }; 

public:
    /*!
    * \class AbstractCreator
    * \brief Abstract class embedded in Factory to link creators 
    * of type Base* <>(const bitpit::Config::Section & xml_root)
    */
    class AbstractCreator {
        public:
            /*! Destructor */
            virtual ~AbstractCreator() {}
            /*! Pre virtual create method
             * \param[in]   xml_root reference to xml-data tree.
             */
            virtual Base* create(const bitpit::Config::Section & xml_root) const = 0;
    };

    /*! Instance the singleton */
    static Factory& instance(){
        static Factory factory;
        return factory;
    }
    /*! Create a Base Class with linked creator
     * \param[in] name name of a creator registered in the singleton
     * \param[in] xml_root reference to info, passed as xml-data tree, referred to the creator
     * \return pointer to the newly instantiated base class
     */ 
    static Base * create(const std::string name, const bitpit::Config::Section & xml_root){
        Factory<Base>& factory = instance();
        if (factory.creators.count(name) == 0) return 0;
        return (factory.creators[name])->create(xml_root);
    }
    
    /*! Remove a creator from registered list
     * \param[in] name name of the creator
     */ 
    void removeCreator(const std::string name){
        if (creators.count(name) > 0) {
            delete creators[name];
            creators.erase(name);
        }
        return;
    }
    
    /*! Add a creator of type AbstractCreator to registered list
     * \param[in] name name of the creator
     * \param[in] creator pointer to the creator
     * \return current size of creators registered in the list
     */ 
    int addCreator(const std::string name, const AbstractCreator* creator){
        removeCreator(name);
        creators[name]= creator;
        return (int) creators.size() + 1;
    }

    /*! Check if a class is already registered 
     * \param[in] name name of the creator
     * \return true if the creator exists
     */
    bool containsCreator(const std::string name){
        return creators.count(name) > 0;
    }

    /*! It sets a default creator for all your possible 
     * registered classes.
     * \param[in] creator pointer to a creator
     */
    int setDefaultCreator(const AbstractCreator* creator){
        defaultCreator = creator;
        return 0;
    }

    /*!\return the whole list of registered creators
     */
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
    const AbstractCreator* defaultCreator; /**< pointer to a default creator. If a registered class has none, use this*/
    std::unordered_map<std::string, const AbstractCreator*> creators; /**< list of registered creators */

    /*!
     * Empty the whole current list of registered creators
     */ 
    void deleteCreators(){
        for(auto & val : creators){
            delete val.second;
        }
        creators.clear();
    };
};

/*!
* \class Creator
* \brief Template class to create an object Base * = new Derived creator class,
* where Derived is a generic derived class of Base. Creator can be any constructor 
* of the derived class who takes as argument a const reference to a xml-data tree 
* bitpit::Config::Section.
*/
template <class Base, class Derived>
class Creator : public Factory<Base>::AbstractCreator {

/*! 
 * \typedef CreateFn
 * for generic Creator function. 
 * All function must have as argument a reference to bitpit::Config::Section object 
 */
typedef Derived* (*CreateFn) (const bitpit::Config::Section & );

public:
    CreateFn createFn; /**< function linked as Creator */
    
    /*! Constructor
     * \param[in] fn CreateFn function
     */
    Creator(CreateFn fn = NULL) : createFn(fn) {}
    /*! Function to create a new object of class Derived, 
     * with custom function of type CreateFn linked as Creator
     * \param[in] xml_root reference of xml-data tree
     * \return pointer to newly created object
     */
    virtual Derived* create(const bitpit::Config::Section & xml_root) const {
        return (createFn ? createFn(xml_root) : new Derived(xml_root));
    };
};

};

/*!
 * \}
 */

/*!
 * \ingroup macro
 * \{
 */

/*!
 * \def REGISTER_NO_UNUSED(Base, Derived, name)
 * Register an executable mimmo class in the Factory singleton, use unused bipit option 
 * internally on the counter of the already registered class to avoid compiler warnings.
 * \param[in]   Base name of Derived's base class
 * \param[in]   Derived name of the class
 * \param[in]   name string to register the creator's name
 */
#define REGISTER_NO_UNUSED(Base, Derived, name) \
/* register a Derived class with its xml default constructor that will be instantiate as Base. Return is unused*/ \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>()); \
BITPIT_UNUSED(factory_##Base##_##Derived);


/*!
 * \def REGISTER(Base, Derived, name)
 *  Register an executable mimmo class in the Factory singleton.
 * \param[in]   Base name of Derived's base class
 * \param[in]   Derived name of the class
 * \param[in]   name string to register the creator's name
 * \return counter of the already registered classes
 */
#define REGISTER(Base, Derived, name) \
/* register a Derived class with its xml default constructor that will be instantiate as Base*/ \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new mimmo::Creator<Base, Derived>());

/*!
 * \def REGISTER_CUSTOM(Base, Derived, name, customCreator)
 * Register an executable mimmo class in the Factory singleton, and passing with it the method
 * also a custom creator method
 * \param[in]   Base name of Derived's base class
 * \param[in]   Derived name of the class
 * \param[in]   name string to register the creator's name
 * \param[in]   customCreator method to create the class
 * \return counter of the already registered classes
 * 
 */
#define REGISTER_CUSTOM(Base, Derived, name, customCreator) \
/* register a Derived class a custom xml constructor/creator method that will be instantiate as Base*/ \
static int factory_##Base##_##Derived = mimmo::Factory<Base>::instance().addCreator(name, new Creator<Base, Derived>(&customCreator));

/*!
 * \def IS_REGISTERED(Base, name)
 * Return if a creator name referencing to a Base class is already registered or not.
 * \param[in]   Base name of Derived's base class
 * \param[in]   name string to register the creator's name
 * \return boolean true/false if the creator name is already registered
 */
#define IS_REGISTERED(Base, name) \
/* find if a class creator is already registered*/ \
mimmo::Factory<Base>::instance().containsCreator(name);


/*!
 * \}
 */ 

#endif

