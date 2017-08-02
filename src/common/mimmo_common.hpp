/*---------------------------------------------------------------------------
    mimmo
  
    Copyright (C) 2015-2017 OPTIMAD engineering Srl
  
    -------------------------------------------------------------------------
    License
    This file is part of mimmo.
  
    mimmo is free software: you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License v3 (LGPL)
    as published by the Free Software Foundation.
  
    mimmo is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
    License for more details.
  
    You should have received a copy of the GNU Lesser General Public License
    along with mimmo. If not, see <http://www.gnu.org/licenses/>.
  
  ----------------------------------------------------------------------------*/



/*!
 * \defgroup common common
 *  Common classes and utils for mimmo
 * \{
 *      \defgroup common_Utils common_Utils
 *       Collection of general methods and functions related to mimmo
 *       \{
 *       \} 
 * \}
 * \defgroup macro macro
 * Preprocessor macros available in mimmo, for executable class registration and factorization
 * \{
 * \}
 * \defgroup typedefs typedefs
 * Custom typedefs available in mimmo
 * \{
 * \}
 * 
 * \defgroup Ports ports
 * Collection of all port types actually available in mimmo.
 * Each port is declared as a costant string M_\<xxx\>, where M_ prefix identifies the mimmo port.
 * Please refer to each Manipulation Block documentation for port usage examples.
 * \{
 * \}
 * 
 * \defgroup PortContainers port_containers
 * Collection of all possible containers for mimmo ports.
 * Each port container is declared as a costant string MC_\<xxx\>, where MC_ prefix identifies the mimmo port container.
 * \{
 * \}
 * 
 * \defgroup PortData port_data
 * Collection of all possible data type contained in mimmo ports.
 * Each port data type is declared as a costant string MD_\<xxx\>, where MD_ prefix identifies the mimmo port data type.
 * 
 * \{
 * \}
 * 
 */


#ifndef __MIMMO_MODULE_COMMON_HPP__
#define __MIMMO_MODULE_COMMON_HPP__

#include "mimmo_version.hpp"

#include "enum.hpp"
#include "factory.hpp"
#include "portManager.hpp"
#include "mimmoTypeDef.hpp"
#include "TrackingPointer.hpp"
#include "customOperators.hpp"



/*!
 * \brief mimmo main namespace
 * General namespace of the API.
 */
namespace mimmo{
}


#endif
