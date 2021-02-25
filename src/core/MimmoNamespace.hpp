/*---------------------------------------------------------------------------*\
*
*  mimmo
*
*  Copyright (C) 2015-2021 OPTIMAD engineering Srl
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
#ifndef __MIMMONAMESPACE_HPP__
#define __MIMMONAMESPACE_HPP__

#include <mimmo_binary_stream.hpp>
#include <logger.hpp>

namespace mimmo{
/*!
 * \class FileDataInfo
 * \ingroup core
 * \brief FileDataInfo is a struct to stock data relative to names of external files.
 *
 * FileDataInfo has three fields: an integer relative to the type of file stored,
 * a string reporting the path to the file and a string reporting the name of the file
 *
 */
struct FileDataInfo{
    int ftype;          /**< file type identifier*/
    std::string fname;  /**< file name*/
    std::string fdir;   /**< file directory*/

    FileDataInfo();
    virtual ~FileDataInfo();
    FileDataInfo(const FileDataInfo & other);
};

}

/*!
 * \ingroup binaryStream
 * \{
 */
mimmo::IBinaryStream& operator>>(mimmo::IBinaryStream &buf, mimmo::FileDataInfo&  element);
mimmo::OBinaryStream& operator<<(mimmo::OBinaryStream &buf, const mimmo::FileDataInfo& element);
/*!
 *\}
 */


 namespace mimmo{

 class BaseManipulation;

/*!
*
*  \brief Utilities to create port connections between executable blocks.
*  \ingroup core
* Here are collected enums and methods to create connections between blocks throughout established ports.
* All the available port types are registered (or need to be if not) in the API singleton mimmo::PortManager
*/
namespace pin{

/*!
*\enum ConnectionType
*\brief Type of allowed connections of the object: bidirectional, only input or only output.
*/
enum class ConnectionType{
    BOTH 		/**<Bidirectional object. It allows both input and output connections.*/,
    BACKWARD 	/**<Uni-directional backward object. It allows only input connections.*/,
    FORWARD 	/**<Uni-directional forwadr object. It allows only output connections.*/
};

typedef	std::string	PortID; /**< mimmo custom definition */


bool addPin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR, bool forced = false);

void removeAllPins(BaseManipulation* objSend, BaseManipulation* objRec);

void removePin(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);

bool checkCompatibility(BaseManipulation* objSend, BaseManipulation* objRec, PortID portS, PortID portR);


} //end of pin namespace


//log variables
extern std::string MIMMO_LOG_FILE;  /**<Name of logger file.*/
extern std::string MIMMO_LOG_DIR;   /**<Directory of logger file.*/

void    setLogger(std::string log);
void    setLoggerDirectory(std::string dir);

void    warningXML(bitpit::Logger* log, std::string name);

//expert variable
extern bool MIMMO_EXPERT; /**<Flag that defines expert mode (true) or safe mode (false).
                                In case of expert mode active the mandatory ports are not checked. */

void setExpertMode(bool flag = true);

}//end namespace mimmo

#endif
