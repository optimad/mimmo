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
\*---------------------------------------------------------------------------*/

#include "VTUGridReader.hpp"

namespace mimmo{

/*!
 * Base Constructor
 */
VTUAbsorbStreamer::VTUAbsorbStreamer() : VTKBaseStreamer(){}

/*!
 * Base Destructor
 */
VTUAbsorbStreamer::~VTUAbsorbStreamer(){}

/*!
 * Absorber of VTU mesh data. Reimplemented from bitpit::VTKBaseStreamer class
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */
void VTUAbsorbStreamer::absorbData(std::fstream &stream, const std::string &name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    bitpit::VTKBaseStreamer::absorbData(stream, name, format, entries, components, datatype);
}
/*!
 * Absorber to get BITPIT_LEGACY from versions lesser then 1.6 release.
 * Absorber of VTU mesh data. Reimplemented from bitpit::VTKBaseStreamer class
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */
void VTUAbsorbStreamer::absorbData(std::fstream &stream, std::string name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    bitpit::VTKBaseStreamer::absorbData(stream, name, format, entries, components, datatype);
}



/*!
 * Base Constructor
 */
VTUGridStreamer::VTUGridStreamer():VTUAbsorbStreamer(){}

/*!
 * Base Destructor
 */
VTUGridStreamer::~VTUGridStreamer(){}

/*!
 * Absorber of VTU mesh data. Reimplemented from bitpit::VTKBaseStreamer class
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */
void VTUGridStreamer::absorbData(std::fstream &stream, const std::string &name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    std::size_t sizeData = std::size_t(entries/components);
    if (name == "Points") {
        //get correct data type
        points.resize(sizeData);
        for (auto & p :points) {
            for(auto &pval : p){
                switch (datatype){
                    case bitpit::VTKDataType::Float32 :
                    {
                        float val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    case bitpit::VTKDataType::Float64:
                    {
                        double val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    default:
                        throw std::runtime_error("VTUGridStreamer::absorbData : Points data format not available");
                    break;
                }

            }
        }
    } else if (name == "offsets") {
        offsets.resize(sizeData);
        for (auto & voffset : offsets) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Offsets data format not available");
                break;
            }
        }
    } else if (name == "types") {
        types.resize(sizeData);
        for (auto & vtype : types) {

            int dumType;
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Types data format not available");
                break;
            }
            //find ElementType
            switch (dumType)  {
                case 1:
                    vtype = bitpit::ElementType::VERTEX;
                    break;

                case 3:
                    vtype = bitpit::ElementType::LINE;
                    break;

                case 5:
                    vtype = bitpit::ElementType::TRIANGLE;
                    break;

                case 7:
                    vtype = bitpit::ElementType::POLYGON;
                    break;

                case 8:
                    vtype = bitpit::ElementType::PIXEL;
                    break;

                case 9:
                    vtype = bitpit::ElementType::QUAD;
                    break;

                case 10:
                    vtype = bitpit::ElementType::TETRA;
                    break;

                case 11:
                    vtype = bitpit::ElementType::VOXEL;
                    break;

                case 12:
                    vtype = bitpit::ElementType::HEXAHEDRON;
                    break;

                case 13:
                    vtype = bitpit::ElementType::WEDGE;
                    break;

                case 14:
                    vtype = bitpit::ElementType::PYRAMID;
                    break;

                case 42:
                    vtype = bitpit::ElementType::POLYHEDRON;
                    break;

                default:
                    vtype = bitpit::ElementType::UNDEFINED;
                    break;

            }
        }
    } else if (name == "connectivity") {
        connectivitylist.resize(sizeData);
        for (auto & vconn : connectivitylist) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Connectivity data format not available");
                    break;
            }
        }
    } else if (name == "faces") {
        faces.resize(sizeData);
        for (auto & vface : faces) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Faces data format not available");
                    break;
            }
        }
    } else if (name == "faceoffsets") {
        faceoffsets.resize(sizeData);
        for (auto & vfoffset : faceoffsets) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Faceoffsets data format not available");
                    break;
            }
        }
    } else if (name == "cellIndex") {
        cellsID.resize(sizeData);
        for (auto & vcid : cellsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : CellIndex data format not available");
                    break;
            }
        }
    } else if (name == "PID") {
        pids.resize(sizeData);
        for (auto & pid : pids) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : PID data format not available");
                    break;
            }
        }
    } else if (name == "vertexIndex") {
        pointsID.resize(sizeData);
        for (auto & vpid : pointsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : VertexIndex data format not available");
                    break;
            }
        }
    }
}

/*!
 * Absorber to get BITPIT_LEGACY from versions lesser then 1.6 release.
 * Absorber of VTU mesh data. Reimplemented from bitpit::VTKBaseStreamer class
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */
void VTUGridStreamer::absorbData(std::fstream &stream, std::string name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    std::size_t sizeData = std::size_t(entries/components);
    if (name == "Points") {
        //get correct data type
        points.resize(sizeData);
        for (auto & p :points) {
            for(auto &pval : p){
                switch (datatype){
                    case bitpit::VTKDataType::Float32 :
                    {
                        float val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    case bitpit::VTKDataType::Float64:
                    {
                        double val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    default:
                        throw std::runtime_error("VTUGridStreamer::absorbData : Points data format not available");
                    break;
                }

            }
        }
    } else if (name == "offsets") {
        offsets.resize(sizeData);
        for (auto & voffset : offsets) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    voffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Offsets data format not available");
                break;
            }
        }
    } else if (name == "types") {
        types.resize(sizeData);
        for (auto & vtype : types) {

            int dumType;
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    dumType = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Types data format not available");
                break;
            }
            //find ElementType
            switch (dumType)  {
                case 1:
                    vtype = bitpit::ElementType::VERTEX;
                    break;

                case 3:
                    vtype = bitpit::ElementType::LINE;
                    break;

                case 5:
                    vtype = bitpit::ElementType::TRIANGLE;
                    break;

                case 7:
                    vtype = bitpit::ElementType::POLYGON;
                    break;

                case 8:
                    vtype = bitpit::ElementType::PIXEL;
                    break;

                case 9:
                    vtype = bitpit::ElementType::QUAD;
                    break;

                case 10:
                    vtype = bitpit::ElementType::TETRA;
                    break;

                case 11:
                    vtype = bitpit::ElementType::VOXEL;
                    break;

                case 12:
                    vtype = bitpit::ElementType::HEXAHEDRON;
                    break;

                case 13:
                    vtype = bitpit::ElementType::WEDGE;
                    break;

                case 14:
                    vtype = bitpit::ElementType::PYRAMID;
                    break;

                case 42:
                    vtype = bitpit::ElementType::POLYHEDRON;
                    break;

                default:
                    vtype = bitpit::ElementType::UNDEFINED;
                    break;

            }
        }
    } else if (name == "connectivity") {
        connectivitylist.resize(sizeData);
        for (auto & vconn : connectivitylist) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vconn = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Connectivity data format not available");
                    break;
            }
        }
    } else if (name == "faces") {
        faces.resize(sizeData);
        for (auto & vface : faces) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vface = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Faces data format not available");
                    break;
            }
        }
    } else if (name == "faceoffsets") {
        faceoffsets.resize(sizeData);
        for (auto & vfoffset : faceoffsets) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vfoffset = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : Faceoffsets data format not available");
                    break;
            }
        }
    } else if (name == "cellIndex") {
        cellsID.resize(sizeData);
        for (auto & vcid : cellsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vcid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : CellIndex data format not available");
                    break;
            }
        }
    } else if (name == "PID") {
        pids.resize(sizeData);
        for (auto & pid : pids) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    pid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : PID data format not available");
                    break;
            }
        }
    } else if (name == "vertexIndex") {
        pointsID.resize(sizeData);
        for (auto & vpid : pointsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUGridStreamer::absorbData : VertexIndex data format not available");
                    break;
            }
        }
    }
}

/*!
 * Use streamer absorbed data, if any, to fill vertices and cells info of target bitpit::PatchKernel container,
 * passed externally. Please notice, container must be empty.
 * \param[in] patch    container for mesh
 */
void VTUGridStreamer::decodeRawData(bitpit::PatchKernel & patch)
{

    //time to check reading result.
    std::size_t nVertices, nCells;
    nVertices = points.size();
    nCells = offsets.size();
    if(nVertices == 0){
        throw std::runtime_error("Error VTUGridStreamer : no point coordinates detected while reading *.vtu file.");
    }
    if(nCells == 0 || types.empty() || connectivitylist.size() < nCells ){
        throw std::runtime_error("Error VTUGridStreamer : no valid connectivity/offsets/types info detected while reading *.vtu file.");
    }

    patch.reserveVertices(nVertices);
    patch.reserveCells(nCells);

    //reading mesh nodes and store it in vertices.
    //check labels if any;
    bool checkPointsID;
    {
        std::unordered_set<long> checkSet(pointsID.begin(), pointsID.end());
        checkPointsID = ( checkSet.size() == nVertices);
    }
    //insert points and recover local/global map of vertices;
    std::unordered_map<int, long> mapVert;
    int counter = 0;
    long idV=0;
    for(const auto & p : points){
        idV = bitpit::Vertex::NULL_ID;
        if(checkPointsID) idV = pointsID[counter];
        bitpit::PatchKernel::VertexIterator it = patch.addVertex(p, idV);
        mapVert[counter] = (*it).getId();
        ++counter;
    }

    //reading mesh connectivity by offsets and store it in cells.
    //check cell labels if any;
    bool checkCellsID, checkFaceOffset, checkPID;
    {
        std::unordered_set<long> checkSet(cellsID.begin(), cellsID.end());
        checkCellsID = ( checkSet.size() == nCells);
        checkFaceOffset = ( (faceoffsets.size() == nCells) && (faces.size()>= nCells) );
        checkPID = (pids.size() == nCells);
    }
    //insert points;
    counter = 0;
    long idC=0;
    int posCellBegin = 0, posFaceBegin=0;
    bitpit::ElementType eltype;
    long PID;
    livector1D conn;
    std::size_t connSize;
    for(const auto & off : offsets){
        idC = bitpit::Cell::NULL_ID;
        PID = 0;
        if(checkCellsID)  {idC= cellsID[counter];}
        if(checkPID)      {PID = pids[counter];}
        eltype = types[counter];

        if(eltype == bitpit::ElementType::POLYHEDRON){
            if(!checkFaceOffset){
                throw std::runtime_error("Error VTUGridStreamer : trying to acquire POLYHEDRON info without faces and faceoffsets data");
            }
            connSize = faceoffsets[counter] - posFaceBegin;
            conn.resize(connSize);
            int loc = 0;
            for(int i=posFaceBegin; i<faceoffsets[counter]; ++i){
                conn[loc] = faces[i];
                ++loc;
            }
            //remap vertices: conn is now written face by face with local vertex indices
            int posfbegin = 1, posfend; //begin from 1- value. 0 value of conn contains the total number fo faces
            while(posfbegin < int(connSize)){
                posfend = posfbegin +conn[posfbegin] + 1;
                for(int i=posfbegin+1; i<posfend; ++i){
                    conn[i] = mapVert[conn[i]];
                }
                posfbegin = posfend;
            }
        }else if(eltype == bitpit::ElementType::POLYGON){
            connSize = off - posCellBegin;
            conn.resize(connSize +1);
            conn[0] = connSize;
            int loc =1;
            for(int i=posCellBegin; i<off; ++i){
                conn[loc] = mapVert[connectivitylist[i]];
                ++loc;
            }
        }else{
            connSize = off - posCellBegin;
            conn.resize(connSize);
            int loc =0;
            for(int i=posCellBegin; i<off; ++i){
                conn[loc] = mapVert[connectivitylist[i]];
                ++loc;
            }
        }

        bitpit::PatchKernel::CellIterator it = patch.addCell(eltype, conn, idC);
        (*it).setPID(PID);

        posCellBegin = off;
        if(checkFaceOffset){
            if(faceoffsets[counter] > 0)    posFaceBegin = faceoffsets[counter];
        }
        ++counter;
    }
}

/*!
 * Base Constructor
 */
VTUPointCloudStreamer::VTUPointCloudStreamer():VTUAbsorbStreamer(){}

/*!
 * Base Destructor
 */
VTUPointCloudStreamer::~VTUPointCloudStreamer(){}

/*!
 * Absorber of VTU mesh data. Reimplemented from bitpit::VTKBaseStreamer class
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for the data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */
void VTUPointCloudStreamer::absorbData(std::fstream &stream, const std::string &name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    std::size_t sizeData = std::size_t(entries/components);
    if (name == "Points") {
        //get correct data type
        points.resize(sizeData);
        for (auto & p :points) {
            for(auto &pval : p){
                switch (datatype){
                    case bitpit::VTKDataType::Float32 :
                    {
                        float val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    case bitpit::VTKDataType::Float64:
                    {
                        double val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    default:
                        throw std::runtime_error("VTUGridStreamer::absorbData : Points data format not available");
                        break;
                }

            }
        }
    }else if (name == "vertexIndex") {
        pointsID.resize(sizeData);
        for (auto & vpid : pointsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUPointCloudStreamer::absorbData : VertexIndex data format not available");
                    break;
            }
        }
    }
}
/*!
 * method to ensure bitpit_legacy with versions lesser then 1.6 release
 * \param[in] stream    stream to read from
 * \param[in] name      name of the geometry field
 * \param[in] format    ASCII or APPENDED
 * \param[in] entries   number of entries for the data container
 * \param[in] components number of components of current data container
 * \param[in] datatype   data format for binary casting
 */

void VTUPointCloudStreamer::absorbData(std::fstream &stream, std::string name, bitpit::VTKFormat format,
                                 uint64_t entries, uint8_t components, bitpit::VTKDataType datatype)
{
    std::size_t sizeData = std::size_t(entries/components);
    if (name == "Points") {
        //get correct data type
        points.resize(sizeData);
        for (auto & p :points) {
            for(auto &pval : p){
                switch (datatype){
                    case bitpit::VTKDataType::Float32 :
                    {
                        float val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    case bitpit::VTKDataType::Float64:
                    {
                        double val;
                        if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                        else                                    bitpit::genericIO::absorbBINARY(stream, val);
                        pval = val;
                    }
                    break;
                    default:
                        throw std::runtime_error("VTUGridStreamer::absorbData : Points data format not available");
                        break;
                }

            }
        }
    }else if (name == "vertexIndex") {
        pointsID.resize(sizeData);
        for (auto & vpid : pointsID) {
            switch(datatype){
                case bitpit::VTKDataType::Int8 :
                {
                    int8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt8 :
                {
                    uint8_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int16 :
                {
                    int16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt16 :
                {
                    uint16_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int32 :
                {
                    int32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt32 :
                {
                    uint32_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::Int64 :
                {
                    int64_t val;
                    if(format == bitpit::VTKFormat::ASCII)  bitpit::genericIO::absorbASCII(stream,val);
                    else                                    bitpit::genericIO::absorbBINARY(stream, val);
                    vpid = val;
                }
                break;
                case bitpit::VTKDataType::UInt64 :
                default:
                    throw std::runtime_error("VTUPointCloudStreamer::absorbData : VertexIndex data format not available");
                    break;
            }
        }
    }
}


/*!
 * Use streamer absorbed data, if any, to fill vertices of target bitpit::PatchKernel container,
 * passed externally. Please notice, container must be empty.
 * \param[in] patch    container for mesh
 */
void VTUPointCloudStreamer::decodeRawData(bitpit::PatchKernel & patch)
{
    //time to check reading result.
    std::size_t nVertices;
    nVertices = points.size();
    if(nVertices == 0){
        throw std::runtime_error("Error VTUPointCloudStreamer : no point coordinates detected while reading *.vtu file.");
    }
    patch.reserveVertices(nVertices);

    //reading mesh nodes and store it in vertices.
    //check labels if any;
    bool checkPointsID;
    {
        std::unordered_set<long> checkSet(pointsID.begin(), pointsID.end());
        checkPointsID = ( checkSet.size() == nVertices);
    }
    //insert points and recover local/global map of vertices;
    std::unordered_map<int, long> mapVert;
    int counter = 0;
    long idV=0;
    for(const auto & p : points){
        idV = bitpit::Vertex::NULL_ID;
        if(checkPointsID) idV = pointsID[counter];
        bitpit::PatchKernel::VertexIterator it = patch.addVertex(p, idV);
        mapVert[counter] = (*it).getId();
        ++counter;
    }
}



/*!
 * Base constructor. Linked reference bitpit::PatchKernel container must be empty. If not,
 * class will provide to destroy its previous contents and fill it with new read values.
 * \param[in] dir   target directory of file to be read
 * \param[in] name  name fo the file to be read.
 * \param[in] patch reference to empty container for storing mesh data.
 * \param[in] streamer streaming class to absorb VTK data.
 * \param[in] eltype [optional] force the elementtype of the grid.
 */
VTUGridReader::VTUGridReader( std::string dir, std::string name, VTUAbsorbStreamer & streamer, bitpit::PatchKernel & patch, bitpit::VTKElementType eltype) :
                              VTKUnstructuredGrid(dir, name, eltype), m_patch(patch), m_streamer(streamer)
{
    for(auto & field : m_geometry){
        field.enable();
        field.setStreamer(streamer);
    }

    addData<long>("vertexIndex", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, &streamer);
    addData<long>("cellIndex", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, &streamer);
    addData<int>("PID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, &streamer);
}

/*!
 * Basic Destructor
 */
VTUGridReader::~VTUGridReader(){}

/*!
 * Read the file. Reimplemented from bitpit::VTKUnstructuredGrid::read().
 */
void VTUGridReader::read(){

    VTKUnstructuredGrid::read();
     //clear target data
    m_patch.reset();
    m_streamer.decodeRawData(m_patch);
}

}
