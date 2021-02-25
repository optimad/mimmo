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
#ifndef __VTUGRIDREADER_HPP__
#define __VTUGRIDREADER_HPP__

#include "mimmoTypeDef.hpp"
#include <bitpit_patchkernel.hpp>

namespace mimmo{

/*!
 * \class VTUAbsorbStreamer
 * \brief Abstract class for custom reader/absorber of *.vtu mesh external files
 * \ingroup core

 * Reader/absorber is focused to mesh data only.
 * Abstract class need to provide:
 *  - A custom implementation of absorbData method to read the data from file (actually void).
 *  - A method to decode raw data acquired by the streamer, in order to fill a
 * custom mesh container of type bitpit::PatchKernel.
 */
class VTUAbsorbStreamer : public bitpit::VTKBaseStreamer{
public:
    VTUAbsorbStreamer();
    virtual ~VTUAbsorbStreamer();
    /*!Copy Constructor */
    VTUAbsorbStreamer(const VTUAbsorbStreamer&) = default;

    virtual void absorbData(std::fstream &stream, const std::string & name, bitpit::VTKFormat format, uint64_t entries, uint8_t components, bitpit::VTKDataType datatype);
    /*! Decode read raw data and fill a bitpit::PatchKernel structure with them */
    virtual void decodeRawData(bitpit::PatchKernel &) = 0;
};

/*!
 * \class VTUGridStreamer
 * \brief Custom mesh/data absorber for unstructured grids given by external files *.vtu
   \ingroup core
 *
 * Read unstructured mesh data only from an external file vtu (any other data will be skipped).
 * Data will be stored in internal structures of the streamer.
 */
class VTUGridStreamer: public VTUAbsorbStreamer{

protected:
    dvecarr3E points;                            /**< list of mesh nodes */
    livector1D connectivitylist;                 /**< one dimensional list of all vertex indices composing connectivity */
    livector1D pointsID;                         /**< labels of mesh nodes, if any*/
    livector1D cellsID;                          /**< labels of mesh cells, if any*/
    livector1D pids;                              /**< cell pid, if any */
    livector1D offsets;                           /**< offset list for decoding connectivity list */
    std::vector<bitpit::ElementType> types;      /**< list of cell types */
    livector1D faces;                             /**< list of face vertex connectivity - for 3D volume polyhedral support */
    livector1D faceoffsets;                       /**< list of offsets to read face vertex connectivity - for 3D volume polyhedral support */
#if MIMMO_ENABLE_MPI
    livector1D cellsRank;                         /**< ranks of mesh cells, if any*/
    livector1D verticesRank;                       /**< ranks of mesh vertices, if any*/
    livector1D cellsGlobalIndex;                   /**< global index of mesh cells, if any*/
#endif

public:
    VTUGridStreamer();
    virtual ~VTUGridStreamer();
    /*! Copy Constructor*/
    VTUGridStreamer(const VTUGridStreamer&) = default;

    virtual void absorbData(std::fstream &stream, const std::string &name, bitpit::VTKFormat format, uint64_t entries, uint8_t components, bitpit::VTKDataType datatype);
    void decodeRawData(bitpit::PatchKernel & patch);
    void readIntegerPod(std::fstream & stream, bitpit::VTKFormat& format, bitpit::VTKDataType& datatype, long &target);
    void readDoublePod(std::fstream & stream, bitpit::VTKFormat& format, bitpit::VTKDataType& datatype, double&target);
};

/*!
 * \class VTUGridReader
 * \brief Custom reader of unstructured grids from external files *.vtu
 * \ingroup core

 * Reader of unstructured grids from external files *.vtu. if successfull reading,
 * store the mesh fields in target bitpit::Patchkernel data structure.
 * Need in construction to specify a streamer of type VTUAbsorbStreamer.
 */
class VTUGridReader: protected bitpit::VTKUnstructuredGrid
{

public:
    VTUGridReader( std::string dir, std::string name, VTUAbsorbStreamer & streamer,
                   bitpit::PatchKernel & patch, bool masterRankOnly = true, bitpit::VTKElementType eltype= bitpit::VTKElementType::UNDEFINED);
    ~VTUGridReader();

    void read();

private:
    bitpit::PatchKernel& m_patch;   /**< reference to patch kernel data structure to fill*/
    VTUAbsorbStreamer & m_streamer; /**< reference to streamer which knows how to read data from file */
#if MIMMO_ENABLE_MPI
    bool m_masterRankOnly; /**< control reading on master rank only **/
#endif
};


};

#endif /* __VTUGRIDREADER_HPP__ */
