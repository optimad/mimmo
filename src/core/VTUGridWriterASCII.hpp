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
#ifndef __VTUGRIDASCIIWRITER_HPP__
#define __VTUGRIDASCIIWRITER_HPP__

#include "bitpit_patchkernel.hpp"
#include "mimmoTypeDef.hpp"

namespace mimmo{

/*!
 * \ingroup core
 */

/*!
 * \class VTUAbsorbStreamer
 * \brief Abstract class for custom ASCII writer/flusher of *.vtu mesh external files
 *
 * Writer/flusher is focused to mesh data only retained by PatchKernel classes.
 * Abstract class need to provide:
 *  - A custom implementation of flushData method to write bitpit::PatchKernel data on target file.
 */
class VTUFlushStreamerASCII : public bitpit::VTKBaseStreamer{
public:
    VTUFlushStreamerASCII();
    virtual ~VTUFlushStreamerASCII();
    /*!Copy Constructor */
    VTUFlushStreamerASCII(const VTUFlushStreamerASCII&) = default;

    void setTargetPatch(bitpit::PatchKernel & patch, bitpit::PiercedStorage<long, long> & vtkVertexMap);
    virtual void flushData(std::fstream &stream, const std::string & name, bitpit::VTKFormat format = bitpit::VTKFormat::ASCII);

private:
    bitpit::PatchKernel * m_patch;
    bitpit::PiercedStorage<long,long> * m_vtkVertexMap;
};

/*!
 * \class VTUGridWriterASCII
 * \brief Custom writer of ASCII unstructured grids to external files *.vtu
 *
 * Writer of bitpit::PatchKernel 's as unstructured grids to external files *.vtu in ASCII format only.
 */
class VTUGridWriterASCII: protected bitpit::VTKUnstructuredGrid
{

public:
    VTUGridWriterASCII(VTUFlushStreamerASCII & streamer,bitpit::PatchKernel & patch,
                        bitpit::VTKElementType eltype= bitpit::VTKElementType::UNDEFINED);
    ~VTUGridWriterASCII();

    void write(const std::string & dir, const std::string & file, bitpit::VTKWriteMode mode= bitpit::VTKWriteMode::DEFAULT);

private:
    bitpit::PatchKernel& m_patch;   /**< reference to patch kernel data structure to fill*/
    VTUFlushStreamerASCII & m_streamer; /**< reference to streamer which knows how to write data to file in ascii */
    bitpit::PiercedStorage<long,long> m_vtkVertexMap;
};


};

#endif /* __VTUGRIDASCIIWRITER_HPP__ */
