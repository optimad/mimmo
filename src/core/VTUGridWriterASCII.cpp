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

#include "VTUGridWriterASCII.hpp"

namespace mimmo{

/*!
 * Base Constructor
 * @param[in] patch target patchkernel
 */
VTUFlushStreamerASCII::VTUFlushStreamerASCII() : VTKBaseStreamer(){
    m_patch = nullptr;
}

/*!
 * Base Destructor
 */
VTUFlushStreamerASCII::~VTUFlushStreamerASCII(){}

/*!
 * Link the target patch having data to be written
 * @param[in] patch target patchkernel
 * @param[in] vtkVertexMap map of vtk vertexindex to be plot.
 */
void VTUFlushStreamerASCII::setTargetPatch(bitpit::PatchKernel & patch, bitpit::PiercedStorage<long, long> & vtkVertexMap){
    m_patch = &patch;
    m_vtkVertexMap = &vtkVertexMap;
}

/*!
 *  Interface for writing data to stream.
 *
 *  @param[in] stream is the stream to write to
 *  @param[in] name is the name of the data to be written. Either user
 *  data or patch data
 *  @param[in] format is the format which must be used. This flusher write always ascii.
 */
void VTUFlushStreamerASCII::flushData(std::fstream &stream, const std::string &name, bitpit::VTKFormat format)
{
	assert(format == bitpit::VTKFormat::ASCII && m_patch != nullptr);
	BITPIT_UNUSED(format);

	if (name == "Points") {
        auto &verts = m_patch->getVertices();
		for (bitpit::PatchKernel::VertexConstIterator itr = m_patch->vertexConstBegin(); itr != m_patch->vertexConstEnd(); ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap->rawAt(vertexRawId);
			if (vertexVTKId != bitpit::Vertex::NULL_ID) {
				const bitpit::Vertex &vertex = verts.rawAt(vertexRawId);
				bitpit::genericIO::flushASCII(stream, 3, vertex.getCoords());
			}
		}
	} else if (name == "offsets") {
		int offset = 0;
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			offset += cell.getVertexCount();
			bitpit::genericIO::flushASCII(stream, offset);
		}
	} else if (name == "types") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::VTKElementType VTKType;
			switch (cell.getType())  {

			case bitpit::ElementType::VERTEX:
				VTKType = bitpit::VTKElementType::VERTEX;
				break;

			case bitpit::ElementType::LINE:
				VTKType = bitpit::VTKElementType::LINE;
				break;

			case bitpit::ElementType::TRIANGLE:
				VTKType = bitpit::VTKElementType::TRIANGLE;
				break;

			case bitpit::ElementType::PIXEL:
				VTKType = bitpit::VTKElementType::PIXEL;
				break;

			case bitpit::ElementType::QUAD:
				VTKType = bitpit::VTKElementType::QUAD;
				break;

			case bitpit::ElementType::POLYGON:
				VTKType = bitpit::VTKElementType::POLYGON;
				break;

			case bitpit::ElementType::TETRA:
				VTKType = bitpit::VTKElementType::TETRA;
				break;

			case bitpit::ElementType::VOXEL:
				VTKType = bitpit::VTKElementType::VOXEL;
				break;

			case bitpit::ElementType::HEXAHEDRON:
				VTKType = bitpit::VTKElementType::HEXAHEDRON;
				break;

			case bitpit::ElementType::WEDGE:
				VTKType = bitpit::VTKElementType::WEDGE;
				break;

			case bitpit::ElementType::PYRAMID:
				VTKType = bitpit::VTKElementType::PYRAMID;
				break;

			case bitpit::ElementType::POLYHEDRON:
				VTKType = bitpit::VTKElementType::POLYHEDRON;
				break;

			default:
				VTKType = bitpit::VTKElementType::UNDEFINED;
				break;

			}

			bitpit::genericIO::flushASCII(stream, (int) VTKType);
		}
	} else if (name == "connectivity") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
			const int nCellVertices = cell.getVertexCount();
			for (int k = 0; k < nCellVertices; ++k) {
				long vertexId = cellVertexIds[k];
				long vtkVertexId = m_vtkVertexMap->at(vertexId);
				bitpit::genericIO::flushASCII(stream, vtkVertexId);
			}
		}
	} else if (name == "faces") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				bitpit::genericIO::flushASCII(stream, (long) 0);
			} else {
				std::vector<long> faceStream = cell.getFaceStream();
				bitpit::Cell::renumberFaceStream(*(m_vtkVertexMap), &faceStream);
				int faceStreamSize = faceStream.size();
				for (int k = 0; k < faceStreamSize; ++k) {
					bitpit::genericIO::flushASCII(stream, faceStream[k]);
				}
			}
		}
	} else if (name == "faceoffsets") {
		int offset = 0;
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			if (cell.getDimension() <= 2 || cell.hasInfo()) {
				offset += 1;
			} else {
				offset += cell.getFaceStreamSize();
			}

			bitpit::genericIO::flushASCII(stream, offset);
		}
	} else if (name == "cellIndex") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::genericIO::flushASCII(stream, cell.getId());
		}
	} else if (name == "PID") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::genericIO::flushASCII(stream, cell.getPID());
		}
	} else if (name == "vertexIndex") {
		for (bitpit::PatchKernel::VertexConstIterator itr = m_patch->vertexConstBegin(); itr != m_patch->vertexConstEnd(); ++itr) {
			std::size_t vertexRawId = itr.getRawIndex();
			long vertexVTKId = m_vtkVertexMap->rawAt(vertexRawId);
			if (vertexVTKId != bitpit::Vertex::NULL_ID) {
				std::size_t vertexId = itr.getId();
				bitpit::genericIO::flushASCII(stream, vertexId);
			}
		}
#if MIMMO_ENABLE_MPI==1
	} else if (name == "cellGlobalIndex") {
		bitpit::PatchNumberingInfo numberingInfo(m_patch);
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::genericIO::flushASCII(stream, numberingInfo.getCellGlobalId(cell.getId()));
		}
	} else if (name == "rank") {
		for (const bitpit::Cell &cell : m_patch->getVTKCellWriteRange()) {
			bitpit::genericIO::flushASCII(stream, m_patch->getCellRank(cell.getId()));
		}
#endif
	}
}


/*!
 * Base constructor. Linked reference bitpit::PatchKernel container must be empty. If not,
 * class will provide to destroy its previous contents and fill it with new read values.
 * \param[in] streamer streaming class to absorb VTK data.
 * \param[in] patch reference to empty container for storing mesh data.
 * \param[in] eltype [optional] force the elementtype of the grid.
 */
VTUGridWriterASCII::VTUGridWriterASCII( VTUFlushStreamerASCII & streamer, bitpit::PatchKernel & patch, bitpit::VTKElementType eltype) :
                              VTKUnstructuredGrid(eltype), m_patch(patch), m_streamer(streamer)
{
    for(auto & field : m_geometry){
        field.enable();
        field.setStreamer(streamer);
    }

    addData<long>("vertexIndex", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::POINT, &streamer);
    addData<long>("cellIndex", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, &streamer);
    addData<int>("PID", bitpit::VTKFieldType::SCALAR, bitpit::VTKLocation::CELL, &streamer);


    // Get VTK cell count
    long vtkCellCount = 0;
    if (m_patch.getVTKWriteTarget() == bitpit::PatchKernel::WriteTarget::WRITE_TARGET_CELLS_ALL) {
         vtkCellCount = m_patch.getCellCount();
 #if MIMMO_ENABLE_MPI==1
    } else if (m_patch.getVTKWriteTarget() == bitpit::PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL) {
         vtkCellCount = m_patch.getInternalCount();
 #endif
    }

     // Set the dimensinos of the mesh
     bitpit::PiercedStorage<long, long> vertexWriteFlag(1, &(m_patch.getVertices()));
     vertexWriteFlag.fill(false);

     bool vtkFaceStreamNeeded = false;
     for (const bitpit::Cell &cell : m_patch.getCells()) {
         if (cell.getDimension() > 2 && !cell.hasInfo()) {
             vtkFaceStreamNeeded = true;
             break;
         }
     }

     long vtkConnectSize    = 0;
     long vtkFaceStreamSize = 0;
     for (const bitpit::Cell &cell : m_patch.getVTKCellWriteRange()) {
         bitpit::ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
         const int nCellVertices = cellVertexIds.size();
         for (int k = 0; k < nCellVertices; ++k) {
             long vertexId = cellVertexIds[k];
             vertexWriteFlag.at(vertexId) = true;
         }

         vtkConnectSize += nCellVertices;
         if (vtkFaceStreamNeeded) {
             if (cell.getDimension() <= 2 || cell.hasInfo()) {
                 vtkFaceStreamSize += 1;
             } else {
                 vtkFaceStreamSize += cell.getFaceStreamSize();
             }
         }
     }

    int vtkVertexCount = 0;
    m_vtkVertexMap.unsetKernel(true);
    m_vtkVertexMap.setStaticKernel(&(m_patch.getVertices()));
    for (bitpit::PatchKernel::VertexConstIterator itr = m_patch.vertexConstBegin(); itr != m_patch.vertexConstEnd(); ++itr) {
         std::size_t vertexRawId = itr.getRawIndex();
         if (vertexWriteFlag.rawAt(vertexRawId)) {
             m_vtkVertexMap.rawAt(vertexRawId) = vtkVertexCount++;
         } else {
             m_vtkVertexMap.rawAt(vertexRawId) = bitpit::Vertex::NULL_ID;
         }
    }

    setDimensions(vtkCellCount, vtkVertexCount, vtkConnectSize, vtkFaceStreamSize);
    setCodex(bitpit::VTKFormat::ASCII);
    m_streamer.setTargetPatch(m_patch, m_vtkVertexMap);

}

/*!
 * Basic Destructor
 */
VTUGridWriterASCII::~VTUGridWriterASCII(){}

/*!
 * Write to file. Overloading of VTK::write();
  *@param[in] filepath to write.
 */
void VTUGridWriterASCII::write(const std::string & filepath,  bitpit::VTKWriteMode mode){
    setName(filepath);
    VTKUnstructuredGrid::write(filepath, mode);
}




}
