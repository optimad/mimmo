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
#ifndef __CGNSPIDEXTRACTOR_HPP__
#define __CGNSPIDEXTRACTOR_HPP__

#include "BaseManipulation.hpp"
#include <set>

namespace mimmo{

/*!
 *  \class    CGNSPidExtractor
 *  \ingroup iocgns
 *  \brief CGNSPidExtractor is the class to extract a 3D surface by means of PIDs from
 *   a whole surface boundary of a volume mesh read from cgns file.
 *
 * The object return the patch with a target pid (or pids) as an independent surface mesh.
 * It homogeneously triangulates the patch if required.
 * Please notice if the extracted sub-patch is forced to be triangulated,
 * no link with the mother boundary surface mesh could be guaranteed anymore.
 *
 * Ports available in CGNSPidExtractor Class :
 *
 * =========================================================

   |                     Port Input ||                                       |
   |------------------|---------------------|----------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | setGeometry         | (SCALAR, MIMMO_)     |


   |               Port Output    ||                                           |
   |------------------|---------------------|------------------------|
   | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | M_GEOM           | getPatch()          | (SCALAR, MIMMO_)       |

 * =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.CGNSPidExtractor</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 *Proper of the class
 * - <B>nPID</B>: number of Geometry PID involved in the extraction;
 * - <B>PID</B>: list of Geometry PID (blank-separated) that need to be extracted(element must be coeherent with nPID);
 * - <B>ForcedToTriangulate</B>: force retriangulation of the extracted patches, boolean 0/1;
 *
 * Geometry has to be mandatorily read or passed through port.
 *
 */
class CGNSPidExtractor: public BaseManipulation{
private:
    bool                         m_force;     /**<force triangulation of the extracted patch if true.*/
    std::set<long>               m_targetpid; /**<list of PID involved in the extraction.*/
    std::unique_ptr<MimmoObject> m_patch;     /**<extracted patch */

public:
    CGNSPidExtractor();
    CGNSPidExtractor(const bitpit::Config::Section & rootXML);
    ~CGNSPidExtractor();

    CGNSPidExtractor(const CGNSPidExtractor & other);
    CGNSPidExtractor & operator=(CGNSPidExtractor other);

    void            buildPorts();

    MimmoObject*    getPatch();
    std::set<long>  whatPIDActive();
    bool            isForcedToTriangulate();

    void            addPID(long val);
    void            setPID(std::vector<long> vval);
    void            setForcedToTriangulate(bool flag);
    void            setGeometry(MimmoObject*);

    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:
    void swap(CGNSPidExtractor & x) noexcept;
    void plotOptionalResults();

};

REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__CGNSPIDEXTRACTOR_HPP__)

REGISTER(BaseManipulation, CGNSPidExtractor,"mimmo.CGNSPidExtractor")
}

#endif /* __CGNSPIDEXTRACTOR_HPP__ */
