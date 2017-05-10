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
 *  \brief CGNSPidExtractor is the class to extract a with pid (or pids) from
 *   a surface boundary mesh read form cgns file.
 *
 * The object return the patch with a target pid (or pids) as an independent surface mesh.
 * It triangulates the patch if required (currently only for elements of type QUAD).
 *
 * Ports available in CGNSPidExtractor Class :
 *
 * =========================================================
 * ~~~
 *  |-----------------------------------------------------------------------|
 *  |                     Port Input                                        |
 *  |-------|------------------|---------------------|----------------------|
 *  |PortID | PortType         | variable/function   | DataTypes            |
 *  |-------|------------------|---------------------|----------------------|
 *  | 99    | M_GEOM           | setGeometry         | (SCALAR, MIMMO_)     |
 *  |-------|------------------|---------------------|----------------------|
 *
 *
 *  |-------------------------------------------------------------------------|
 *  |               Port Output                                               |
 *  |-------|------------------|---------------------|------------------------|
 *  |PortID | PortType         | variable/function   | DataTypes              |
 *  |-------|------------------|---------------------|------------------------|
 *  | 99    | M_GEOM           | getPatch()          | (SCALAR, MIMMO_)       |
 *  |-------|------------------|---------------------|------------------------|
 * ~~~
 * =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.CGNSPidExtractor</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>nPID</B>: number of Geometry PID involved in the extraction;
 * - <B>PID</B>: list of Geometry PID (blank-separated) that need to be extracted(element must be coeherent with nPID);
 * - <B>ForcedToTriangulate</B>: force retriangulation of the extracted patches, boolean 0/1;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Geometry has to be mandatorily read or passed through port.
 *
 */
class CGNSPidExtractor: public BaseManipulation{
private:
    bool                         m_force;     /**<force triangulation of the extracted patch if true.*/
    std::set<short>              m_targetpid; /**<list of PID involved in the extraction.*/
    std::unique_ptr<MimmoObject> m_patch;     /**<extracted patch */

public:
    CGNSPidExtractor();
    CGNSPidExtractor(const bitpit::Config::Section & rootXML);
    ~CGNSPidExtractor();

    CGNSPidExtractor(const CGNSPidExtractor & other);
    CGNSPidExtractor & operator=(const CGNSPidExtractor & other);

    void            buildPorts();

    MimmoObject*    getPatch();
    std::set<short> whatPIDActive();
    bool            isForcedToTriangulate();

    void            addPID(short val);
    void            setPID(std::vector<short int> vval);
    void            setForcedToTriangulate(bool flag);
    void            setGeometry(MimmoObject*);

    void            execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

    void plotOptionalResults();

};

REGISTER(BaseManipulation, CGNSPidExtractor,"mimmo.CGNSPidExtractor")
}

#endif /* __CGNSPIDEXTRACTOR_HPP__ */
