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
#ifndef __MULTIAPPLYDEFORMATION_HPP__
#define __MULTIAPPLYDEFORMATION_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{
/*!
 * \class MultiApply
 * \ingroup manipulators
 * \brief MultiApply is the class that applies one or more deformations resulting from one or more manipulation objects 
 * to different geometries.
 *
 * MultiApply is derived from BaseManipulation class. It uses a map list to retrack each deformation field w.r.t 
 * its associated geometry.
 * The deformation displacements and pointers to geometry structure have to be passed to the class as a 
 * std::pair<MimmoObject*, std::vector<std::array<double,3>>* > or a list of pair (std::unordered_map) through port linking.
 * Different lists added as input to the class in different moments will be appended
 * After the execution of an object MultiApply, the original geometries linked will be 
 * modified (if pointing to non void object).
 *
 * \n
 * Ports available in MultiApply Class :
 *
 *    =========================================================
 * ~~~
     |---------------------------------------------------------------------|
     |                  Port Input                                         |
     |-------|-----------|-------------------|-----------------------------|
     |PortID | PortType  | variable/function | dataType                    |
     |-------|-----------|-------------------|-----------------------------|
     | 107   | M_UMGEOVFD| setInputList      | (UMAP, MIMMO_VECARR3FLOAT_) |
     |-------|-----------|-------------------|-----------------------------|


     |--------------------------------------------------------|
     |                  Port Output                           |
     |-------|----------|-------------------|-----------------|
     |PortID | PortType | variable/function | dataType        |
     |-------|----------|-------------------|-----------------|
     |-------|----------|-------------------|-----------------|
  ~~~
 *    =========================================================
 * \n
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.MultiApply</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>RefreshGeometryTrees</B>: 0/1 to force update of trees for current linked geometries;
 *
 * Geometries input has to be mandatorily passed through port.
 *
 */
class MultiApply: public BaseManipulation{
public:

    std::unordered_map<MimmoObject*, dvecarr3E*> m_input; /**< container for vector field of displacements to apply */
    bool m_force; /**< member to force rebuilding of trees of a geometrical object*/

    MultiApply();
    MultiApply(const bitpit::Config::Section & rootXML);
    ~MultiApply();

    MultiApply(const MultiApply & other);
    MultiApply & operator=(const MultiApply & other);

    void buildPorts();

    bool getRefreshGeometryTrees();
    void setRefreshGeometryTrees(bool force);

    void addInput(std::pair<MimmoObject*, dvecarr3E*> input);
    void setInputList(std::unordered_map<MimmoObject*, dvecarr3E*> input);

    void clearList();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");
};

REGISTER(BaseManipulation, MultiApply, "mimmo.MultiApply")

};

#endif /* __MULTIAPPLYDEFORMATION_HPP__ */
