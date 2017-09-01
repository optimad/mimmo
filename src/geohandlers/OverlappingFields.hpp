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

#ifndef __OVERLAPPINGFIELDS_HPP__
#define __OVERLAPPINGFIELDS_HPP__

#include "ReconstructFields.hpp"

namespace mimmo{

/*!
 * \class OverlapScalarFields
 * \ingroup geohandlers
 * \brief Manipulator that overlaps multiple scalar fields defined on the same mesh objects
 * 
 * Class/BaseManipulation Object overlapping one or more double scalar
 * fields associated to a MimmoObject geometry.
 * The class handles more possible overlapping referring
 * to different geoemtries at the same time.
 * It returns a list of
 * overlapped fields associated to their geometry. 
 * 
 * Ports available in OverlapScalarFields Class :
 * 
 *    =========================================================
 *
     |                   Port Input |||                                                  |
     |-------|----------------|--------------------|----------------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 81    | M_PAIRSCAFIELD | setAddDataField    | (PAIR, MIMMO_VECFLOAT_)          |
     | 106   | M_UMGEOSFD     | setDataFieldMap    | (UMAP, MIMMO_VECFLOAT_)          |
     | 200   | M_VECPAIRSF    | setDataFieldList   | (VECTOR, PAIRMIMMO_VECFLOAT_)    |

     |             Port Output    |||                                                |
     |-------|----------------|--------------------|------------------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 106   | M_UMGEOSFD     | getDataFieldMap    | (UMAP, MIMMO_VECFLOAT_)      |
     | 200   | M_VECPAIRSF    | getDataFieldList   | (VECTOR, PAIRMIMMO_VECFLOAT_)|

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation
 * - <B>ClassName</B>: name of the class as <tt>mimmo.OverlapScalarFields</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: plot optional results in execution;
 * - <B>OutputPlot</B>: path to store optional results.
 *
 * Proper of the class:
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class OverlapScalarFields: public mimmo::BaseManipulation {

private:

    OverlapMethod m_overlapCriterium;                                                   /**<Overlap Method */
    std::unordered_map < mimmo::MimmoObject*, std::vector<dvector1D *> > m_originals;   /**<List of input geometries and related field pointers. */
    std::unordered_map < mimmo::MimmoObject*, dvector1D > m_results;                    /**<List of output geometry pointers and overlapped fields. */

public:
    OverlapScalarFields();
    OverlapScalarFields(const bitpit::Config::Section & rootXML);
    virtual ~OverlapScalarFields();
    OverlapScalarFields(const OverlapScalarFields & other);

    OverlapMethod            getOverlapCriteriumENUM();
    int                     getOverlapCriterium();
    dvector1D                 getResultData(mimmo::MimmoObject * patch );
    int                     getNEffectiveFields();
    int                        getNLinkedFields();

    std::vector<MimmoObject*>                            whichGeometriesLinked();
    std::unordered_map<MimmoObject*, dvector1D* >        getDataFieldMap();
    std::vector<std::pair<MimmoObject*, dvector1D* > >    getDataFieldList();

    void         setOverlapCriteriumENUM( OverlapMethod);
    void         setOverlapCriterium( int);
    void        setAddDataField( std::pair<MimmoObject*, dvector1D*> field );
    void         setDataFieldMap(std::unordered_map<MimmoObject*, dvector1D*> fieldMap );
    void         setDataFieldList(std::vector<std::pair<MimmoObject*, dvector1D*> > fieldList );

    void        removeData(mimmo::MimmoObject* );
    void        removeAllData();

    void        buildPorts();

    void clear();

    void    plotData(std::string dir, std::string name, bool flag, int counter, MimmoObject * geo);
    void    plotAllData(std::string dir, std::string name, bool flag);

    void        execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

    virtual void plotOptionalResults();

private:
    double     overlapFields(dvector1D & locField);
};

/*!
 * 
 * \class OverlapVectorFields
 * \ingroup geohandlers
 * \brief Manipulator that overlaps multiple vector fields defined on the same mesh objects
 * 
 * Class/BaseManipulation Object overlapping one or more 3 double elements vector fields associated to a MimmoObject geometry.
 * The class handles more possible overlapping referring to different geometries at the same time.
 * It returns a list of
 * overlapped fields associated to their geometry.
 * 
 * Ports available in OverlapVectorFields Class :
 * 
 *    =========================================================
 *
     |                   Port Input    |||                                              |
     |-------|----------------|------------------|-----------------------------------|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 80    | M_PAIRVECFIELD | setAddDataField  | (PAIR, MIMMO_VECARR3FLOAT_)       |
     | 107   | M_UMGEOSVD     | setDataFieldMap  | (UMAP, MIMMO_VECARR3FLOAT_)       |
     | 201   | M_VECPAIRVF    | setDataFieldList | (VECTOR, PAIRMIMMO_VECARR3EFLOAT_)|

     |             Port Output      |||                                                   |
     |-------|----------------|--------------------|-----------------------------------|
     |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | 107   | M_UMGEOVFD     | getDataFieldMap    | (UMAP, MIMMO_VECARR3EFLOAT_)      |
     | 201   | M_VECPAIRVF    | getDataFieldList   | (VECTOR, PAIRMIMMO_VECARR3EFLOAT_)|

 *    =========================================================
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation
 * - <B>ClassName</B>: name of the class as <tt>mimmo.OverlapScalarFields</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: plot optional results in execution;
 * - <B>OutputPlot</B>: path to store optional results.
 *
 * Proper of the class:
 * - <B>OverlapCriterium</B>: set how to treat fields in the overlapped region 1-MaxVal, 2-MinVal, 3-AverageVal, 4-Summing;
 *
 * Fields and Geometry have to be mandatorily passed through port.
 *
 */
class OverlapVectorFields: public mimmo::BaseManipulation {

private:

    OverlapMethod m_overlapCriterium;                                                   /**<Overlap Method */
    std::unordered_map < mimmo::MimmoObject*, std::vector<dvecarr3E *> > m_originals;   /**<List of input geometries and related field pointers. */
    std::unordered_map < mimmo::MimmoObject*, dvecarr3E > m_results;                    /**<List of output geometry pointers and overlapped fields. */

public:
    OverlapVectorFields();
    OverlapVectorFields(const bitpit::Config::Section & rootXML);
    virtual ~OverlapVectorFields();
    OverlapVectorFields(const OverlapVectorFields & other);

    OverlapMethod            getOverlapCriteriumENUM();
    int                     getOverlapCriterium();
    dvecarr3E                 getResultData(mimmo::MimmoObject * patch );
    int                     getNEffectiveFields();
    int                        getNLinkedFields();

    std::vector<MimmoObject*>                            whichGeometriesLinked();
    std::unordered_map<MimmoObject*, dvecarr3E* >        getDataFieldMap();
    std::vector<std::pair<MimmoObject*, dvecarr3E* > >    getDataFieldList();

    void         setOverlapCriteriumENUM( OverlapMethod);
    void         setOverlapCriterium( int);
    void        setAddDataField( std::pair<MimmoObject*, dvecarr3E*> field );
    void         setDataFieldMap(std::unordered_map<MimmoObject*, dvecarr3E*> fieldMap );
    void         setDataFieldList(std::vector<std::pair<MimmoObject*, dvecarr3E*> > fieldList );

    void        removeData(mimmo::MimmoObject* );
    void        removeAllData();

    void        buildPorts();

    void clear();

    void    plotData(std::string dir, std::string name, bool flag, int counter, MimmoObject * geo);
    void    plotAllData(std::string dir, std::string name, bool flag);

    void        execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name="");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name="");

protected:

    virtual void plotOptionalResults();

private:
    darray3E     overlapFields(dvecarr3E & locField);
};

REGISTER(BaseManipulation, OverlapScalarFields, "mimmo.OverlapScalarFields")
REGISTER(BaseManipulation, OverlapVectorFields, "mimmo.OverlapVectorFields")

};

#endif /* __OVERLAPPINGFIELDS_HPP__ */
