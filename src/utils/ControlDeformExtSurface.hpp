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
#ifndef __CONTROLDEFORMEXTSURFACE_HPP__
#define __CONTROLDEFORMEXTSURFACE_HPP__

#include "BaseManipulation.hpp"
#include "MimmoGeometry.hpp"

namespace mimmo{

/*!
 * \class ControlDeformExtSurface
 * \ingroup utils
 * \brief ControlDeformExtSurface is a class that check a deformation field,
          associated to a MimmoObject geometry, for eventual collisions/penetrations
          w.r.t. one or more external constraint surface meshes.
 *
 * ControlDeformExtSurface is derived from BaseManipulation class.
 * It needs one or more external surface meshes, representing the
   constraint of your deformed object.
 * It returns a double value V, namely the maximum signed distance from constraint
   surfaces amongst all deformed geometry nodes, reporting how much the current
   deformation violates the constraint itself.
 * if V > 0 a violation occurs. If V=0, a contact occurs, otherwise if V<0
   no violation occurs.
 * Need to link the deformation field in exam, as MimmoPiercedVector referred to
   target geometry and defined to POINT as location.
 * Class absorbs/flushes its parameters from/to xml dictionaries
 *
 * \n
 * Ports available in ControlDeformExtSurface Class :
 *
 *    =========================================================

    | Port Input | | |
    |-|-|-|
    | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
    | M_GDISPLS| setDefField       | (MC_SCALAR, MD_MPVECARR3FLOAT_)       |
    | M_GEOM   | setGeometry       | (MC_SCALAR, MD_MIMMO_)        |


    |Port Output | | |
    |-|-|-|
    | <B>PortType</B> | <B>variable/function</B> |<B>DataType</B>              |
    | M_SCALARFIELD | getViolationField | (MC_SCALAR, MD_MPVECFLOAT_)             |
    | M_VALUED      | getViolation      | (MC_SCALAR, MD_FLOAT)             |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.ControlDeformExtSurface</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>Files</B>: external constraint surfaces list of file: \n\n
 *   <tt> <B>\<Files\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<file0\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<fullpath\></B> full path to file <B>\</fullpath\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<tag\></B> tag extension of your file (see MimmoGeometry general doxy) <B>\</tag\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<tolerance\></B> offset to control violation (float value) <B>\</tolerance\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\</file0\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\<file1\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<fullpath\></B> full path to file <B>\</fullpath\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<tag\></B> tag extension of your file (see MimmoGeometry general doxy) <B>\</tag\></B> \n
 *   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<B>\<tolerance\></B> offset to control violation (float value) <B>\</tolerance\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>\</file1\></B> \n
 *   &nbsp;&nbsp;&nbsp;<B>...</B> \n
 *   &nbsp;&nbsp;&nbsp;<B>...</B> \n
 *   <B>\</Files\></B> </tt> \n\n
 * - <B>BGDetails</B>: OPTIONAL define spacing of background grid, dividing diagonal of AABBox containing geometries by this factor (int);
 *
 * Geometry and deformation field have to be mandatorily passed through port.
 *
 */
class ControlDeformExtSurface: public BaseManipulation{
public:
    /*!
    \ingroup typedefs
    \{
    */
    typedef std::unordered_map<std::string, std::pair<double, int> > fileListWithType; /**< custom typedef for a structure hosting a list of filenames with a double and int values attached on it*/
    /*
    \}
    */

private:
    fileListWithType        m_geolist;        /**< list of file for geometrical proximity check*/
    dmpvector1D             m_violationField; /**<Violation Field as distance from constraint */
    dmpvecarr3E             m_defField;       /**<Deformation field*/
    int                     m_cellBackground; /**< Number of cells N to determine background grid spacing */
    std::unordered_set<int> m_allowed;        /**< list of currently file format supported by the class*/

public:
    ControlDeformExtSurface();
    ControlDeformExtSurface(const bitpit::Config::Section & rootXML);
    virtual ~ControlDeformExtSurface();

    ControlDeformExtSurface(const ControlDeformExtSurface & other);
    ControlDeformExtSurface & operator=(ControlDeformExtSurface other);

    void    buildPorts();

    double                        getViolation();
    dmpvector1D *                 getViolationField();
    double                        getToleranceWithinViolation(std::string);
    int                           getBackgroundDetails();
    const    fileListWithType &   getFiles() const;

    void     setDefField(dmpvecarr3E *field);
    void     setGeometry(MimmoObject * geo);
    void     setBackgroundDetails(int nCell=50);
    void     setFiles(fileListWithType list );
    void     addFile(std::string file, double tol, int format);
    void     addFile(std::string file, double tol, FileType format);
    void     removeFile(std::string);
    void     removeFiles();

    void    clear();
    void    execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

protected:
    void plotOptionalResults();
    void swap(ControlDeformExtSurface & x) noexcept;

private:
    void readGeometries(std::vector<std::unique_ptr<MimmoGeometry> > & extGeo, std::vector<double> & tols);
    svector1D extractInfo(std::string file);
    double evaluateSignedDistance(darray3E &point, mimmo::MimmoObject * geo, long & id, darray3E & normal, double &initRadius);
    void writeLog();
};

REGISTER_PORT(M_GDISPLS, MC_SCALAR, MD_MPVECARR3FLOAT_,__CONTROLDEFORMEXTSURFACE_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_,__CONTROLDEFORMEXTSURFACE_HPP__)
REGISTER_PORT(M_VALUED, MC_SCALAR, MD_FLOAT,__CONTROLDEFORMEXTSURFACE_HPP__)
REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__CONTROLDEFORMEXTSURFACE_HPP__)


REGISTER(BaseManipulation, ControlDeformExtSurface,"mimmo.ControlDeformExtSurface")

}

#endif /* __CONTROLDEFORMEXTSURFACE_HPP__ */
