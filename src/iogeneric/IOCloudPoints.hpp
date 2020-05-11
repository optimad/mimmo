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

#ifndef __IOCLOUDPOINTS_HPP__
#define __IOCLOUDPOINTS_HPP__

#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class IOCloudPoints
 * \ingroup iogeneric
 * \brief IOCloudPoints is the class to read from file a set of cloud 3D points w/ attached
 * a scalar field of floats and/or a vector field of floats
 *
 * The only admissible File format is an ascii list of values, organized as follow:
 *
 * <tt>
 * <B>$POINT</B>    l1  0.0 0.0  1.0 \n
 * <B>$POINT</B>   l2 -1.0 0.12 0.0 \n
 * <B>...</B>
 * </tt>
 *
 * where <B>$POINT</B> keyword identify the row relative to a single point,
   l1, l2,... the unique int label associated to the point and the following 3 floats
   represent the point coordinate. If <B>$POINT</B> is missing, the point will not be read.
 * After all points declaration, to set a scalar value on a point define:
 *
 * <tt>
 * <B>$SCALARF</B> l1 12.0 \n
 * <B>$SCALARF</B> l2 -4.232 \n
 * <B>...</B>
 * </tt>
 *
 * where l1, l2, are still the unique labels of points.
   Similarly for vector values on points, define:
 *
 * <tt>
 * <B>$VECTORF</B> l1 3.0 2.1 3.3 \n
 * <B>$VECTORF</B> l2 -4.2 0.0 0.0 \n
 * <B>...</B>
 * </tt>
 *
 * Missing keywords or point without field defined will be considered
   at values {0.0} or {0.0,0.0,0.0};
 *
 * IOCloudPoints is derived from BaseManipulation class. The class working in
   both Read and Write mode, that is can read from or write to file, provided
   that its format requirements are met.
 * When in write mode the class can generate a template file for both scalar and
   vector fields, that can be filled in a second moment for different purposes.
 * The layout of this file will be:
 *
 * <tt>
 * <B>$SCALARF</B>    l1  {sl1} \n
 * <B>$VECTORF</B> l2  {xl2} {yl2}  {zl2} \n
 * <B>...</B>
 * </tt>
 *
 * where {xxx} uniquely naming the component of displacement.
 *
 * The point cloud is provided as a MimmoObject point cloud.
 * The geometry can be stored internally or given by an external block by set Geometry method/port.
 *
 *
 * \n
 * Ports available in IOCloudPoints Class :
 *
 *    =========================================================

     |                 Port Input   ||                                     |
     |---------------|-------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM        | setGeometry       | (MC_SCALAR,MD_MIMMO_)     |
     | M_SCALARFIELD | setScalarField | (MC_SCALAR,MD_MPVECFLOAT_)     |
     | M_VECTORFIELD | setVectorField | (MC_SCALAR,MD_MPVECARR3FLOAT_) |


     |              Port Output  ||                                        |
     |---------------|-------------------|-----------------------|
     | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
     | M_GEOM        | getGeometry       | (MC_SCALAR,MD_MIMMO_)     |
     | M_SCALARFIELD | getScalarField | (MC_SCALAR,MD_MPVECFLOAT_)     |
     | M_VECTORFIELD | getVectorField | (MC_SCALAR,MD_MPVECARR3FLOAT_) |

 *    =========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * Inherited from BaseManipulation:
 * - <B>ClassName</B>: name of the class as <tt>mimmo.IOCloudPoints</tt>;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>PlotInExecution</B>: boolean 0/1 print optional results of the class;
 * - <B>OutputPlot</B>: target directory for optional results writing.
 *
 * Proper of the class:
 * - <B>IOmode</B>: 1/0 enable Read and Write mode,respectively;
 * - <B>ReadDir</B>: path to input directory in read mode;
 * - <B>ReadFilename</B>: name of input file with tag extension in read mode;
 * - <B>WriteDir</B>: path to output directory in write mode;
 * - <B>WriteFilename</B>: name of output file with tag extension in write mode;
 * - <B>Template</B>: 0/1 option to activate writing file in template mode;
 *
 */
class IOCloudPoints: public BaseManipulation{
protected:
    //members
    bool            m_read;         /**<True if in Read mode, False if in Write mode.*/
    std::string     m_dir;          /**<Directory path for I/O*/
    std::string     m_filename;     /**<I/O filename with extension tag*/
    bool            m_template;     /**<True/False enable the writing template mode */

    dmpvector1D     m_scalarfield;  /**< MimmoPiercedVector scalar field */
    dmpvecarr3E     m_vectorfield;  /**< MimmoPiercedVector vector field */

    bool                            m_isInternal;   /**< Flag for internal instantiated MimmoObject */
    std::unique_ptr<MimmoObject>    m_intgeo;       /**< Pointer to internal allocated geometry, if any */

    livector1D      m_labels;   /**< Labels associated to displacement, for internal use. */
    dvecarr3E       m_points;   /**< Cloud points list, for internal use. */

public:
    IOCloudPoints(bool readMode = true);
    virtual ~IOCloudPoints();
    IOCloudPoints(const bitpit::Config::Section & rootXML);
    IOCloudPoints(const IOCloudPoints & other);
    IOCloudPoints & operator=(IOCloudPoints other);

    void buildPorts();

    bool            isTemplate();

    MimmoObject* getGeometry();
    dmpvector1D* getScalarField();
    dmpvecarr3E* getVectorField();

    void setReadDir(std::string dir);
    void setReadFilename(std::string filename);
    void setWriteDir(std::string dir);
    void setWriteFilename(std::string filename);
    void setTemplate(bool flag);

    void setGeometry(MimmoObject* geometry);
    void setScalarField(dmpvector1D* scalarfield);
    void setVectorField(dmpvecarr3E* vectorfield);

    void clear();

    void execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");
    void plotOptionalResults();

protected:
    void swap(IOCloudPoints & x) noexcept;

private:
    virtual void read();
    virtual void write();
};

REGISTER_PORT(M_SCALARFIELD, MC_SCALAR, MD_MPVECFLOAT_,__IOCLOUDPOINTS_HPP__)
REGISTER_PORT(M_VECTORFIELD, MC_SCALAR, MD_MPVECARR3FLOAT_,__IOCLOUDPOINTS_HPP__)
REGISTER_PORT(M_GEOM, MC_SCALAR, MD_MIMMO_ ,__IOCLOUDPOINTS_HPP__)

REGISTER(BaseManipulation, IOCloudPoints, "mimmo.IOCloudPoints")
}

#endif /* __IOCLOUDPOINTS_HPP__ */
