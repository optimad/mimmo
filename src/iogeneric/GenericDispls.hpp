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
#ifndef __GENERICDISPLS_HPP__
#define __GENERICDISPLS_HPP__

#include <string>
#include "BaseManipulation.hpp"

namespace mimmo{

/*!
 * \class GenericDispls
 * \ingroup iogeneric
 * \brief GenericDispls is the class to read from file an initial set of displacements 
 * as a generic vector field of floats or write it to file
 * 
 * The only admissible File format is an ascii list of values, organized as follow:
 * 
 * <tt>
 * $DISPL   l1  0.0 0.0  1.0 \n
 * $DISPL   l2 -1.0 0.12 0.0 \n
 * ...
 * </tt>
 * 
 * The $DISPL keyword identify the value, l1, l2 the unique int label associated to the displacement and the 
 * following 3 vector coordinate represents the entity of the displacement. If $DISPL is missing, the value will not be read.
 * 
 * GenericDispls is derived from BaseManipulation class. The class working in both Read and Write mode, that is can 
 * read displacement values from file (written in the proper format) or write them on it;
 * When in write mode the class can generate a template file for displacements, that can be filled in a second moment for different purposes.
 * The layout of this file will be:
 * 
 * <tt>
 * $DISPL   l1  {xl1} {yl1}  {zl1} \n
 * $DISPL   l2  {xl2} {yl2}  {zl2} \n
 * ...
 * </tt>
 * 
 * where {xxx} uniquely naming the component of displacement
 *
 * \n
 * Ports available in GenericDispls Class :
 *
 * =========================================================

   |                 Port Input   |||                                  |
   |-------|------------|-------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 10    | M_DISPLS   | setDispl          | (VECARR3, FLOAT)      |
   | 18    | M_VECTORLI | setLabels         | (VECTOR, LONG)        |
   | 31    | M_VALUEI   | setNDispl         | (SCALAR, INT)         |


   |              Port Output  |||                                        |
   |-------|---------------|-------------------|-----------------------|
    |<B>PortID</B> | <B>PortType</B>   | <B>variable/function</B>  |<B>DataType</B> |
   | 10    | M_DISPLS      | getDispl          | (VECARR3, FLOAT)      |
   | 18    | M_VECTORLI    | getLabels         | (VECTOR, LONG)        |
   | 31    | M_VALUEI      | getNDispl         | (SCALAR, INT)         |

 *
 *=========================================================
 * \n
 *
 * The xml available parameters, sections and subsections are the following :
 *
 * - <B>ClassName</B>: name of the class as <tt>mimmo.GenericDispls</tt>;
 * - <B>IOmode</B>: 1/0 enable Read and Write mode,respectively;
 * - <B>Priority</B>: uint marking priority in multi-chain execution;
 * - <B>ReadDir</B>: path to input directory in read mode;
 * - <B>ReadFilename</B>: name of input file with tag extension in read mode;
 * - <B>WriteDir</B>: path to output directory in write mode;
 * - <B>WriteFilename</B>: name of output file with tag extension in write mode;
 * - <B>NDispl</B>: path to your current file dataTAG;
 * - <B>Template</B>: path to your current file data;
 *
 */
class GenericDispls: public BaseManipulation{
protected:
    bool            m_read;     /**<True if in Read mode, False if in Write mode.*/
    std::string     m_dir;      /**<Source/Destination directory path*/
    std::string     m_filename; /**<Source/Destination filename with extension tag*/
    int             m_nDispl;   /**<Number of displacement hold by the class */
    dvecarr3E       m_displ;    /**<Displacement list*/
    livector1D      m_labels;   /**<Labels associated to displacement */
    bool            m_template; /**<True/False enable the writing template mode */

public:
    GenericDispls(bool readMode = true);
    virtual ~GenericDispls();
    GenericDispls(const bitpit::Config::Section & rootXML);
    GenericDispls(const GenericDispls & other);
    GenericDispls & operator=(const GenericDispls & other);

    void buildPorts();

    int         getNDispl();
    dvecarr3E   getDispl();
    livector1D  getLabels();
    bool        isTemplate();

    void setReadDir(std::string dir);
    void setReadFilename(std::string filename);
    void setWriteDir(std::string dir);
    void setWriteFilename(std::string filename);
    void setNDispl(int nD);
    void setLabels(livector1D labels);
    void setDispl(dvecarr3E displs);
    void setTemplate(bool flag);

    void    clear();
    void    execute();

    virtual void absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name = "");
    virtual void flushSectionXML(bitpit::Config::Section & slotXML, std::string name= "");

private:
    virtual void read();
    virtual void write();
};

REGISTER(BaseManipulation, GenericDispls, "mimmo.GenericDispls")
}

#endif /* __GENERICDISPLS_HPP__ */
