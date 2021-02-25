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
 \ *---------------------------------------------------------------------------*/

#include "AABBox.hpp"

namespace mimmo{

/*! Base Constructor. Doing nothing.*/
AABBox::AABBox(){
    m_name = "mimmo.AABBox";
    m_origin.fill(0.0);
    m_span.fill(1.0);
    int counter = 0;
    for(auto &val : m_axes)    {
        val.fill(0.0);
        val[counter] = 1.0;
        ++counter;
    }
    m_writeInfo = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
AABBox::AABBox(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.AABBox";
    m_origin.fill(0.0);
    m_span.fill(1.0);
    int counter = 0;
    for(auto &val : m_axes)    {
        val.fill(0.0);
        val[counter] = 1.0;
        ++counter;
    }
    m_writeInfo = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.AABBox"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*! Destructor */
AABBox::~AABBox(){};

/*! Copy Constructor
 *\param[in] other AABBox where copy from
 */
AABBox::AABBox(const AABBox & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_span = other.m_span;
    m_axes = other.m_axes;
    m_listgeo = other.m_listgeo;
    m_writeInfo = other.m_writeInfo;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void AABBox::swap(AABBox & x) noexcept
{
    std::swap(m_origin, x.m_origin);
    std::swap(m_span, x.m_span);
    std::swap(m_axes, x.m_axes);
    std::swap(m_listgeo, x.m_listgeo);
    std::swap(m_writeInfo, x.m_writeInfo);

    BaseManipulation::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
AABBox::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoSharedPointer<MimmoObject>, AABBox>(this, &mimmo::AABBox::setGeometry, M_GEOM,true,1));
    built = (built && createPortIn<std::vector<MimmoSharedPointer<MimmoObject> >, AABBox>(this, &mimmo::AABBox::setGeometries, M_VECGEOM, true,1));
    built = (built && createPortIn<dmatrix33E, AABBox>(this, &mimmo::AABBox::setAxes, M_AXES));
    built = (built && createPortOut<darray3E, AABBox>(this, &mimmo::AABBox::getOrigin, M_POINT));
    built = (built && createPortOut<dmatrix33E, AABBox>(this, &mimmo::AABBox::getAxes, M_AXES));
    built = (built && createPortOut<darray3E, AABBox>(this, &mimmo::AABBox::getSpan, M_SPAN));

    m_arePortsBuilt = built;
};

/*!Clean all stuffs in your class */
void
AABBox::clearAABBox(){
    clear(); //base manipulation stuff clear
    m_origin.fill(0.0);
    m_span.fill(1.0);
    m_listgeo.clear();
    int counter = 0;
    for(auto &val : m_axes)    {
        val.fill(0.0);
        val[counter] = 1.0;
        ++counter;
    }
     m_writeInfo = false;
};

/*!
 * \return the origin of the AABB.
 */
darray3E
AABBox::getOrigin(){
    return(m_origin);
}

/*!
 * \return the span of the AABB.
 */
darray3E
AABBox::getSpan(){
    return(m_span);
}


/*!
 * \return the axes sdr of the AABB (default or provided by the User).
 */
dmatrix33E
AABBox::getAxes(){

    return(m_axes);
}


/*!
 * \return current list of linked geometry for aabbox calculation
 */
std::vector<MimmoSharedPointer<MimmoObject> >
AABBox::getGeometries(){
    std::vector<MimmoSharedPointer<MimmoObject> > result;
    result.reserve(m_listgeo.size());
    for(auto pp: m_listgeo)
        result.push_back(pp.first);
    return result;
};


/*!
 * Set the list of target geometries once and for all, and erase any pre-existent list.
 * \param[in] listgeo list of pointers to target geometries
 */
void
AABBox::setGeometries(std::vector<MimmoSharedPointer<MimmoObject> > listgeo){
    m_listgeo.clear();
    for( auto val : listgeo){
        setGeometry(val);
    }
};

/*!
 * * Add your target geometry to the list of target geometries.
 * Reimplemented from BaseManipulation::setGeometry().
 * \param[in] geo pointer to target geometry
 */
void
AABBox::setGeometry(MimmoSharedPointer<MimmoObject> geo){

    if (geo == nullptr)    {
        (*m_log)<<"warning: "<<m_name<<" not valid Geometry pointer. Doing nothing"<<std::endl;
        return;
    }

    m_listgeo.insert(std::make_pair(geo, geo->getType()));
};


/*!
 * set the new s.d.r. of axes to calculate AABB.
   \param[in] axes new sdr of axes (different from canonical x,y,z)
 */
void
AABBox::setAxes(dmatrix33E axes){
    m_axes = axes;
}

/*!
 * Force class to write AABB calculated origin, axes and span to file(in the plotOptionalResults directory)
 * \param[in] flag true, write the AABB info on file.
 */
void
AABBox::setWriteInfo(bool flag){
    m_writeInfo = flag;
}

/*! Plot the AABB as a structured grid to *vtu file.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary    boolean flag for 0-"ascii" or 1-"appended" writing
 */
void
AABBox::plot(std::string directory, std::string filename,int counter, bool binary){


    dvecarr3E activeP(8);

    activeP[0] =  - 0.5 * m_span;
    activeP[6] =    0.5 * m_span;

    activeP[1] = activeP[0]; activeP[1][0] += m_span[0];
    activeP[3] = activeP[0]; activeP[3][1] += m_span[1];
    activeP[2] = activeP[6]; activeP[2][2] += -1.0*m_span[2];

    activeP[7] = activeP[6]; activeP[7][0] += -1.0*m_span[0];
    activeP[5] = activeP[6]; activeP[5][1] += -1.0*m_span[1];
    activeP[4] = activeP[0]; activeP[4][2] += m_span[2];

    darray3E temp;
    dmatrix33E    inv = inverse(m_axes);
    for(auto &val : activeP){
        for(int i=0; i<3; ++i){
            temp[i] = dotProduct(val, inv[i]);
        }
        val = temp + m_origin;
    }

    ivector2D activeConn(1);
    for(int i=0; i<8; ++i)    activeConn[0].push_back(i);

    bitpit::VTKFormat codex = bitpit::VTKFormat::ASCII;
    if(binary){codex=bitpit::VTKFormat::APPENDED;}
    bitpit::VTKElementType elDM = bitpit::VTKElementType::HEXAHEDRON;
    bitpit::VTKUnstructuredGrid vtk(directory, filename, elDM);
    vtk.setGeomData( bitpit::VTKUnstructuredField::POINTS, activeP) ;
    vtk.setGeomData( bitpit::VTKUnstructuredField::CONNECTIVITY, activeConn) ;
    vtk.setDimensions(1, 8);
    vtk.setCodex(codex);
    if(counter>=0){vtk.setCounter(counter);}

    vtk.write();
};


/*!Execute your object, calculate the AABBox of your geometry.
 * Implementation of pure virtual BaseManipulation::execute
 */
void
AABBox::execute(){

    // check if the list is empty.
    if(m_listgeo.empty()){
        m_span.fill(0.0);
        m_origin.fill(0.0);
        return;
    };

    //check sdr;
    {
        double normx1= norm2(m_axes[0]);
        double normx2= norm2(m_axes[1]);
        double tol = std::numeric_limits<double>::min();
        if (normx1 > tol && normx2 > tol){
            m_axes[0] /= normx1;
            m_axes[1] /= normx2;
        }else{
            m_axes[0] = {{0,0,1}};
            m_axes[1] = {{0,1,0}};
        }

        m_axes[2] = crossProduct(m_axes[0], m_axes[1]);
    }

    darray3E pmin, pmax;
    pmin.fill(std::numeric_limits<double>::max());
    pmax.fill(-1.0*std::numeric_limits<double>::max());
    double val;

    for(auto ptr : getGeometries()){
        for(auto & vert: ptr->getVertices()){
            darray3E coord = vert.getCoords();
            for(int i=0;i<3; ++i){
                val = dotProduct(coord, m_axes[i]);
                pmin[i] = std::fmin(pmin[i], val);
                pmax[i] = std::fmax(pmax[i], val);
            }
        }
    }

#if MIMMO_ENABLE_MPI
    MPI_Allreduce(MPI_IN_PLACE, pmin.data(), 3, MPI_DOUBLE, MPI_MIN, m_communicator);
    MPI_Allreduce(MPI_IN_PLACE, pmax.data(), 3, MPI_DOUBLE, MPI_MAX, m_communicator);
#endif

    dmatrix33E inv = inverse(m_axes);
    m_span = pmax - pmin;
    //check if one of the span goes to 0;
    double avg_span = 0.0;
    for(auto & val: m_span)    avg_span+=val;
    avg_span /= 3.0;

    for(auto &val : m_span)    {
        val = std::fmax(val, 1.E-04*avg_span);
    }

    darray3E originLoc = 0.5*(pmin+pmax);
    for(int i=0; i<3; ++i){
        m_origin[i] = dotProduct(originLoc, inv[i]);
    }


    if (m_writeInfo){

        std::ofstream out;
        out.open(m_outputPlot+"/"+m_name+std::to_string(getId())+"_INFO.dat");
        if(out.is_open()){
            out<<"AABBox "<<std::to_string(getId())<<" info:"<<std::endl;
            out<<std::endl;
            out<<std::endl;
            out<<"Origin: "<<std::scientific<<m_origin<<std::endl;
            out<<std::endl;
            out<<"Axis 0: "<<std::scientific<<m_axes[0]<<std::endl;
            out<<"Axis 1: "<<std::scientific<<m_axes[1]<<std::endl;
            out<<"Axis 2: "<<std::scientific<<m_axes[2]<<std::endl;
            out<<std::endl;
            out<<"Span:   "<<std::scientific<<m_span<<std::endl;

            out.close();
        }
    }
};

/*!
 * Plot Optional results of the class,
 * that is the oriented bounding box as *.vtu mesh
 */
void
AABBox::plotOptionalResults(){
    std::string dir = m_outputPlot ;
    std::string nameGrid  = m_name;
    plot(dir, nameGrid, getId(), true );
}

/*!
 * It sets infos reading from a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
AABBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("WriteInfo")){
        std::string input = slotXML.get("WriteInfo");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setWriteInfo(value);
    }

    if(slotXML.hasSection("Axes")){

        const bitpit::Config::Section & axesXML = slotXML.getSection("Axes");
        dmatrix33E axes;
        for(int i=0; i<3; ++i){
            axes[i].fill(0.0);
            axes[i][i] = 1.0;
        }

        if(axesXML.hasOption("axis0")){
            std::string input = axesXML.get("axis0");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[0][0]>>axes[0][1]>>axes[0][2];
            }
        }

        if(axesXML.hasOption("axis1")){
            std::string input = axesXML.get("axis1");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[1][0]>>axes[1][1]>>axes[1][2];
            }
        }

        if(axesXML.hasOption("axis2")){
            std::string input = axesXML.get("axis2");
            input = bitpit::utils::string::trim(input);
            if(!input.empty()){
                std::stringstream ss(input);
                ss>>axes[2][0]>>axes[2][1]>>axes[2][2];
            }
        }

        setAxes(axes);
    }
}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
AABBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("WriteInfo", std::to_string((int)m_writeInfo));
    {
        dmatrix33E axes = getAxes();
        bitpit::Config::Section & axesXML = slotXML.addSection("Axes");

        for(int i=0; i<3; ++i){
            std::string name = "axis"+std::to_string(i);
            std::stringstream ss;
            ss<<std::scientific<<axes[i][0]<<'\t'<<axes[i][1]<<'\t'<<axes[i][2];
            axesXML.set(name, ss.str());
        }
    }
};

/*!
    Transpose a 3x3 double matrix
    \param[in] mat target matrix
    \return new matrix transposed
*/
dmatrix33E AABBox::transpose(const dmatrix33E & mat){
    dmatrix33E out;

    for(std::size_t i=0; i<3; ++i){
        for(std::size_t j=0; j<3; ++j){
            out[j][i] = mat[i][j];
        }
    }
    return out;
}

/*!
    Invert a 3x3 double matrix
    \param[in] mat target matrix
    \return new matrix transposed
*/
dmatrix33E AABBox::inverse(const dmatrix33E & mat){
    dmatrix33E out;

    double det = mat[0][0] * (mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2]) -
                 mat[0][1] * (mat[1][0]*mat[2][2] - mat[2][0]*mat[1][2]) +
                 mat[0][2] * (mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1]);

    out[0][0] = (mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])/det;
    out[0][1] = (mat[0][2]*mat[2][1] - mat[2][2]*mat[0][1])/det;
    out[0][2] = (mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2])/det;
    out[1][0] = (mat[1][2]*mat[2][0] - mat[2][2]*mat[1][0])/det;
    out[1][1] = (mat[0][0]*mat[2][2] - mat[2][0]*mat[0][2])/det;
    out[1][2] = (mat[0][2]*mat[1][0] - mat[1][2]*mat[0][0])/det;
    out[2][0] = (mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1])/det;
    out[2][1] = (mat[0][1]*mat[2][0] - mat[2][1]*mat[0][0])/det;
    out[2][2] = (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1])/det;


    return out;
}
}
