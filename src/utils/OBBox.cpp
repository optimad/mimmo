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

#include "OBBox.hpp"
#include <lapacke.h>

#include <chrono>

using namespace std::chrono;

using namespace std;
namespace mimmo{

/*! Basic Constructor. Doing nothing.*/
OBBox::OBBox(){
    m_name = "mimmo.OBBox";
    m_origin.fill(0.0);
    m_span.fill(1.0);
    int counter = 0;
    for(auto &val : m_axes)    {
        val.fill(0.0);
        val[counter] = 1.0;
        ++counter;
    }
    m_forceAABB = false;
    m_writeInfo = false;
};

/*!
 * Custom constructor reading xml data
 * \param[in] rootXML reference to your xml tree section
 */
OBBox::OBBox(const bitpit::Config::Section & rootXML){

    m_name = "mimmo.OBBox";
    m_origin.fill(0.0);
    m_span.fill(1.0);
    int counter = 0;
    for(auto &val : m_axes)    {
        val.fill(0.0);
        val[counter] = 1.0;
        ++counter;
    }
    m_forceAABB = false;
    m_writeInfo = false;

    std::string fallback_name = "ClassNONE";
    std::string input = rootXML.get("ClassName", fallback_name);
    input = bitpit::utils::string::trim(input);
    if(input == "mimmo.OBBox"){
        absorbSectionXML(rootXML);
    }else{
        warningXML(m_log, m_name);
    };
}

/*! Destructor */
OBBox::~OBBox(){};

/*! Copy Constructor
 *\param[in] other OBBox where copy from
 */
OBBox::OBBox(const OBBox & other):BaseManipulation(other){
    m_origin = other.m_origin;
    m_span = other.m_span;
    m_axes = other.m_axes;
    m_listgeo = other.m_listgeo;
    m_forceAABB = other.m_forceAABB;
    m_writeInfo = other.m_writeInfo;
};

/*!
 * Swap function
 * \param[in] x object to be swapped
 */
void OBBox::swap(OBBox & x) noexcept
{
    std::swap(m_origin, x.m_origin);
    std::swap(m_span, x.m_span);
    std::swap(m_axes, x.m_axes);
    std::swap(m_listgeo, x.m_listgeo);
    std::swap(m_forceAABB, x.m_forceAABB);
    std::swap(m_writeInfo, x.m_writeInfo);

    BaseManipulation::swap(x);
}
/*! It builds the input/output ports of the object
 */
void
OBBox::buildPorts(){

    bool built = true;
    built = (built && createPortIn<MimmoObject*, OBBox>(this, &mimmo::OBBox::setGeometry, M_GEOM,true,1));
    built = (built && createPortIn<std::vector<MimmoObject*>, OBBox>(this, &mimmo::OBBox::setGeometries, M_VECGEOM, true,1));
    built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getOrigin, M_POINT));
    built = (built && createPortOut<dmatrix33E, OBBox>(this, &mimmo::OBBox::getAxes, M_AXES));
    built = (built && createPortOut<darray3E, OBBox>(this, &mimmo::OBBox::getSpan, M_SPAN));

    m_arePortsBuilt = built;
};

/*!Clean all stuffs in your class */
void
OBBox::clearOBBox(){
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
    m_forceAABB = false;
    m_writeInfo = false;
};

/*!
 * Return the origin of the OBB.
 * \return Number of control nodes
 */
darray3E
OBBox::getOrigin(){
    return(m_origin);
}

/*!
 * Return the span of the OBB.
 * \return Number of control nodes
 */
darray3E
OBBox::getSpan(){
    return(m_span);
}


/*!
 * Return the oriented axes of the OBB.
 * \return Number of control nodes
 */
dmatrix33E
OBBox::getAxes(){

    return(m_axes);
}

/*!
 * \return current list of linnked geoemtry for obbox calculation
 */
std::vector<MimmoObject*>
OBBox::getGeometries(){
    std::vector<MimmoObject*> result;
    result.reserve(m_listgeo.size());
    for(auto pp: m_listgeo)
        result.push_back(pp.first);
    return result;
};

/*!
 * \return true if the class is set to evaluate a standard Axis Aligned Bounding Box
 */
bool
OBBox::isForcedAABB(){
    return m_forceAABB;
}

/*!
 * Set the list of target geometries once and for all, and erase any pre-existent list.
 * Not supported volumetric tessellations(type =2).
 * \param[in] listgeo list of pointers to target geometries
 */
void
OBBox::setGeometries(std::vector<MimmoObject *> listgeo){
    m_listgeo.clear();
    for( auto val : listgeo){
        setGeometry(val);
    }
};

/*!
 * * Add your target geometry to the list of target geometries. Not supported volumetric tessellations(type =2).
 * Reimplemented from BaseManipulation::setGeometry().
 * \param[in] geo pointer to target geometry
 */
void
OBBox::setGeometry(MimmoObject * geo){

    if (geo == NULL)    {
        (*m_log)<<"warning: "<<m_name<<" not valid Geometry pointer. Doing nothing"<<std::endl;
        return;
    }

    if (geo->isEmpty())    {
        (*m_log)<<"warning: "<<m_name<<" empty Geometry set. Doing nothing"<<std::endl;
        return;
    }

    if (geo->getType() == 2 )    {
        (*m_log)<<"warning: "<<m_name<<" does not support volumetric tessellation. Geometry not set"<<std::endl;
        return;
    }
    m_listgeo.insert(std::make_pair(geo, geo->getType()));
};

/*!
 * Force class to evaluate a standard Axis Aligned Bounding Box instead of the Oriented one.
 * \param[in] flag true, force AABB evaluation.
 */
void
OBBox::setForceAABB(bool flag){
    m_forceAABB = flag;
}

/*!
 * Force class to write OBB calculated origin, axes and span to file(in the plotOptionalResults directory)
 * \param[in] flag true, write the OBB info on file.
 */
void
OBBox::setWriteInfo(bool flag){
    m_writeInfo = flag;
}

/*! Plot the OBB as a structured grid to *vtu file.
 * \param[in] directory output directory
 * \param[in] filename  output filename w/out tag
 * \param[in] counter   integer identifier of the file
 * \param[in] binary    boolean flag for 0-"ascii" or 1-"appended" writing
 */
void
OBBox::plot(std::string directory, std::string filename,int counter, bool binary){


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


/*!Execute your object, calculate the OBBox of your geometry.
 * If forced externally, evaluate the AABB, no matter what.
 * Implementation of pure virtual BaseManipulation::execute
 */
void
OBBox::execute(){

    darray3E pmin, pmax;
    pmin.fill(1.e18);
    pmax.fill(-1.e18);
    double val;
    dmatrix33E eye;
    {
        int count = 0;
        for(auto & local: eye) {
            local.fill(0);
            local[count] = 1.0;
            m_axes[count] = local;
            ++count;
        }
    }
    // check if the list is empty.
    if(m_listgeo.empty()){
        m_span.fill(0.0);
        m_origin.fill(0.0);
        return;
    };

    std::unordered_map<MimmoObject*, int>::iterator itB = m_listgeo.begin();
    //if one geometry at least is a cloud point, solve all them as cloud points.
    bool allCloud = false;
    while(itB != m_listgeo.end() && !allCloud){
        allCloud = (itB->second == 3 || itB->second == 4);
        itB++;
    }

    if(!m_forceAABB){

        dmatrix33E covariance;
        darray3E spectrum;

        if(allCloud){
            covariance = evaluatePointsCovarianceMatrix(getGeometries());
        }else{
            covariance = evaluateElementsCovarianceMatrix(getGeometries());
        }
        m_axes = eigenVectors(covariance, spectrum);
        adjustBasis(m_axes, spectrum);
    }


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
            out<<"OBBox "<<std::to_string(getId())<<" info:"<<std::endl;
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
OBBox::plotOptionalResults(){
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
OBBox::absorbSectionXML(const bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::absorbSectionXML(slotXML, name);

    if(slotXML.hasOption("ForceAABB")){
        std::string input = slotXML.get("ForceAABB");
        input = bitpit::utils::string::trim(input);
        bool value = false;
        if(!input.empty()){
            std::stringstream ss(input);
            ss >> value;
        }
        setForceAABB(value);
    }
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

}

/*!
 * It sets infos from class members in a XML bitpit::Config::section.
 * \param[in] slotXML bitpit::Config::Section of XML file
 * \param[in] name   name associated to the slot
 */
void
OBBox::flushSectionXML(bitpit::Config::Section & slotXML, std::string name){

    BITPIT_UNUSED(name);
    BaseManipulation::flushSectionXML(slotXML, name);

    slotXML.set("ForceAABB", std::to_string((int)m_forceAABB));
    slotXML.set("WriteInfo", std::to_string((int)m_writeInfo));

};


/*!
 * Assembly covariance matrix of various target geometries treated as cloud points.
 * \param[in] list    list of geometries
 * \return covariance matrix;
 */
dmatrix33E
OBBox::evaluatePointsCovarianceMatrix(std::vector<MimmoObject*> list){

    //evaluate the mass center first;
    darray3E masscenter = {{0.0,0.0,0.0}};
    int countVert = 0;
    for(auto geo: list){
        for(auto & vert: geo->getVertices()){
            masscenter += vert.getCoords();
            ++countVert;
        }
    }

    masscenter /= double(countVert);

    //assembly covariance matrix node by node.
    dmatrix33E loc;
    dmatrix33E covariance;
    for(auto & val:covariance)    val.fill(0.0);

    darray3E temp;
    for(auto geo: list){
        for(auto & val:loc)    val.fill(0.0);
        //evaluate local contributes to covariance
        for(auto & vert: geo->getVertices()){
            temp = vert.getCoords() - masscenter;
            for(int j=0; j<3; ++j){
                for(int k=j; k<3; ++k){
                    loc[j][k] += temp[j]*temp[k];
                }
            }
        }
        //store temporarely in covContributes
        int counter = 0;
        for(auto & val: covariance){
            val += loc[counter];
            ++counter;
        }
    }

    int counter = 0;
    for(auto & val: covariance){
        val /= double(countVert);
        ++counter;
    }

    covariance[1][0] = covariance[0][1];
    covariance[2][0] = covariance[0][2];
    covariance[2][1] = covariance[1][2];

    return covariance;
};


/*!
 * Assembly covariance matrix of various target surface geometries treated as tesselations.
 * Here are excluded 3DCurve, Point Clouds and Volume meshes
 * \param[in] list    list of geometries
 * \return covariance matrix;
 */
dmatrix33E
OBBox::evaluateElementsCovarianceMatrix(std::vector<MimmoObject*> list){

    //You need to evaluate the mass center and the 3x3 matrix of the total moments
    //obtained as the sum of each element moments.
    // Non triangular cells are approximated as a representative triangle formed by
    // the first set of 3 non-aligned vertices.


    darray3E masscenter = {{0.0,0.0,0.0}};
    //assembly moments matrix cell by cell.
    dmatrix33E moments;
    for(auto & val:moments)    val.fill(0.0);

    darray3E p,q,r, centroid;
    double areatri, areatot(0.0);

    for(auto geo: list){
        for(auto & cell: geo->getCells()){

            std::array<long,3> vids = get3RepPoints(cell.getId(), geo);
            if(vids[0] < 0){
                continue;
            }
            p = geo->getVertexCoords(vids[0]);
            q = geo->getVertexCoords(vids[1]);
            r = geo->getVertexCoords(vids[2]);

            areatri = 0.5*norm2(crossProduct(q-p, r-p));
            centroid = (p+q+r)/3.0;
            masscenter += areatri*centroid;
            areatot += areatri;

            //assemby moments
            // on-diagonal terms
            moments[0][0] += areatri*(9.0*centroid[0]*centroid[0] + p[0]*p[0] + q[0]*q[0] + r[0]*r[0])/12.0;
            moments[1][1] += areatri*(9.0*centroid[1]*centroid[1] + p[1]*p[1] + q[1]*q[1] + r[1]*r[1])/12.0;
            moments[2][2] += areatri*(9.0*centroid[2]*centroid[2] + p[2]*p[2] + q[2]*q[2] + r[2]*r[2])/12.0;

            // off-diagonal terms
            moments[0][1] += areatri*(9.0*centroid[0]*centroid[1] + p[0]*p[1] + q[0]*q[1] + r[0]*r[1])/12.0;
            moments[0][2] += areatri*(9.0*centroid[0]*centroid[2] + p[0]*p[2] + q[0]*q[2] + r[0]*r[2])/12.0;
            moments[1][2] += areatri*(9.0*centroid[1]*centroid[2] + p[1]*p[2] + q[1]*q[2] + r[1]*r[2])/12.0;
        }
    }// end of cell by cell loop.

    //normalize masscenter
    masscenter /= areatot;

    // moments are symmetric
    moments[1][0] = moments[0][1];
    moments[2][0] = moments[0][2];
    moments[2][1] = moments[1][2];

    //get final covariance matrix
    dmatrix33E covariance;
    for(int i=0; i<3; ++i ){
        for(int j=0; j<3; ++j ){
            covariance[i][j] = moments[i][j]/areatot - masscenter[i]*masscenter[j];
        }
    }
    return covariance;
};

/*!
 * Calculates and returns the eigenVectors and eigenvalues of a 3x3 matrix.
 * eigvec and eigenval are in descending order (eigvec from max to min eigenval they refers)
 * \param[in] matrix    target matrix
 * \param[out]    eigenvalues eigenvalues of the matrix
 * \return    matrix of eigenvectors by column
 */
dmatrix33E
OBBox::eigenVectors( dmatrix33E & matrix, darray3E & eigenvalues){

    dmatrix33E result;
    double * a = new double [9];
    double * s = new double[3];

    int k=0;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            //transpose
            a[k] = matrix[j][i];
            k++;
        }
    }

    LAPACKE_dsyev( LAPACK_COL_MAJOR, 'V','U', 3, a, 3, s);
    //eigenvec are returned in a column wise, and with eigenvals s in ascending order.

    //rearrange the eigvec
    for (int i=0; i<9; i++){
        result[i/3][i%3] = a[i];
    }
    for(int i=0; i<3; ++i)    {
        result[i] /= norm2(result[i]);
        eigenvalues[i] = s[i];
    }

    //change them in descending order.
    std::swap(result[0], result[2]);
    std::swap(eigenvalues[0], eigenvalues[2]);

    // clean allocations.
    delete [] a; a = NULL;
    delete [] s; s = NULL;

    return result;
}

/*!
 * Adjust basis of eigenvectors, given is eigenvalues spectrum. In order to best fit the
 * 3D shape.  The rules are the following
 * 1) three real coincident eigenvalues, return fundamental axis as eigenvectors
 * 2) three real different eigenvalues, order in max, mid, min axis
 * 3) two coincident eigenvalues, one not, find best fit to shape for eigVec associated to coincident eigenvalues
 * \param[in,out]    eigVec basis of eigenvectors by rows
 * \param[in]    eigVal eigenvalues associated to eigVec
 */
void
OBBox::adjustBasis(dmatrix33E & eigVec, darray3E & eigVal){

    darray3E diff;
    double tol = std::numeric_limits<double>::min();

    diff[0] = std::abs(eigVal[1]-eigVal[0]);
    diff[1] = std::abs(eigVal[2]-eigVal[1]);
    diff[2] = std::abs(eigVal[0]-eigVal[2]);

    //3 real quasi coincident eigenval-> return fundamental axis.
    if(diff[0]<=tol && diff[1]<=tol && diff[2]<= tol){
        int counter=0;
        for(auto & val : eigVec){
            val.fill(0.0);
            val[counter] = 1.0;
            ++counter;
        }
        return;
    }

    // 3 real distinct eigenval
    if(diff[0]>tol && diff[1]>tol && diff[2]>tol){
        //you need a real reference triad here. Check it out
        if(dotProduct(eigVec[2], crossProduct(eigVec[0], eigVec[1])) < 0.0){
            eigVec[2] *= -1.0;
        };
        return;
    };

    int guess = 0, third =1, stable=2;

    if(diff[1] <= tol)    {
        guess = 1;
        third = 2;
        stable= 0;
    }

    eigVec[third] = crossProduct(eigVec[stable],eigVec[guess]);
    double normv = norm2(eigVec[third]);
    if(normv > std::numeric_limits<double>::min()){
        eigVec[third] /= normv;
    }

    //you need a real reference triad here. Check it out
    if(dotProduct(eigVec[2], crossProduct(eigVec[0], eigVec[1])) < 0.0){
        eigVec[2] *= -1.0;
    };

};


/*!
    Transpose a 3x3 double matrix
    \param[in] mat target matrix
    \return new matrix transposed
*/
dmatrix33E OBBox::transpose(const dmatrix33E & mat){
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
dmatrix33E OBBox::inverse(const dmatrix33E & mat){
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


/*!
    Get first 3, non-aligned representative point ids of a 2D cell of a certain geometry.
    Works only for surface mesh (mimmoObject type 2);
    \param[in] cellID id of the target cell.
    \param[in] geo pointer to reference geometry (must be always of mimmoObject type 1)
    \return ids of vertices composing the representative triangle.Return -1,-1,-1 for unsupported cell elements
*/
std::array<long,3> OBBox::get3RepPoints(long cellID, MimmoObject * geo){
    if(!geo) return{{-1,-1,-1}};
    bitpit::Cell & cell = geo->getCells().at(cellID);
    std::array<long,3> result;

    switch(cell.getType()){
        case bitpit::ElementType::TRIANGLE:
            {
                long * conn = cell.getConnect();
                result[0] = conn[0];
                result[1] = conn[1];
                result[2] = conn[2];
            }
        break;
        case bitpit::ElementType::QUAD:
        case bitpit::ElementType::POLYGON:
            {
                bitpit::ConstProxyVector<long> vertids = cell.getVertexIds();
                std::size_t ssz= vertids.size();

                //find a set of 3 points whose cross product is non zero, if you can.
                darray3E p0,s1,s2;
                result[0] = vertids[0];
                p0 = geo->getVertexCoords(result[0]);

                std::size_t cand = 2;
                bool found = false;

                while(cand < ssz && !found){
                    result[1] = vertids[cand-1];
                    result[2] = vertids[cand];
                    s1 = geo->getVertexCoords(result[1]) - p0;
                    s2 = geo->getVertexCoords(result[2]) - p0;
                    found = (norm2(crossProduct(s1,s2)) > 0.0);
                    ++cand;
                }
            }
        break;
        default:
            result.fill(-1);
        break;
    }

    return result;
}


}
