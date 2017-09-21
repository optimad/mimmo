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
#include "LinearAlgebra.hpp"
#include "lapacke.h"

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
    dmatrix33E    trasp = bitpit::linearalgebra::transpose(m_axes);
    for(auto &val : activeP){


        for(int i=0; i<3; ++i){
            temp[i] = dotProduct(val, trasp[i]);
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
 * If the OBBox is not fit enough, return to its AABB version.
 * If forced externally, evaluate the AABB, no matter what. 
 * Implementation of pure virtual BaseManipulation::execute
 */
void
OBBox::execute(){
    if(m_listgeo.empty()){
        throw std::runtime_error(m_name + " : found empty list of geometries");
    };

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

    std::unordered_map<MimmoObject*, int>::iterator itB = m_listgeo.begin();
    while(itB != m_listgeo.end() && itB->second != 3){
        itB++;
    }
    bool allCloud = (itB != m_listgeo.end());

    if(!m_forceAABB){

        dmatrix33E covariance, covContributes;
        darray3E spectrum;
        darray3E etaPoint;


        etaPoint = evaluateMassCenter(getGeometries(), allCloud);
        assemblyCovContributes(getGeometries(), allCloud, covContributes); 
        covariance = evaluateCovarianceMatrix(covContributes, etaPoint);

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

    dmatrix33E trasp = bitpit::linearalgebra::transpose(m_axes);
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
        m_origin[i] = dotProduct(originLoc, trasp[i]);
    }

    double volOBB = m_span[0]*m_span[1]*m_span[2];

    if(!m_forceAABB){ //check if AABB get a better result
        pmin.fill(1.e18);
        pmax.fill(-1.e18);

        for(auto ptr : getGeometries()){
            for(auto & vert: ptr->getVertices()){
                darray3E coord = vert.getCoords();
                for(int i=0;i<3; ++i){
                    pmin[i] = std::fmin(pmin[i], coord[i]);
                    pmax[i] = std::fmax(pmax[i], coord[i]);
                }
            }
        }

        darray3E span2 = pmax - pmin;
        darray3E orig = 0.5*(pmin+pmax);
        //check if one of the span goes to 0;
        avg_span = 0.0;
        for(auto & val: span2)    avg_span+=val;
        avg_span /= 3.0;

        for(auto &val : span2)    {
            val = std::fmax(val, 1.E-04*avg_span);
        }

        double volAAA =span2[0]*span2[1]*span2[2];


        if(volAAA <= volOBB){
            for(int i=0; i<3; ++i)    m_axes[i] = eye[i];
            m_span = span2;
            m_origin = orig;
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
    plot(dir, nameGrid, getClassCounter(), true );
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

};

/*!
 *\return assembled covariance matrix, form collection of covariance contributes of all geometries and
 *estimation of t heri global center mass are provided.
 *\param[in] covContributes covariance contributes from all geometries linked
 *\param[in] centermass global centermass of group of geometries   
 * 
 */
dmatrix33E 
OBBox::evaluateCovarianceMatrix(dmatrix33E & covContributes, darray3E & centermass){
    dmatrix33E result;

    for(auto & val: result) val.fill(0);

    for(int i=0; i<3; ++i){
        for(int j=i; j<3; ++j){
            result[i][j] = covContributes[i][j] - centermass[i]*centermass[j];
        }
    }

    result[1][0] = result[0][1];
    result[2][0] = result[0][2];
    result[2][1] = result[1][2];
    return result;
}


/*! 
 * Assembly covariance contributes of various target geometries, that can be used to evaluate covariance matrix
 * The method intrinsecally distinguish between cloud point and tessellation, according to target MimmoObject geometry type.
 * \param[out] covContributes matrix to store all covariance contributes
 * \param[in] list    list of geometries
 * \param[in] flag boolean, if true, force all geometries to be treated as point clouds
 */
void
OBBox::assemblyCovContributes(std::vector<MimmoObject*> list, bool flag, dmatrix33E & covContributes){

    dmatrix33E loc;
    for(auto & val:covContributes)    val.fill(0.0);
    int size, counter;
    dvector1D area;
    darray3E temp;
    double areaTot = 0.0;

    if(flag){
        for(auto geo: list){
            for(auto & val:loc)    val.fill(0.0);
            size = geo->getNVertex();

            //evaluate local contributes to covariance
            for(auto & vert: geo->getVertices()){
                temp = vert.getCoords();
                for(int j=0; j<3; ++j){
                    for(int k=j; k<3; ++k){
                        loc[j][k] += temp[j]*temp[k];
                    }
                }
            }
            //store temporarely in covContributes
            counter = 0;
            for(auto & val: covContributes){
                val += loc[counter];
                ++counter;
            }

            areaTot += double(size);
        }
    }else{
        for(auto geo: list){
            for(auto & val:loc)    val.fill(0.0);
            size = geo->getNCells();
            bitpit::SurfUnstructured * tri = static_cast<bitpit::SurfUnstructured * >(geo->getPatch());
            area.resize(size);
            counter = 0;
            for(auto &cell : tri->getCells()){
                area[counter] = tri->evalCellArea(cell.getId());
                areaTot +=area[counter];
                ++counter;
            }

            counter = 0;
            darray3E centroid;
            for(auto & cell: tri->getCells()){
                centroid = tri->evalCellCentroid(cell.getId());

                for(int j=0; j<3; ++j){
                    for(int k=j; k<3; ++k){
                        loc[j][k] += area[counter]*centroid[j]*centroid[k];
                    }
                }
                ++counter;
            }

            //store temporarely in covContributes
            counter = 0;
            for(auto & val: covContributes){
                val += loc[counter];
                ++counter;
            }
        }
    }

    covContributes[1][0] = covContributes[0][1];
    covContributes[2][0] = covContributes[0][2];
    covContributes[2][1] = covContributes[1][2];

    counter = 0;
    for(auto & val: covContributes){
        val /= areaTot;
        ++counter;
    }
};


/*!
 * \return global mass center of the whole list of geometries.
 * \param[in] list of target geometries
 * \param[in] flag boolean, if true, force all geometries to be treated as point clouds
 */
darray3E
OBBox::evaluateMassCenter(std::vector<MimmoObject *> list, bool flag){

    long size = 0;
    darray3E centermass = {{0,0,0}}, eta;
    double areaTot = 0.0;
    int counter;

    if(flag){
        for(auto geo : list){
            size = geo->getNVertex();
            eta.fill(0.0);
            //evaluate eta;
            for(auto & vert: geo->getVertices()){
                 eta += vert.getCoords();
            }
            centermass += eta;
            areaTot += double(size);
        }
    }else{
        for(auto geo : list){
            size = geo->getNCells();
            bitpit::SurfUnstructured * tri = static_cast<bitpit::SurfUnstructured * >(geo->getPatch());
            dvector1D area(size);
            counter = 0;
            for(auto &cell : tri->getCells()){
                area[counter] = tri->evalCellArea(cell.getId());
                areaTot +=area[counter];
                ++counter;
            }

            //evaluate eta;
            eta.fill(0.0);
            counter = 0;
            for(auto & cell: tri->getCells()){
                eta += tri->evalCellCentroid(cell.getId())*area[counter];
                ++counter;
            }
            centermass += eta;
        }
    }

    centermass /=areaTot;
    return centermass;
};

/*!
 * Calculates and returns the eigenVectors and eigenvalues of a 3x3 matrix.
 * \param[in] matrix    target matrix
 * \param[out]    eigenvalues eigenvalues of the matrix
 * \return    matrix of eigenvectors by column
 */
dmatrix33E
OBBox::eigenVectors( dmatrix33E & matrix, darray3E & eigenvalues){

    dmatrix33E result;
    double * a = new double [9];
    double * u = new double [9];
    double * vt = new double [9];
    double * s = new double[3];
    double * superb = new double[3];

    int k=0;
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            //transpose
            a[k] = matrix[j][i];
            k++;
        }
    }

    LAPACKE_dgesvd( LAPACK_COL_MAJOR, 'N','S', 3, 3, a, 3, s, u, 3, vt, 3, superb);

    //solution norm
    for (int i=0; i<9; i++){
        result[i/3][i%3] = vt[i];
    }
    for(int i=0; i<2; ++i)    {
        result[i] /= norm2(result[i]);
        eigenvalues[i] = s[i];
    }
    result[2] = crossProduct(result[0],result[1]);
    eigenvalues[2] = s[2];

    delete [] a; a = NULL;
    delete [] u; u = NULL;
    delete [] vt; vt = NULL;
    delete [] s; s = NULL;
    delete [] superb; superb = NULL;

    return result;
}

/*!
 * Adjust basis of eigenvectors, given is eigenvalues spectrum. In order to best fit the
 * 3D shape.  The rules are the following
 * 1) three real coincident eigenvalues, return fundamental axis as eigenvectors
 * 2) three real different eigenvalues, do nothing
 * 3) two coincident eigenvalues, one not, find best fit to shape for eigVec associated to coincident eigenvalues.
 * \param[in,out]    eigVec basis of eigenvectors by rows
 * \param[in]    eigVal eigenvalues associated to eigVec
 */
void
OBBox::adjustBasis(dmatrix33E & eigVec, darray3E & eigVal){

    darray3E diff;
    double tol = 1.0E-3;

    diff[0] = std::abs(eigVal[1]-eigVal[0]);
    diff[1] = std::abs(eigVal[2]-eigVal[1]);
    diff[2] = std::abs(eigVal[0]-eigVal[2]);

    if(diff[0]>tol && diff[1]>tol && diff[2]>tol)    return;
    if(diff[0]<=tol && diff[1]<=tol && diff[2]<= tol){
        int counter=0;
        for(auto & val : eigVec){
            val.fill(0.0);
            val[counter] = 1.0;
            ++counter;
        }
        return;
    }

    int guess = 0, third =1, stable=2;

    if(diff[1] <= tol)    {
        guess = 1;
        third = 2;
        stable= 0;
    }

    if(diff[2] <= tol)    {
        guess = 2;
        third = 0;
        stable= 1;
    }

    dmatrix33E axes; //, trasp;
    int counter=0;
    for(const auto & val: eigVec){
        axes[counter] = val;
        ++counter;
    }

    darray3E    guessVec = axes[guess];
    darray3E    refVec   = axes[stable];

    int nstage = 15;
    int niterations = 8;

    double distance = M_PI/(2.0*(nstage));
    double volume, theta, val;
    std::map<double, double>    mapval;
    darray3E pmin, pmax, temp;

    for(int i=0; i<nstage; ++i){

        theta = i*distance;

        axes[guess] = guessVec*std::cos(theta) + std::sin(theta)*crossProduct(refVec, guessVec);
        axes[third] = crossProduct(axes[stable],axes[guess]);

        pmin.fill(1.e18);
        pmax.fill(-1.e18);

        for(auto geo : getGeometries()){
            for(auto & vert: geo->getVertices()){
                darray3E coord = vert.getCoords();
                for(int i=0;i<3; ++i){
                    val = dotProduct(coord, axes[i]); //trasp[i]);
                    pmin[i] = std::fmin(pmin[i], val);
                    pmax[i] = std::fmax(pmax[i], val);
                }
            }
        }

        temp = pmax - pmin;
        volume = temp[0]*temp[1]*temp[2];

        mapval[volume] = theta;
    }

    int it =0;
    while (it < niterations){


        distance /= 2.0;

        double thetadx = mapval.begin()->second + distance;
        if(thetadx >= M_PI/2.0)    thetadx += -1.0*M_PI/2.0;

        double thetasx = mapval.begin()->second - distance;
        if(thetasx <= 0.0)    thetasx += M_PI/2.0;

        //evaluate on the right
        axes[guess] = guessVec*std::cos(thetadx) + std::sin(thetadx)*crossProduct(refVec, guessVec);
        axes[third] = crossProduct(axes[stable],axes[guess]);


        //trasp = bitpit::linearalgebra::transpose(axes);
        pmin.fill(1.e18);
        pmax.fill(-1.e18);
        for(auto geo : getGeometries()){
            for(auto & vert: geo->getVertices()){
                darray3E coord = vert.getCoords();
                for(int i=0;i<3; ++i){
                    val = dotProduct(coord, axes[i]); //trasp[i]);
                    pmin[i] = std::fmin(pmin[i], val);
                    pmax[i] = std::fmax(pmax[i], val);
                }
            }
        }

        temp = pmax - pmin;
        volume = temp[0]*temp[1]*temp[2];

        mapval[volume] = thetadx;

        //evaluate on the left
        axes[guess] = guessVec*std::cos(thetasx) + std::sin(thetasx)*crossProduct(refVec, guessVec);
        axes[third] = crossProduct(axes[stable],axes[guess]);


        pmin.fill(1.e18);
        pmax.fill(-1.e18);
        for(auto geo : getGeometries()){
            for(auto & vert: geo->getVertices()){
                darray3E coord = vert.getCoords();
                for(int i=0;i<3; ++i){
                    val = dotProduct(coord, axes[i]); //trasp[i]);
                    pmin[i] = std::fmin(pmin[i], val);
                    pmax[i] = std::fmax(pmax[i], val);
                }
            }
        }

        temp = pmax - pmin;
        volume = temp[0]*temp[1]*temp[2];

        mapval[volume] = thetasx;
        ++it;
    }

    theta = mapval.begin()->second;
    eigVec[guess] = guessVec*std::cos(theta) + std::sin(theta)*crossProduct(refVec, guessVec);
    eigVec[third] = crossProduct(eigVec[stable],eigVec[guess]);

};

}
