/*---------------------------------------------------------------------------*\
 *
 *  MiMMO
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of MiMMO.
 *
 *  MiMMO is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  MiMMO is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with MiMMO. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#include "RotationBox.hpp"

///*!Default constructor of RotationBox
// */
RotationBox::RotationBox(darray3E origin, darray3E direction){
	m_ndeg = 1;
	m_origin = origin;
	m_direction = direction;
};
//
///*!Default destructor of RotationBox
// */
RotationBox::~RotationBox(){};

///*!Copy constructor of RotationBox.
// */
RotationBox::RotationBox(const RotationBox & other):BaseManipulation(other){
	m_origin = other.m_origin;
	m_direction = other.m_direction;
};

/*!Assignement operator of RotationBox.
 */
RotationBox & RotationBox::operator=(const RotationBox & other){
	*(static_cast<BaseManipulation*> (this)) = *(static_cast<const BaseManipulation*> (&other));
	m_origin = other.m_origin;
	m_direction = other.m_direction;
	return(*this);
};

void
RotationBox::setAxis(darray3E origin, darray3E direction){
	m_origin = origin;
	m_direction = direction;
}

void
RotationBox::setOrigin(darray3E origin){
	m_origin = origin;
}

void
RotationBox::setDirection(darray3E direction){
	m_direction = direction;
	for (int i=0; i<3; i++)
		m_direction[i] /= norm2(m_direction);
}

void
RotationBox::setRotation(double alpha){
	m_displ.resize(1);
	m_displ[0] = { {alpha, 0 , 0} };
}

void
RotationBox::useInfo(){
	m_axes.resize(m_info->m_naxes);
	for (int i=0; i<m_info->m_naxes; i++){
		for (int j=0; j<m_info->m_naxes; j++){
			m_axes[i][j] = m_info->m_axes[i][j];
		}
	}
}

/*!Execution command. It modifies the coordinates of the origin given by the child manipulation object
 * with the rotation conditions. After exec() the original origin will be permanently modified.
 * Set the translated origin only for one child (the first one) and it has to be a FFDLattice
 * (static cast to use setOrigin method of basic shape).
 */
void
RotationBox::execute(){
	dmatrix44E D, Rx, Ry, Rz, Rym1, Rxm1, Dm1, T;
	dvecarr4E P2(m_axes.size()), P1(m_axes.size());
	double V = sqrt(m_direction[1]*m_direction[1] + m_direction[2]*m_direction[2]);
	double L = sqrt(m_direction[0]*m_direction[0] + m_direction[1]*m_direction[1] + m_direction[2]*m_direction[2]);

	bvector1D perpend(3);

	for (int i=0;i<3; i++){
		if (m_direction[i] == 0){
			perpend[i] = true;
		}
		std::cout << "perpend " << i << " : " << perpend[i] << std::endl;
	}

	std::cout << "V : " << V << std::endl;
	std::cout << "L : " << L << std::endl;


	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			D[i][j] = Rx[i][j] = Ry[i][j] = Rz[i][j] = Rym1[i][j] = Rxm1[i][j] = Dm1[i][j] = T[i][j] = 0.0;
		}
	}
	for (int i=0; i<3; i++){
		D[i][i] = 1.0;
		D[i][3] = -m_origin[i];
		Rx[i][i] = 1.0;
		Ry[i][i] = 1.0;
		Rz[i][i] = 1.0;
		Dm1[i][i] = 1.0;
		Dm1[i][3] = m_origin[i];
		Rxm1[i][i] = 1.0;
		Rym1[i][i] = 1.0;
		for (int j=0; j<m_axes.size(); j++){
			P1[j][i] = m_axes[j][i];
		}
	}
	D[3][3] = 1.0;
	Rx[3][3] = 1.0;
	Ry[3][3] = 1.0;
	Rz[3][3] = 1.0;
	Dm1[3][3] = 1.0;
	Rxm1[3][3] = 1.0;
	Rym1[3][3] = 1.0;
	for (int j=0; j<m_axes.size(); j++){
		P1[j][3] = 1.0;
	}

	Rx[1][1] = m_direction[2]/V;
	Rx[2][1] = m_direction[1]/V;
	Rx[1][2] = -m_direction[1]/V;
	Rx[2][2] = m_direction[2]/V;

	Rxm1[1][1] = m_direction[2]/V;
	Rxm1[2][1] = -m_direction[1]/V;
	Rxm1[1][2] = m_direction[1]/V;
	Rxm1[2][2] = m_direction[2]/V;

	Ry[0][0] = V/L;
	Ry[0][2] = -m_direction[0]/L;
	Ry[2][0] = m_direction[0]/L;
	Ry[2][2] = V/L;

	Rym1[0][0] = V/L;
	Rym1[0][2] = m_direction[0]/L;
	Rym1[2][0] = -m_direction[0]/L;
	Rym1[2][2] = V/L;

	Rz[0][0] = cos(m_displ[0][0]);
	Rz[0][1] = -sin(m_displ[0][0]);
	Rz[1][0] = sin(m_displ[0][0]);
	Rz[1][1] = cos(m_displ[0][0]);


	T = Dm1;
	T = matMul(T,Rxm1);

	T = matMul(T, Rym1);
	T = matMul(T, Rz);
	T = matMul(T, Ry);
	T = matMul(T, Rx);
	T = matMul(T, D);

	std::cout << T << std::endl;
	std::cout << P1 << std::endl;

	for (int i=0; i<m_axes.size(); i++){
		for (int j=0; j<4; j++){
			P2[i][j] += dotProduct(T[j],P1[i]);
		}
	}
	std::cout << "P2 : " << P2 << std::endl;

	dvecarr3E rotated(m_axes.size());

	for (int i=0; i<m_axes.size(); i++){
		for (int j=0; j<3; j++){
			rotated[i][j] = P2[j][i];
		}
	}

	if (m_child[0] != NULL){
		static_cast<FFDLattice*>(m_child[0])->getShape()->setRefSystem(rotated[0], rotated[1], rotated[2]);
	}
	return;
};


dmatrix44E
RotationBox::matMul(dmatrix44E & a, dmatrix44E & b){
	dmatrix44E res;
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			res[i][j] = 0;
			for (int k=0; k<4; k++){
				res[i][j] += a[i][k]*b[k][j];
			}
		}
	}
	return(res);
}

