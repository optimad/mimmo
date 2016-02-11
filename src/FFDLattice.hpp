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
#ifndef __FFDLATTICE_HPP__
#define __FFDLATTICE_HPP__

/*!
 *	\date			09/feb/2016
 *	\authors		Rocco Arpa
 *	\authors		Edoardo Lombardi
 *
 *	\brief Free Form Deformation of a 3D surface and cloud points.
 *
 *	Bla Bla
 *
 */
class FFDLattice: public BaseManipulation{

private:
	//members
	dvecarr3E		m_result;

public:
	FFDLattice();
	FFDLattice(MimmoObject* geometry);
	~FFDLattice();

	//internal methods
	int 		accessPoint(int i, int j, int k);
	ivector1D 	accessPoint(int idx);
	int 		accessCell(int i, int j, int k);
	ivector1D	accessCell(int idx);
	ivector1D 	getNeighs(int i, int j, int k);
	ivector1D 	getCellOwner(darray3E & coord);

	//TODO VERIFICARE METODI PER KNOTS E GRADO POLINOMI
	dvecarr3E	getResult();
	ivector1D	getDimension();
	ivector1D 	getKnotsDimension();
	darray3E	getOrigin();
	darray3E 	getSpan();
	int  		getKnotInterval(double, int);
	double 		getKnotValue(int, int);
	int 		getKnotIndex(int,int);

	void		setDisplacements(dvecarr3E & displ);
	void 		setGeometry(MimmoObject* geometry);
	void		setDimension(ivector1D dimension);
	void		setDimension(ivector1D dimension, ivector1D knotsDimension);
	void		setKnotsDimension(ivector1D knotsDimension);
	void		setOrigin(darray3E origin);
	void		setSpan(darray3E span);

	void		plotGrid(std::string directory, std::string filename, bool ascii, bool deformed);
	void		plotCloud(std::string directory, std::string filename, bool ascii, bool deformed);

	//relationship methods
protected:
	void 		recoverInfo();

public:
	void 		exec();
	darray3E 	apply(darray3E & point);
	dvecarr3E 	apply(dvecarr3E & point);

};

#endif /* __BASEMANIPULATION_HPP__ */
