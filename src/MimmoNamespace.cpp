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
#include "MimmoNamespace.hpp"

using namespace std;
using namespace mimmo;

/*!It works for subpatch.
 *
 */
map<int,int>
mimmo::mapVertex(MimmoObject* obj1, MimmoObject* obj2){

	map<int,int> mapvertex;

	dvecarr3E vs1 = obj1->getVertex();
	dvecarr3E vs2 = obj2->getVertex();
	ivector1D conn1, conn2;

	bool check;
	darray3E v1, v2;
	int i1, i2, nv1, nv2;
	nv1 = obj1->getNVertex();
	nv2 = obj2->getNVertex();

	for (i1=0; i1<nv1; i1++){
		check = false;
		i2 = 0;
		while (!check && i2 < nv2){
			if (norm2(vs1[i1]-vs2[i2]) <= 1.0e-05){
				check = true;
				mapvertex[i1] = i2;
			}
			i2++;
		}
	}
	return mapvertex;
}

