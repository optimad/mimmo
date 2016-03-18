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
#ifndef __MIMMONAMESPACE_HPP__
#define __MIMMONAMESPACE_HPP__

namespace mimmo{

namespace pin{


template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL& (G::*fget) (), void (S::*fset) (VAL*));

template<typename OO, typename G, typename OI, typename S, typename VAL>
void addPin(OO* objSend, OI* objRec, VAL* (G::*fget) (), void (S::*fset) (VAL*));

template<typename T, typename U, typename VAL>
std::function<VAL(void)> pinGet(VAL (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<VAL&(void)> pinGetR(VAL& (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<VAL*(void)> pinGetR(VAL* (T::*fget) (), U* obj);

template<typename T, typename U, typename VAL>
std::function<void(VAL)> pinSet(void (T::*fset) (VAL), U* obj);

template<typename T, typename U, typename VAL>
std::function<void(VAL)> pinSet(void (T::*fset) (VAL*), U* obj);

//TODO Implementation of remove pin methods
template<typename OO, typename G, typename OI, typename S, typename VAL>
void removeAllPins(OO* objSend, OI* objRec);

template<typename OO, typename G, typename OI, typename S, typename VAL>
void removePin(OO* objSend, OI* objRec, VAL (G::*fget) (), void (S::*fset) (VAL));


};



};

#include "MimmoNamespace.tpp"





#endif
