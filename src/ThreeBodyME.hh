///////////////////////////////////////////////////////////////////////////////////
//    ThreeBodyME.hh, part of  imsrg++
//    Copyright (C) 2018  Ragnar Stroberg
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
///////////////////////////////////////////////////////////////////////////////////

#ifndef ThreeBodyME_h
#define ThreeBodyME_h 1

#include "ModelSpace.hh"
#include <fstream>
#include <unordered_map>

//typedef double ThreeBME_type;
//typedef float ThreeBME_type;

/// The three-body piece of an operator, stored in nested vectors.
/// The 3BMEs are stored in unnormalized JT coupled form
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f \rangle_{JT} \f$.
/// To minimize the number of stored matrix elements, only elements with
/// \f$ a\geq b \geq c, a\geq d\geq e \geq f \f$ are stored.
/// The other combinations are obtained on the fly by GetME().
/// The storage format is MatEl[{a,b,c,d,e,f,J,Jab,Jde}][T_index] =
/// \f$ \langle (abJ_{ab}t_{ab})c | V | (deJ_{de}t_{de})f  \rangle_{JT} \f$.
class ThreeBodyME
{
 public:
  typedef float ME_type;
  ModelSpace * modelspace;
  std::vector<ME_type> MatEl;
  std::unordered_map<size_t, size_t> OrbitIndexHash; 
  int E3max;
  int emax; // usually, this should be the emax of the modelspace, but we might want something smaller.
  int herm; // +1 for hermitian, -1 for anti-hermitian
  size_t total_dimension;

  int rank_J=0;
  int rank_T=0;
  int parity=0;
  int ISOSPIN_BLOCK_DIMENSION=5;

  bool is_allocated = false;

  const static int ABC;
  const static int BCA;
  const static int CAB;
  const static int ACB;
  const static int BAC;
  const static int CBA;
  
  ~ThreeBodyME();
  ThreeBodyME();
  ThreeBodyME(ModelSpace*);
  ThreeBodyME(ModelSpace*, int rank_J, int rank_T, int parity);
  ThreeBodyME(ModelSpace* ms, int e3max);
  ThreeBodyME(ModelSpace* ms, int e3max, int rank_J, int rank_T, int parity);
//  ThreeBodyME(ThreeBodyME tbme);
  ThreeBodyME(const ThreeBodyME& Tbme);


  ThreeBodyME& operator*=(const double);
  ThreeBodyME& operator+=(const ThreeBodyME&);
  ThreeBodyME& operator-=(const ThreeBodyME&);

  size_t KeyHash(size_t,size_t,size_t,size_t,size_t,size_t) const;
  void KeyUnhash(size_t& key, size_t& a, size_t& b, size_t& c, size_t& d, size_t& e, size_t& f) const;
  void Allocate();

  void SetModelSpace(ModelSpace *ms){modelspace = ms;};
  ModelSpace* GetModelSpace(){return modelspace;};

//// Three body setter getters
//  std::vector<std::pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n) const;
  void AddToME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ME_type V);
  void   SetME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ME_type V);
  ME_type GetME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n) const;
//  void AddToME_pn(int Jab_in, int Jde_in, int twoJ, int i, int j, int k, int l, int m, int n, ME_type V);
  void   SetME_pn(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoT, int i, int j, int k, int l, int m, int n, ME_type V);
  ME_type GetME_pn(int Jab_in, int Jde_in, int twoJ, int i, int j, int k, int l, int m, int n) const;


  std::vector<std::pair<size_t,double>> AccessME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;
  void AddToME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V);
  void   SetME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f, ME_type V);
  ME_type GetME(int Jab_in, int Jde_in, int twoJ, int tab_in, int tde_in, int twoTabc, int twoTdef, int a, int b, int c, int d, int e, int f) const;

///// Some other three body methods

  int SortOrbits(int a_in, int b_in, int c_in, int& a,int& b,int& c) const;
  double RecouplingCoefficient(int recoupling_case, double ja, double jb, double jc, int Jab_in, int Jab, int twoJ) const;
  void SetE3max(int e){E3max = e;};
  int GetE3max(){return E3max;};
  int Getemax(){return emax;};
  void Setemax(int e){emax=  e;};
  void SetHermitian(){herm = +1;};
  void SetAntiHermitian(){herm = -1;};

  double Norm() const;

  void Erase(); // set all three-body terms to zero
  void Deallocate();
  size_t size(){return total_dimension * sizeof(ME_type);};


  void WriteBinary(std::ofstream&);
  void ReadBinary(std::ifstream&);

};


#endif
