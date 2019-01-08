#ifndef NaturalOrbital_h
#define NaturalOrbital_h

#include "HartreeFock.hh"
#include <armadillo>
#include <map>
#define OCC_CUT 1e-6

typedef unsigned long long int index_t;
class HFMBPT : public HartreeFock
{
  public:
    arma::mat C_HO2NAT; // transforamtion coefficients, 1st index ho, 2nd index NAT basis
    arma::mat C_HF2NAT; // transforamtion coefficients, 1st index hf, 2nd index NAT basis
    arma::vec Occ; // Occupation number
    vector<index_t> holes;     // in the reference Slater determinant
    vector<index_t> particles; // above the reference Slater determinant

    ~HFMBPT();
    HFMBPT(Operator& hbare); // same as HartreeFock constructor
    void GetNaturalOrbital(string mode);
    void GetDensityMatrix();
    void DensityMatrixPP(Operator& H);
    void DensityMatrixHH(Operator& H);
    void DensityMatrixPH(Operator& H);
    void DiagRho();
    void PrintOccupation();
    Operator TransformHFToNATBasis(Operator& OpIn);
    Operator TransformHOToNATBasis(Operator& OpIn);
    Operator GetNormalOrderedHNAT();
};
#endif
