
#ifndef TMtrans_h
#define TMtrans_h 1

#include "ModelSpace.hh"
#include "Operator"


class RelCMChannel
{
  public:
    int J;
    int Parity;
    int Tz;
    int Nstates;
    std::vector<std::array<int,6>> NLnlSJrel;
    std::map<std::array<int,6>,int> IndexList;
    ~RelCMChannel();
    RelCMChannel(int j, int p, int tz, int nmax);
    int GetNumberStates() const {return Nstates;};
    int GetRelCMStateIndex(int N, int L, int n, int l, int S, int Jrel) const {return IndexList[{N,L,n,l,S,Jrel}];};
}

class OpRelCM
{
  public:
    int rank_J;
    int rank_P;
    int rank_Tz;
    ModelSpace * modelspace;
    std::map<std::array<int,3>,int> IndexList;
    std::vector<RelCMChannel> RelCMChannels;
    std::map<int,arma::mat> Tcoefs;
    std::map<std::array<int,2>,arma::mat> MatEl;

    ~OpRelCM();
    OpRelCM(ModelSpace&, int rankJ, int rankP, int rankTz);
    RelCMChannel & GetRelCMChannel(int ch) const {return (RelCMChannel &) RelCMChannels[ch];};
    arma::mat SetTcoefsChannel(TwoBodyChannel &, RelCmChannel &);
    double Tcoef(int N, int L, int n, int l, int S, int Jrel, int na, int la, int ja, int nb, int lb, int jb, int J, int Tz);
    void TransLabToRelCM(Operator &);
    Operator TransRelCMToLab(ModelSpace &);
    int GetRelCMChannelIndex(int j, int p, int tz) const {return IndexList[{j,p,tz}];};
    double GetME(int ch_bra, int ch_ket, int bra, int ket) const {return MatEl[{ch_bra,ch_ket}](bra,ket);};
    double GetME(int Nbra, int Lbra, int nbra, int lbra, int Sbra, int Jrbra, int Jbra, int Tzbra,
        int Nket, int Lket, int nket, int lket, int Sket, int Jrket, int Jket, int Tzket);
    double GetME(int Nbra, int Lbra, int nbra, int lbra, int Sbra, int Jrbra,
        int Nket, int Lket, int nket, int lket, int Sket, int Jrket, int J, int Tz);
}

#endif
