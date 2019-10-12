#include "OpRelCM.hh"
RelCMChannel::~RelCMChannel()
{}

RelCMChannel::RelCMChannel(int j, int p, int tz, int nmax)
  : J(j), Parity(p), Tz(tz)
{
  for (int L=0; L<=nmax; L++){
    for (int l=0; l<=(nmax-L); l++){
      if((l+L)%2 != Parity) continue;
      for (int S : {0,1}){
        for (int Jrel=std::abs(l-S); Jrel<=(l+S); Jrel++){
          if(std::abs(Jrel-L) > j) continue;
          if(Jrel+L < j) continue;
          for (int N=0; N<=(nmax-L-l)/2; N++){
            for (int n=0; n<=(nmax-2*N-L-l)/2; n++){
              NLnlSJrel.push_back({N,L,n,l,S,Jrel});
            }
          }
        }
      }
    }
  }
  Nstates = NLnlSJrel.size();
}

OpRelCM::~OpRelCM()
{}

OpRelCM::OpRelCM(ModelSpace & ms, int rankJ, int rankP, int rankTz)
  : modelspace(&ms), rank_J(rankJ), rank_P(rankP), rank_Tz(rankTz)
{
  int E2max = modelspace->GetE2max();
  for (auto & tbc: modelspace->TwoBodyChannels){
    RelCMChannel tmp(tbc.J, tbc.Parity, tbc.Tz, e2max);
    RelCMChannels.push_back(tmp);
    IndexList[{tbc.J, tbc.Parity, tbc.Tz}] = modelspace->GetTwoBodyChannelIndex(tbc.J, tbc.Parity, tbc.Tz);
  }

  for (int ch=0; ch < modelspace->GetNumberTwoBodyChannels(); ch++){
    TwoBodyChannel & tbc = modelspace->GetTwoBodyChannel(ch);
    RelCMChannel & relcm = GetRelCMChannel(ch);
    Tcoefs[ch] = SetTcoefsChannel(tbc, relcm);
  }

  for (int ch_bra=0; ch_bra<modelsapce->GetNumberTwoBodyChannels(); ch_bra++){
    RelCMChannel & relcm_bra = GetRelCMChannel(ch_bra);
    int Jbra = relcm_bra.J;
    int Pbra = relcm_bra.Parity;
    int Zbra = relcm_bra.Tz;

    for (int ch_ket=0; ch_ket<modelsapce->GetNumberTwoBodyChannels(); ch_ket++){
      RelCMChannel & relcm_ket = GetRelCMChannel(ch_ket);
      int Jket = relcm_ket.J;
      int Pket = relcm_ket.Parity;
      int Zket = relcm_ket.Tz;
      if( std::abs(Jbra-Jket) > rank_J) continue;
      if(          Jbra+Jket  < rank_J) continue;
      if( (Pbra + Pket + rank_P)%2 == 1) continue;
      if( std::abs(Zbra-Zket) != rank_Tz) continue;
      MatEl[{ch_bra,ch_ket}] = arma::mat(relcm_bra.GetNumberStates(), relcm_ket.GetNumberStates(), arma::fill::zero);
    }
  }
}

void OpRelCM::TransLabToRelCM(Operator & op)
{
  for (auto & it: MatEl){
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    MatEl[{ch_bra,ch_ket}] = Tcoefs[ch_bra] * op.TwoBody.MatEl[{ch_bra,ch_ket}] * Tcoefs[ch_ket].t();
  }
}

Operator OpRelCM::TransRelCMToLab(ModelSpace & ms);
{
  Operator op = Operator(ms, rank_J, rank_P, rank_Tz);
  for (auto & it: MatEl){
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    op.TwoBody.MatEl[{ch_bra,ch_ket}] = Tcoefs[ch_bra].t() * MatEl[{ch_bra,ch_ket}] * Tcoefs[ch_ket];
  }
  return op;
}

arma::mat OpRelCM::SetTcoefsChannel(TwoBodyChannel & tbc, RelCmChannel & relcm)
{
  arma::mat Mat(relcm.GetNumberStates(), tbc.GetNumberKets(), arma::fill::zeros);
  int J = tbc.J;
  int Tz = tbc.Tz;
  for (int lab=0; lab<tbc.GetNumberKets(); lab++){
    Ket & ket = tbc.GetKet(lab);
    Orbit & oa = modelspace->GetOrbit(ket.p);
    Orbit & ob = modelspace->GetOrbit(ket.q);
    bool ex_ab = false;
    int na = oa.n; int la = oa.l; int ja = oa.j2; int tza = oa.tz2;
    int nb = ob.n; int lb = ob.l; int jb = ob.j2; int tzb = ob.tz2;
    if(tza == 1 and tzb ==-1){
      ex_ab = true;
      int na = ob.n; int la = ob.l; int ja = ob.j2;
      int nb = oa.n; int lb = oa.l; int jb = oa.j2;
    }

    for (int rel=0; rel<relcm.GetNumberStates(); rel++){
      std::array<int,6> tmp = relcm.NLnlSJrel[rel];
      int N = tmp[0]; int L = tmp[1];
      int n = tmp[2]; int l = tmp[3]; int S = tmp[4]; int Jrel = tmp[5];
      double tc = Tcoef(N, L, n, l, S, Jrel, na, la, ja, nb, lb, jb, J, Tz) / sqrt(1.0+ket.delta_pq());
      if(ex_ab) tc *= ket.Phase(J);
      Mat(rel,lab) = tc;
    }
  }
  return Mat;
}

double OpRelCM::Tcoef(int N, int L, int n, int l, int S, int Jrel, int na, int la, int ja, int nb, int lb, int jb, int J, int Tz)
{
  double factor = 0.0;
  if(std::abs(Tz) == 1) factor = SQRT2 * ( 1 - (la+lb+S-L)%2 );
  if(Tz == 0) factor = 1.0;
  if(std::abs(factor) < 1.e-10) return 0.0;

  double tc = 0.0;
  lam_min = std::max ( std::max( std::abs(la-lb), std::abs(L-l) ), std::abs(J-S) );
  lam_max = std::min ( std::min(          la+lb,           L+l  ),          J+S  );
  for (int lam=lam_min; lam<=lam_max; lam++){
    tc += (double)(1 - 2*(L+l+s+J)%2) * (double)(2*lam+1) * sqrt( (2*S+1) * (ja+1) * (jb+1) * (2*Jrel+1) ) *
      AngMom::NineJ(la, lb, lam, 0.5, 0.5, S, ja*0.5, jb*0.5, J) *
      AngMom::SixJ(L, l, lam, S, J, Jrel) *
      AngMom::Moshinsky(N, L, n, l, na, la, nb, lb, lam);
  }
  return tc * factor;
}

double OpRelCM::GetME(int Nbra, int Lbra, int nbra, int lbra, int Sbra, int Jrbra, int Jbra, int Tzbra,
    int Nket, int Lket, int nket, int lket, int Sket, int Jrket, int Jket, int Tzket);
{
  int Pbra = (Lbra + lbra)%2;
  int Pket = (Lket + lket)%2;
  if( std::abs(Jbra-Jket) > rank_J) return 0.0;
  if(          Jbra+Jket  < rank_J) return 0.0;
  if( (Pbra + Pket + rank_P)%2 == 1) return 0.0;
  if( std::abs(Zbra-Zket) != rank_Tz) return 0.0;
  int ch_bra = GetRelCMChannelIndex(Jbra, Pbra, Tzbra);
  int ch_ket = GetRelCMChannelIndex(Jket, Pket, Tzket);
  RelCMChannel & relcm_bra = GetRelCMChannel(ch_bra);
  RelCMChannel & relcm_ket = GetRelCMChannel(ch_ket);
  int bra = relcm_bra.GetRelCMStateIndex(Nbra, Lbra, nbra, lbra, Sbra, Jrbra);
  int ket = relcm_ket.GetRelCMStateIndex(Nket, Lket, nket, lket, Sket, Jrket);
  return GetME(ch_bra, ch_ket, bra, ket);
}

double OpRelCM::GetME(int Nbra, int Lbra, int nbra, int lbra, int Sbra, int Jrbra,
    int Nket, int Lket, int nket, int lket, int Sket, int Jrket, int J, int Tz);
{
  int Pbra = (Lbra + lbra)%2;
  int Pket = (Lket + lket)%2;
  if( (Pbra + Pket)%2 == 1) return 0.0;
  int ch = GetRelCMChannelIndex(J, Pket, Tz);
  RelCMChannel & relcm = GetRelCMChannel(ch);
  int bra = relcm.GetRelCMStateIndex(Nbra, Lbra, nbra, lbra, Sbra, Jrbra);
  int ket = relcm.GetRelCMStateIndex(Nket, Lket, nket, lket, Sket, Jrket);
  return GetME(ch, ch, bra, ket);
}
