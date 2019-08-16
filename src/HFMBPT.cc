#include <iomanip>
#include "HFMBPT.hh"
#include "HartreeFock.hh"

HFMBPT::~HFMBPT()
{}

HFMBPT::HFMBPT(Operator& hbare)
  : HartreeFock(hbare)
{}

// post Hartree-Fock method
void HFMBPT::GetNaturalOrbital()
{
  int norbits = HartreeFock::modelspace->GetNumberOrbits();
  int A = HartreeFock::modelspace->GetTargetMass();
  C_HF2NAT = arma::mat(norbits,norbits,arma::fill::eye);
  Occ      = arma::vec(norbits,arma::fill::zeros);
  GetDensityMatrix();
  DiagRho();
  double AfromTr = 0.0;
  for(int i=0; i< norbits; ++i)
  {
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    AfromTr += rho(i,i) * (oi.j2+1);
  }

  if(abs(AfromTr - A) > 1e-8)
  {
    std::cout << "Warning: Mass != Tr(rho)" << std::endl;
    exit(0);
  }
  C_HO2NAT = C * C_HF2NAT;
  // use fractional occupation
  //for(int i=0; i< norbits; ++i)
  //{
  //  Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
  //  oi.occ = Occ(i);
  //}
  //
}

void HFMBPT::DiagRho()
{
  for (auto& it : Hbare.OneBodyChannels)
  {
    arma::uvec orbvec(it.second);
    arma::uvec orbvec_d = sort(orbvec, "descend");
    arma::mat rho_ch = rho.submat(orbvec, orbvec);
    arma::mat vec;
    arma::vec eig;
    bool success = false;
    success = arma::eig_sym(eig, vec, rho_ch);
    if(not success)
    {
      std::cout << "Error in diagonalization of density matrix" << std::endl;
      std::cout << "Density Matrix:" << std::endl;
      rho_ch.print();
      exit(0);
    }
    Occ(orbvec_d) = eig;
    C_HF2NAT.submat(orbvec, orbvec_d) = vec;
  }
}

Operator HFMBPT::TransformHFToNATBasis( Operator& OpHF)
{
  Operator OpNAT(OpHF);
  OpNAT.OneBody = C_HF2NAT.t() * OpHF.OneBody * C_HF2NAT;

  for (auto& it : OpHF.TwoBody.MatEl )
  {
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    TwoBodyChannel& tbc_bra = OpNAT.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = OpNAT.modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = it.second.n_rows;
    int nkets = it.second.n_cols;
    arma::mat Dbra(nbras,nbras);
    arma::mat Dket(nkets,nkets);

    for (int i = 0; i<nkets; ++i)
    {
      Ket & ket_hf = tbc_ket.GetKet(i);
      for (int j=0; j<nkets; ++j)
      {
        Ket & ket_nat = tbc_ket.GetKet(j);
        Dket(i,j) = C_HF2NAT(ket_hf.p,ket_nat.p) * C_HF2NAT(ket_hf.q,ket_nat.q);
        if(ket_hf.p != ket_hf.q)
        {
          Dket(i,j) += C_HF2NAT(ket_hf.q, ket_nat.p) * C_HF2NAT(ket_hf.p, ket_nat.q) *
            ket_hf.Phase(tbc_ket.J);
        }
        if (ket_hf.p==ket_hf.q)    Dket(i,j) *= SQRT2;
        if (ket_nat.p==ket_nat.q)    Dket(i,j) /= SQRT2;
      }
    }
    if (ch_bra == ch_ket) {
      Dbra = Dket.t();
    }
    else
    {
      for (int i=0; i<nbras; ++i)
      {
        Ket & bra_nat = tbc_bra.GetKet(i);
        for (int j=0; j<nbras; ++j)
        {
          Ket & bra_hf = tbc_bra.GetKet(j);
          Dbra(i,j) = C_HF2NAT(bra_hf.p,bra_nat.p) * C_HF2NAT(bra_hf.q,bra_nat.q);
          if (bra_hf.p!=bra_hf.q)
          {
            Dbra(i,j) += C_HF2NAT(bra_hf.q, bra_nat.p) * C_HF2NAT(bra_hf.p, bra_nat.q)
              * bra_hf.Phase(tbc_bra.J);
          }
          if (bra_hf.p==bra_hf.q)    Dbra(i,j) *= SQRT2;
          if (bra_nat.p==bra_nat.q)    Dbra(i,j) /= SQRT2;
        }
      }
    }
    auto& IN  =  it.second;
    auto& OUT =  OpNAT.TwoBody.GetMatrix(ch_bra,ch_ket);
    OUT  =    Dbra * IN * Dket;
  }
  return OpNAT;
}

Operator HFMBPT::TransformHOToNATBasis( Operator& OpHO)
{
  Operator OpNAT(OpHO);
  OpNAT.OneBody = C_HO2NAT.t() * OpHO.OneBody * C_HO2NAT;

  for (auto& it : OpHO.TwoBody.MatEl )
  {
    int ch_bra = it.first[0];
    int ch_ket = it.first[1];
    TwoBodyChannel& tbc_bra = OpNAT.modelspace->GetTwoBodyChannel(ch_bra);
    TwoBodyChannel& tbc_ket = OpNAT.modelspace->GetTwoBodyChannel(ch_ket);
    int nbras = it.second.n_rows;
    int nkets = it.second.n_cols;
    arma::mat Dbra(nbras,nbras);
    arma::mat Dket(nkets,nkets);

    for (int i = 0; i<nkets; ++i)
    {
      Ket & ket_ho = tbc_ket.GetKet(i);
      for (int j=0; j<nkets; ++j)
      {
        Ket & ket_nat = tbc_ket.GetKet(j);
        Dket(i,j) = C_HO2NAT(ket_ho.p,ket_nat.p) * C_HO2NAT(ket_ho.q,ket_nat.q);
        if(ket_ho.p != ket_ho.q)
        {
          Dket(i,j) += C_HO2NAT(ket_ho.q, ket_nat.p) * C_HO2NAT(ket_ho.p, ket_nat.q) *
            ket_ho.Phase(tbc_ket.J);

        }
        if (ket_ho.p==ket_ho.q)    Dket(i,j) *= SQRT2;
        if (ket_nat.p==ket_nat.q)    Dket(i,j) /= SQRT2;
      }
    }
    if (ch_bra == ch_ket) {
      Dbra = Dket.t();
    }
    else
    {
      for (int i=0; i<nbras; ++i)
      {
        Ket & bra_nat = tbc_bra.GetKet(i);
        for (int j=0; j<nbras; ++j)
        {
          Ket & bra_ho = tbc_bra.GetKet(j);
          Dbra(i,j) = C_HO2NAT(bra_ho.p,bra_nat.p) * C_HO2NAT(bra_ho.q,bra_nat.q);
          if (bra_ho.p!=bra_ho.q)
          {
            Dbra(i,j) += C_HO2NAT(bra_ho.q,bra_nat.p) * C_HO2NAT(bra_ho.p,bra_nat.q) *
              bra_ho.Phase(tbc_bra.J);
          }
          if (bra_ho.p==bra_ho.q)    Dbra(i,j) *= SQRT2;
          if (bra_nat.p==bra_nat.q)    Dbra(i,j) /= SQRT2;
        }
      }
    }

    auto& IN  =  it.second;
    auto& OUT =  OpNAT.TwoBody.GetMatrix(ch_bra,ch_ket);
    OUT  =    Dbra * IN * Dket;
   }
   return OpNAT;
}

Operator HFMBPT::GetNormalOrderedHNAT()
{
  double start_time = omp_get_wtime();
  std::cout << "Getting normal-ordered H in NAT basis" << std::endl;
  arma::mat rho_swap = rho;
  arma::mat tmp = C_HO2NAT.cols(holeorbs);
  rho = (tmp.each_row() % hole_occ) * tmp.t();

  // fractional occupation
  // rho = C_HO2NAT * diagmat(Occ) * C_HO2NAT.t();
  //

  UpdateF();
  CalcEHF();
  std::cout << std::fixed <<  std::setprecision(7);
  std::cout << "e1Nat = " << e1hf << std::endl;
  std::cout << "e2Nat = " << e2hf << std::endl;
  std::cout << "e3Nat = " << e3hf << std::endl;
  std::cout << "E_Nat = "  << EHF  << std::endl;

  Operator HNO = Operator(*HartreeFock::modelspace,0,0,0,2);
  HNO.ZeroBody = EHF;
  HNO.OneBody = C_HO2NAT.t() * F * C_HO2NAT;

  int nchan = HartreeFock::modelspace->GetNumberTwoBodyChannels();
  int norb = HartreeFock::modelspace->GetNumberOrbits();
  for (int ch=0; ch<nchan; ++ch)
  {
    TwoBodyChannel& tbc = HartreeFock::modelspace->GetTwoBodyChannel(ch);
    int J = tbc.J;
    int npq = tbc.GetNumberKets();

    arma::mat D(npq,npq,arma::fill::zeros);
    arma::mat V3NO(npq,npq,arma::fill::zeros);
#pragma omp parallel for schedule(dynamic,1)
    for (int i=0; i<npq; ++i)
    {
      Ket & bra = tbc.GetKet(i);
      int e2bra = 2*bra.op->n + bra.op->l + 2*bra.oq->n + bra.oq->l;
      for (int j=0; j<npq; ++j)
      {
        Ket & ket = tbc.GetKet(j);
        int e2ket = 2*ket.op->n + ket.op->l + 2*ket.oq->n + ket.oq->l;
        D(i,j) = C_HO2NAT(bra.p,ket.p) * C_HO2NAT(bra.q,ket.q);
        if (bra.p!=bra.q)
        {
          D(i,j) += C_HO2NAT(bra.q,ket.p) * C_HO2NAT(bra.p,ket.q) * bra.Phase(J);
        }
        if (bra.p==bra.q)    D(i,j) *= SQRT2;
        if (ket.p==ket.q)    D(i,j) /= SQRT2;

        // Now generate the NO2B part of the 3N interaction
        if (Hbare.GetParticleRank()<3) continue;
        if (i>j) continue;
        for (int a=0; a<norb; ++a)
        {
          Orbit & oa = HartreeFock::modelspace->GetOrbit(a);
          if ( 2*oa.n+oa.l+e2bra > Hbare.GetE3max() ) continue;
          for (int b : Hbare.OneBodyChannels.at({oa.l,oa.j2,oa.tz2}))
          {
            Orbit & ob = HartreeFock::modelspace->GetOrbit(b);
            if ( 2*ob.n+ob.l+e2ket > Hbare.GetE3max() ) continue;
            int J3min = abs(2*J-oa.j2);
            int J3max = 2*J + oa.j2;
            for (int J3=J3min; J3<=J3max; J3+=2)
            {
              V3NO(i,j) += rho(a,b) * (J3+1) * Hbare.ThreeBody.GetME_pn(J,J,J3,bra.p,bra.q,a,ket.p,ket.q,b);
            }
          }
        }
        V3NO(i,j) /= (2*J+1);
        if (bra.p==bra.q)  V3NO(i,j) /= SQRT2;
        if (ket.p==ket.q)  V3NO(i,j) /= SQRT2;
        V3NO(j,i) = V3NO(i,j);
      }
    }
    auto& V2  =  Hbare.TwoBody.GetMatrix(ch);
    auto& OUT =  HNO.TwoBody.GetMatrix(ch);
    OUT  =    D.t() * (V2 + V3NO) * D;
  }

  rho = rho_swap;
  profiler.timer["HFMBT_GetNormalOrderedHNO"] += omp_get_wtime() - start_time;
  return HNO;
}

void HFMBPT::GetDensityMatrix()
{
  Operator Hhf = HartreeFock::GetNormalOrderedH();
  Operator& H(Hhf);
  int norbits = HartreeFock::modelspace->GetNumberOrbits();
  rho      = arma::mat(norbits,norbits,arma::fill::zeros);
  double t_start = omp_get_wtime();
  int norbits = HartreeFock::modelspace->GetNumberOrbits();
  rho = arma::mat(norbits,norbits,arma::fill::zeros);
  DensityMatrixPP(H);
  DensityMatrixHH(H);
  DensityMatrixPH(H);
  profiler.timer["HFMBPT DensityMatrix"] += omp_get_wtime() - t_start;
}

void HFMBPT::PrintOccupation()
{

  for (int i=0; i<HartreeFock::modelspace->GetNumberOrbits(); ++i)
  {
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
    std::cout << std::fixed << std::setw(4) << oi.n << std::setw(4) << oi.l <<
      std::setw(4) << oi.j2 << std::setw(4) << oi.tz2 << "   " <<
      std::setw(8) << Occ(i) << std::endl;
  }
}

void HFMBPT::DensityMatrixPP(Operator& H)
{
  for (auto& a : HartreeFock::modelspace->particles) {
    double ea = H.OneBody(a,a);
    Orbit& oa = HartreeFock::modelspace->GetOrbit(a);

    for (auto& b : HartreeFock::modelspace->particles) {
      if(b > a) continue;
      double eb = H.OneBody(b,b);
      Orbit& ob = HartreeFock::modelspace->GetOrbit(b);
      if(oa.j2 != ob.j2) continue;
      if(oa.l != ob.l) continue;
      if(oa.tz2 != ob.tz2) continue;

      double r = 0.0;
      for(auto& c : HartreeFock::modelspace->particles){
        double ec = H.OneBody(c,c);
        Orbit& oc = HartreeFock::modelspace->GetOrbit(c);

        for(auto& i : HartreeFock::modelspace->holes){
          double ei = H.OneBody(i,i);
          Orbit& oi = HartreeFock::modelspace->GetOrbit(i);

          for(auto& j : HartreeFock::modelspace->holes){
            double ej = H.OneBody(j,j);
            Orbit& oj = HartreeFock::modelspace->GetOrbit(j);

            double e_acij = ea + ec - ei - ej;
            double e_bcij = eb + ec - ei - ej;
            if(e_acij < 1.e-8) continue;
            if(e_bcij < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-oc.j2), std::max(std::abs(oi.j2-oj.j2), std::abs(ob.j2-oc.j2)))/2;
            int Jmax = std::min(         oa.j2+oc.j2,  std::min(         oi.j2+oj.j2,           ob.j2+oc.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J){
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,a,c,i,j) *
                H.TwoBody.GetTBME_J(J,i,j,b,c);
            }
            r += tbme / (e_acij * e_bcij);
          }
        }
      }
      rho(a,b) = r * 0.5 / (oa.j2+1);
      rho(b,a) = r * 0.5 / (oa.j2+1);
    }
  }
}

void HFMBPT::DensityMatrixHH(Operator& H)
{
  for (auto& i : HartreeFock::modelspace->holes) {
    double ei = H.OneBody(i,i);
    Orbit& oi = HartreeFock::modelspace->GetOrbit(i);

    for (auto& j : HartreeFock::modelspace->holes) {
      if(j > i) continue;
      double ej = H.OneBody(j,j);
      Orbit& oj = HartreeFock::modelspace->GetOrbit(j);
      if(oi.j2  != oj.j2) continue;
      if(oi.l   != oj.l) continue;
      if(oi.tz2 != oj.tz2) continue;

      double r = 0.0;
      for(auto& a : HartreeFock::modelspace->particles){
        double ea = H.OneBody(a,a);
        Orbit& oa = HartreeFock::modelspace->GetOrbit(a);

        for(auto& b : HartreeFock::modelspace->particles){
          double eb = H.OneBody(b,b);
          Orbit& ob = HartreeFock::modelspace->GetOrbit(b);

          for(auto& k : HartreeFock::modelspace->holes){
            double ek = H.OneBody(k,k);
            Orbit& ok = HartreeFock::modelspace->GetOrbit(k);

            double e_abik = ea + eb - ei - ek;
            double e_abjk = ea + eb - ek - ej;
            if(e_abik < 1.e-8) continue;
            if(e_abjk < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-ob.j2), std::max(std::abs(oi.j2-ok.j2), std::abs(oj.j2-ok.j2)))/2;
            int Jmax = std::min(         oa.j2+ob.j2,  std::min(         oi.j2+ok.j2,           oj.j2+ok.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J){
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,a,b,i,k) *
                H.TwoBody.GetTBME_J(J,j,k,a,b);
            }
            r += tbme / (e_abik * e_abjk);
          }
        }
      }
      rho(i,j) = - r * 0.5 / (oi.j2+1);
      rho(j,i) = - r * 0.5 / (oi.j2+1);
    }
    rho(i,i) += oi.occ;
  }
}

void HFMBPT::DensityMatrixPH(Operator& H)
{
  for (auto& a : HartreeFock::modelspace->particles) {
    double ea = H.OneBody(a,a);
    Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
    for (auto& i : HartreeFock::modelspace->holes) {
      double ei = H.OneBody(i,i);
      Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
      if(oa.j2  != oi.j2) continue;
      if(oa.l   != oi.l) continue;
      if(oa.tz2 != oi.tz2) continue;

      double r = 0.0;
      for(auto& b : HartreeFock::modelspace->particles){
        double eb = H.OneBody(b,b);
        Orbit& ob = HartreeFock::modelspace->GetOrbit(b);

        for(auto& c : HartreeFock::modelspace->particles){
          double ec = H.OneBody(c,c);
          Orbit& oc = HartreeFock::modelspace->GetOrbit(c);

          for(auto& j : HartreeFock::modelspace->holes){
            double ej = H.OneBody(j,j);
            Orbit& oj = HartreeFock::modelspace->GetOrbit(j);

            double e_ai = ea - ei;
            double e_bcij = eb + ec - ei - ej;
            if(e_ai < 1.e-8) continue;
            if(e_bcij < 1.e-8) continue;
            int Jmin = std::max(std::abs(oa.j2-oj.j2), std::max(std::abs(ob.j2-oc.j2), std::abs(oi.j2-oj.j2)))/2;
            int Jmax = std::min(         oa.j2+oj.j2,  std::min(         ob.j2+oc.j2,           oi.j2+oj.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J){
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,a,j,b,c) *
                H.TwoBody.GetTBME_J(J,b,c,i,j);
            }
            r += tbme / (e_ai * e_bcij);
          }
        }
      }
      rho(a,i) += r * 0.5 / (oa.j2+1);
      rho(i,a) += r * 0.5 / (oa.j2+1);
    }
  }

  for (auto& a : HartreeFock::modelspace->particles) {
    double ea = H.OneBody(a,a);
    Orbit& oa = HartreeFock::modelspace->GetOrbit(a);
    for (auto& i : HartreeFock::modelspace->holes) {
      double ei = H.OneBody(i,i);
      Orbit& oi = HartreeFock::modelspace->GetOrbit(i);
      if(oa.j2  != oi.j2) continue;
      if(oa.l   != oi.l) continue;
      if(oa.tz2 != oi.tz2) continue;

      double r = 0.0;
      for(auto& b : HartreeFock::modelspace->particles){
        double eb = H.OneBody(b,b);
        Orbit& ob = HartreeFock::modelspace->GetOrbit(b);

        for(auto& j : HartreeFock::modelspace->holes){
          double ej = H.OneBody(j,j);
          Orbit& oj = HartreeFock::modelspace->GetOrbit(j);

          for(auto& k : HartreeFock::modelspace->holes){
            double ek = H.OneBody(k,k);
            Orbit& ok = HartreeFock::modelspace->GetOrbit(k);

            double e_ai = ea - ei;
            double e_abkj = ea + eb - ek - ej;
            if(e_ai < 1.e-8) continue;
            if(e_abkj < 1.e-8) continue;
            int Jmin = std::max(std::abs(ok.j2-oj.j2), std::max(std::abs(oi.j2-ob.j2), std::abs(oa.j2-ob.j2)))/2;
            int Jmax = std::min(         ok.j2+oj.j2,  std::min(         oi.j2+ob.j2,           oa.j2+ob.j2))/2;

            double tbme = 0.0;
            for(int J = Jmin; J <= Jmax; ++J){
              tbme += (2*J+1) * H.TwoBody.GetTBME_J(J,k,j,i,b) *
                H.TwoBody.GetTBME_J(J,a,b,k,j);
            }
            r += tbme / (e_ai * e_abkj);
          }
        }
      }
      rho(a,i) -= r * 0.5 / (oa.j2+1);
      rho(i,a) -= r * 0.5 / (oa.j2+1);
    }
  }
}
