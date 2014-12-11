
#include "HartreeFock.hh"
#include "ModelSpace.hh"

using namespace std;

HartreeFock::HartreeFock(Operator& hbare) : Hbare(hbare)
{
   tolerance = 1e-10;
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();

   C    = arma::mat(norbits,norbits,arma::fill::eye);
   Vab  = arma::mat(norbits,norbits);
   H    = arma::mat(norbits,norbits);
   Vmon = arma::mat(2*nKets,2*nKets);
   prev_energies = arma::vec(norbits,arma::fill::zeros);

   t = Hbare.OneBody;
   energies = t.diag();
   UpdateDensityMatrix();
   BuildMonopoleV();
}



void HartreeFock::Solve()
{
   int iter = 0; // counter so we don't go on forever
   int maxiter = 1000;

   while( (! CheckConvergence()) && iter < maxiter)
   {
      Diagonalize();  // Diagonalize the Fock matrix
      ReorderCoefficients(); // Reorder C so that new ordering corresponds to original ordering
      UpdateDensityMatrix();  // 1 body density matrix, used in UpdateH()
      UpdateH();  // Update the Fock matrix
      iter++;
   }

   if (iter==maxiter)
   {
      cerr << "Warning: Hartree-Fock calculation didn't converge after " << maxiter << " iterations." << endl;
   }
   // Calculate the HF energy.
   CalcEHF();
}



void HartreeFock::CalcEHF()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   EHF = 0;
   int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int i=0;i<norbits;i++)
   {
      for (int j=0;j<norbits;j++)
      {
         EHF += 0.5 * rho(i,j) * (ms->GetOrbit(i).j2+1) * (H(i,j)+t(i,j));
      }
   }
}


///**********************************************************************************
// Diagonalize() -- [See Suhonen eq. 4.85]
// Diagonalize <a|H|b> and put the
// eigenvectors in C(i,alpha) = <alpha|i>
// and eigenvalues in the vector energies.
// Save the last vector of energies to check
// for convergence.
// Different channels are diagonalized independently.
// This guarantees that J,Tz, and parity remain good. 
///***********************************************************************************
void HartreeFock::Diagonalize()
{
   prev_energies = energies;
   vector<int> oit_list;
   arma::mat H_ch;
   arma::mat C_ch;
   arma::vec E_ch;
   bool success;
   int norbits = Hbare.GetModelSpace()->GetNumberOrbits();
   for (int p = 0; p<=1;p++)
   {
      for (int Tz = -1; Tz<=1; Tz+=2)
      {
          for (int J=0;J<JMAX;J++)
          {

            // Find all the SP orbits that have quantum numbers J,p,Tz
            // and store them in a list
             oit_list.resize(0);
             for (int a=0;a<norbits;a++)
             {
                //Orbit &oa = Hbare->GetModelSpace()->GetOrbit(a);
                Orbit &oa = Hbare.GetModelSpace()->GetOrbit(a);
                if (oa.j2==J and oa.tz2==Tz and (oa.l%2)==p)
                {
                   oit_list.push_back(a);
                }
             }
             // Now create submatrices corresponding to just these orbits
             int norb = oit_list.size();
             if (norb < 1) continue;
             H_ch = arma::mat(norb,norb,arma::fill::zeros);
             C_ch = arma::mat(norb,norb,arma::fill::zeros);
             E_ch = arma::vec(norb,arma::fill::zeros);
             for (int a=0;a<norb;a++)
             {
                for (int b=0;b<norb;b++)
                {
                   H_ch(a,b) = H(oit_list[a],oit_list[b]);
                }
             }
             // Diagonalize the submatrix
             success = arma::eig_sym(E_ch, C_ch, H_ch);
             // Update the full overlap matrix C and energy vector
             for (int a=0;a<norb;a++)
             {
                energies( oit_list[a] ) = E_ch(a);
                for (int b=0;b<norb;b++)
                {
                   C(oit_list[a],oit_list[b]) = C_ch(a,b);
                }
             }

          } // for J...
      } // For Tz ...
   } // For p...
}

//**************************************************
// BuildMonopoleV()
// Construct the unnormalized monople Hamiltonian
// <ab|V_mon|cd> = Sum_J (2J+1) <ab|V|cd>_J.
//  I calculate and store each term, rather than
//  repeatedly calculating things.
//   To facilitate looping, the matrix has 4 blocks.
//   If a<=b and c<=d, the matrix looks like:
//
//         [ <ab|V|cd> ... <ab|V|dc> ]
//    V  = |     :             :     |
//         [ <ba|V|cd> ... <ba|V|dc> ]
//         
//**************************************************
void HartreeFock::BuildMonopoleV()
{
   Vmon.zeros();
   int nKets = Hbare.GetModelSpace()->GetNumberKets();
   for (int ibra=0;ibra<nKets;++ibra)
   {
      Ket & bra = Hbare.GetModelSpace()->GetKet(ibra);
      int a = bra.p;
      int b = bra.q;
      Orbit & oa = Hbare.GetModelSpace()->GetOrbit(a);
      Orbit & ob = Hbare.GetModelSpace()->GetOrbit(b);
      double norm = (oa.j2+1)*(ob.j2+1);
      for (int iket=0;iket<nKets;++iket)
      {
         Ket & ket = Hbare.GetModelSpace()->GetKet(iket);
         int c = ket.p;
         int d = ket.q;
         Vmon(ibra,iket)             = Hbare.GetTBMEmonopole(a,b,c,d) * norm;
         Vmon(ibra+nKets,iket)       = Hbare.GetTBMEmonopole(b,a,c,d) * norm;
         Vmon(ibra,iket+nKets)       = Hbare.GetTBMEmonopole(a,b,d,c) * norm;
         Vmon(ibra+nKets,iket+nKets) = Hbare.GetTBMEmonopole(b,a,d,c) * norm;
      }
   }
   
}



//**************************************************************************
// 1-body density matrix 
// <i|rho|j> = Sum_beta <i|beta> <j|beta>
// where beta runs over HF orbits in the core (i.e. below the fermi surface)
//**************************************************************************
void HartreeFock::UpdateDensityMatrix()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();

   // Pcore is a projector onto core orbits.
   arma::mat Pcore = arma::mat(norbits,norbits,arma::fill::zeros);
   for (int &beta : ms->holes)
   {
       Pcore(beta,beta) = 1;
   }

   rho = C * Pcore * C.t();
}


//*********************************************************************
// UpdateH() -- [See Suhonen eq 4.85]
// <a|H|b> = <a|t|b>  +  Sum_ij <i|rho|j> <ai|V_mon|bj>
// * H is the fock matrix, to be diagonalized
// * t is the kinetic energy
// * rho is the density matrix defined in UpdateDensityMatrix()
// * V_mon is the monopole component of the 2-body interaction.
//*********************************************************************
void HartreeFock::UpdateH()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   int nKets = ms->GetNumberKets();
   int bra, ket;
   Vab.zeros();
   for (int a=0;a<norbits;a++)
   {
      for (int b=a;b<norbits;b++)
      {
         if (ms->GetOrbit(a).j2  != ms->GetOrbit(b).j2)  continue;
         if (ms->GetOrbit(a).tz2 != ms->GetOrbit(b).tz2) continue;
         if (ms->GetOrbit(a).l   != ms->GetOrbit(b).l)   continue;
         for (int i=0;i<norbits;i++)
         {
            // The monopoles are listed for |ab> fist with a<=b, then for a>=b
            // so if a>b, add nKets.
            bra = ms->GetKetIndex(min(a,i),max(a,i));
            if (a>i) bra += nKets;
            for (int j=0;j<norbits;j++)
            {
               ket = ms->GetKetIndex(min(b,j),max(b,j));
               if (b>j) ket += nKets;

               Vab(a,b) += rho(i,j)*Vmon(min(bra,ket),max(bra,ket)); // <i|rho|j> * <ai|Vmon|bj>
           }
         }
         Vab(a,b) /= (ms->GetOrbit(a).j2+1); // divide by 2ja+1
         Vab(b,a) = Vab(a,b);  // Hermitian & real => symmetric
      }
   }
   H = t + Vab;

}


//********************************************************
// Check for convergence using difference in s.p. energies
// between iterations.
//********************************************************
bool HartreeFock::CheckConvergence()
{
   arma::vec de = energies - prev_energies;
   double ediff = sqrt(arma::dot(de,de)) / energies.size();
   return (ediff < tolerance);
}



void HartreeFock::PrintOrbits()
{
  ModelSpace * ms = Hbare.GetModelSpace();
  int norbits = ms->GetNumberOrbits();
  for (int a=0;a<norbits;a++)
  {
     Orbit & oa = ms->GetOrbit(a);
     cout << a << ": E= " << oa.spe << " L=" << oa.l
          << " J=" << oa.j2 << " N=" << oa.n << " Tz=" << oa.tz2
          << " ph=" << oa.ph << " io=" << oa.io << endl;
  }

}

//**********************************************************************
// Eigenvectors/values come out of the diagonalization energy-ordered.
// We want them ordered corresponding to the input ordering, i.e. we want
// the matrix C to be maximal and positive along the diagonal.
//**********************************************************************
void HartreeFock::ReorderCoefficients()
{
   ModelSpace * ms = Hbare.GetModelSpace();
   int norbits = ms->GetNumberOrbits();
   arma::mat C_tmp = C;
   arma::vec e_tmp = energies;
   for (int i=0;i<norbits;i++)
   {
      float fmax = 0.0;
      int kmax;
      int sign = 1;
      for (int k=0;k<norbits;k++)
      {
         if ( abs(C_tmp(i,k)) > fmax )
         {
            fmax = abs(C_tmp(i,k));
            kmax = k;
            if (C_tmp(i,k) < 0) sign = -1;
         }
      }
      // make sure we have a positive coefficient for the biggest term
      C.col(kmax) = C_tmp.col(i) * sign;
      energies(kmax) = e_tmp(i);
   }
}




//**************************************************************************
// Takes in an operator expressed in the basis of the original Hamiltonian,
// and returns that operator in the Hartree-Fock basis.
//**************************************************************************
Operator HartreeFock::TransformToHFBasis( Operator& OpIn)
{
   Operator OpHF = OpIn;

   //Update the one-body part by multiplying by the matrix C(i,a) = <i|a>.
   OpHF.OneBody = C.t() * OpIn.OneBody * C;

   //Update the two-body part by multiplying by the matrix D(ij,ab) = <ij|ab>
   // for each channel J,p,Tz. Most of the effort here is in constructing D.
   int nchan = OpIn.GetModelSpace()->GetNumberTwoBodyChannels();
   for (int ch=0;ch<nchan;++ch)
   {
      TwoBodyChannel& tbc = OpIn.GetModelSpace()->GetTwoBodyChannel(ch);
      int npq = tbc.GetNumberKets();
      if (npq<1) continue;

      arma::mat D     = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ab> = <ji|ba>
      arma::mat Dexch = arma::mat(npq,npq,arma::fill::zeros);  // <ij|ba> = <ji|ab>

      // loop over all possible original basis configurations <pq| in this J,p,Tz channel.
      // and all possible HF configurations |p'q'> in this J,p,Tz channel                                    
      // bra is in the original basis, ket is in the HF basis                                              
      // i and j are the indices of the matrix D for this channel                    
      for (int i=0; i<npq; ++i)    
      {
         Ket & bra = tbc.GetKet(i);   
         for (int j=0; j<npq; ++j)    
         {
            Ket & ket = tbc.GetKet(j); // 
            double normfactor = sqrt((1.0+ket.delta_pq())/(1.0+bra.delta_pq()));
            D(i,j) = C(bra.p,ket.p) * C(bra.q,ket.q) * normfactor;
            Dexch(i,j) = C(bra.p,ket.q) * C(bra.q,ket.p) * ket.Phase(tbc.J) * (1-ket.delta_pq()) * normfactor;
         }
      }

     // Do all the matrix multiplication in one expression so Armadillo can do optimizations.
     OpHF.TwoBody[ch]   = D.t()      * OpIn.TwoBody[ch] * D
                        + Dexch.t()  * OpIn.TwoBody[ch] * D
                        + D.t()      * OpIn.TwoBody[ch] * Dexch
                        + Dexch.t()  * OpIn.TwoBody[ch] * Dexch;

   }
   return OpHF;
}

