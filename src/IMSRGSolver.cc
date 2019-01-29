
#include "IMSRGSolver.hh"
#include "Commutator.hh"
#include "Operator.hh"
#include <iomanip>

#ifndef NO_ODE
#include <boost/numeric/odeint.hpp>
#endif

IMSRGSolver::~IMSRGSolver()
{
  CleanupScratch();
}

IMSRGSolver::IMSRGSolver()
    : rw(NULL),s(0),ds(0.1),ds_max(0.5),
     norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
     flowfile(""), n_omega_written(0),max_omega_written(50),magnus_adaptive(true),hunter_gatherer(false)
     ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{}

// Constructor
IMSRGSolver::IMSRGSolver( Operator &H_in)
   : modelspace(H_in.GetModelSpace()),rw(NULL), H_0(&H_in), FlowingOps(1,H_in), Eta(H_in),
    istep(0), s(0),ds(0.1),ds_max(0.5),
    smax(2.0), norm_domega(0.1), omega_norm_max(2.0),eta_criterion(1e-6),method("magnus_euler"),
    flowfile(""), n_omega_written(0),max_omega_written(500),magnus_adaptive(true),hunter_gatherer(false)
    ,ode_monitor(*this),ode_mode("H"),ode_e_abs(1e-6),ode_e_rel(1e-6)
{
   Eta.Erase();
   Eta.SetAntiHermitian();
//   Omega.push_back( Eta);
   Omega.emplace_back( Eta);
}


void IMSRGSolver::NewOmega()
{
  H_saved = FlowingOps[0];
  std::cout << "pushing back another Omega. Omega.size = " << Omega.size()
            << " , operator size = " << Omega.front().Size()/1024./1024. << " MB"
            << ",  memory usage = " << profiler.CheckMem()["RSS"]/1024./1024. << " GB"
            << std::endl;
  if ((rw != NULL) and (rw->GetScratchDir() !=""))
  {

//    char tmp[512];
//    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), n_omega_written);
//    std::string fname(tmp);
//    std::ofstream ofs(fname, std::ios::binary);
    std::ostringstream filename;
    filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << n_omega_written;
    std::ofstream ofs(filename.str(), std::ios::binary);
    Omega.back().WriteBinary(ofs);
    if (Omega.back().GetModelSpace() != Eta.GetModelSpace()) Omega.back() = Eta;
    n_omega_written++;
//    std::cout << "Omega written to file " << fname << "  written " << n_omega_written << " so far." << std::endl;
    std::cout << "Omega written to file " << filename.str() << "  written " << n_omega_written << " so far." << std::endl;
    if (n_omega_written > max_omega_written)
    {
      std::cout << "n_omega_written > max_omega_written.  (" << n_omega_written << " > " << max_omega_written
                << " ) deleting OMEGA files and calling terminate." << std::endl;
      CleanupScratch();
      std::terminate();
    }
  }
  else
  {
    Omega.emplace_back(Eta);
  }
  Omega.back().Erase();

}


// Use a hunter-gatherer approach to finding Omega.
// We only store two Omega operators, the "hunter" and the "gatherer".
// The hunter is updated by the generator eta with the BCH formula
// exp[Omega_h(s+ds)] = exp[eta(s)]*exp[Omega_h(s)]
// When the norm of the hunter gets to a threshold set by omega_norm_max,
// we gather it, also using the BCH formula
// exp[Omega_g] = exp[Omega_h] * exp[Omega_g].
// We then clear out the hunter, and update H(s) according to the gathered Omega.
// And we continue hunting, using H(s) = exp[Omega_g] H(0) exp[-Omega_g] as our starting point.
// This is essentially a compromise between the "split" and "no split" approaches, combining
// the advantages of both. As we hunt, we have a relatively small Omega, so we don't need to
// evaluate as many nested commutators, but in the end we have just one Omega so that if
// we want to transform some operator we don't need to do a bunch of transformations.
void IMSRGSolver::GatherOmega()
{
  std::cout << "gathering Omega. "  << std::endl;
  if (Omega.size()<2)
  {
    auto& last = Omega.back();
    Omega.emplace_back( last );
  }
  // the last omega in the list is the hunter. the one just preceeding it is the gatherer.
  auto& hunter = Omega.back();
  auto& gatherer = Omega[ Omega.size()-2];
  if ( hunter.Norm() > 1e-6 )
  {
    gatherer = Commutator::BCH_Product( hunter, gatherer);
  }
  hunter.Erase();
  H_saved = *H_0;
  for (size_t i=0; i<Omega.size()-1; i++)
  {
    H_saved = Commutator::BCH_Transform( H_saved, Omega[i] );
  }
}



void IMSRGSolver::SetHin( Operator & H_in)
{
   modelspace = H_in.GetModelSpace();
   H_0 = &H_in;
//   H_s = H_in;
   FlowingOps[0] = H_in;
   Eta = H_in;
   Eta.Erase();
   Eta.SetAntiHermitian();
   if (Omega.back().Norm() > 1e-6)
   {
    NewOmega();
   }
   else
   {
     Omega.back() = Eta;
   }
}



void IMSRGSolver::SetOmega(size_t i, Operator& om)
{
 if ((i+1)> Omega.size())
 {
  Omega.resize(i+1);
 }
 Omega[i] = om;
}

void IMSRGSolver::Reset()
{
   s=0;
   Eta.Erase();
   Omega.resize(0);
   NewOmega();
}

void IMSRGSolver::SetGenerator(std::string gen)
{
  generator.SetType(gen);
  if (Omega.back().Norm() > 1e-6)
  {
    Eta.Erase();
    NewOmega();
  }
}

void IMSRGSolver::SetFlowFile(std::string str)
{
   flowfile = str;
   std::ofstream flowf;
   if (flowfile != "")
   {
      flowf.open(flowfile,std::ofstream::out);
      flowf.close();
   }
}


void IMSRGSolver::Solve()
{
  if (s<1e-4)
   WriteFlowStatusHeader(std::cout);


  if (method == "magnus_euler" or method =="magnus")
    Solve_magnus_euler();
  else if (method == "magnus_modified_euler")
    Solve_magnus_modified_euler();
  else if (method == "flow_adaptive" or method == "flow")
    Solve_ode_adaptive();
  else if (method == "magnus_adaptive")
    Solve_ode_magnus();
  else if (method == "flow_euler")
    Solve_ode();
  else if (method == "restore_4th_order")
  {
    FlowingOps.emplace_back( Operator( *(FlowingOps[0].GetModelSpace()), 0,0,0,1));
    Solve_ode_adaptive();
  }
  else
    std::cout << "IMSRGSolver: I don't know method " << method << std::endl;
}

void IMSRGSolver::UpdateEta()
{
   generator.Update(&FlowingOps[0],&Eta);
}


void IMSRGSolver::Solve_magnus_euler()
{
   istep = 0;
   generator.Update(&FlowingOps[0],&Eta);

   if (generator.GetType() == "shell-model-atan")
   {
     generator.SetDenominatorCutoff(1.0); // do we need this?
   }

    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(std::cout);

   for (istep=1;s<smax;++istep)
   {

      double norm_eta = Eta.Norm();
      if (norm_eta < eta_criterion )
      {
        break;
      }
      double norm_omega = Omega.back().Norm();
      if (norm_omega > omega_norm_max)
      {
        if (hunter_gatherer)
        {
          GatherOmega();
        }
        else
        {
          NewOmega();
        }
        norm_omega = 0;
      }
      // ds should never be more than 1, as this is over-rotating
      if (magnus_adaptive)
         ds = std::min( std::min( std::min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), omega_norm_max/norm_eta), ds_max);
      ds = std::min(ds,smax-s);
//      if (s+ds > smax) ds = smax-s;
      s += ds;
      Eta *= ds; // Here's the Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
//      Omega.back() = Eta.BCH_Product( Omega.back() );
      Omega.back() = Commutator::BCH_Product( Eta, Omega.back() );

      // transformed Hamiltonian H_s = exp(Omega) H_0 exp(-Omega)
      if ((Omega.size()+n_omega_written)<2)
      {
//        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
        FlowingOps[0] = Commutator::BCH_Transform( *H_0, Omega.back() );
      }
      else
      {
//        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
        FlowingOps[0] = Commutator::BCH_Transform( H_saved, Omega.back() );
      }

      if (norm_eta<1.0 and generator.GetType() == "shell-model-atan")
      {
        generator.SetDenominatorCutoff(1e-6);
      }

      generator.Update(&FlowingOps[0],&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(std::cout);
//      profiler.PrintMemory();

   }

}


void IMSRGSolver::Solve_magnus_modified_euler()
{
   istep = 0;
   generator.Update(&FlowingOps[0],&Eta);

   Operator H_temp;
    // Write details of the flow
   WriteFlowStatus(flowfile);
   WriteFlowStatus(std::cout);

   for (istep=1;s<smax;++istep)
   {
      double norm_eta = Eta.Norm();
      double norm_omega = Omega.back().Norm();
      if (norm_omega > omega_norm_max)
      {
        NewOmega();
        norm_omega = 0;
      }
      // ds should never be more than 1, as this is over-rotating
      ds = std::min( std::min(norm_domega/norm_eta, norm_domega / norm_eta / (norm_omega+1.0e-9)), ds_max);
      if (s+ds > smax) ds = smax-s;
      s += ds;

      H_temp = FlowingOps[0] + ds * Commutator::Commutator(Eta,FlowingOps[0]);
      generator.AddToEta(&H_temp,&Eta);

      Eta *= ds*0.5; // Here's the modified Euler step.

      // accumulated generator (aka Magnus operator) exp(Omega) = exp(dOmega) * exp(Omega_last)
//      Omega.back() = Eta.BCH_Product( Omega.back() );
      Omega.back() = Commutator::BCH_Product( Eta, Omega.back() );

      if ((Omega.size()+n_omega_written)<2)
      {
//        FlowingOps[0] = H_0->BCH_Transform( Omega.back() );
        FlowingOps[0] = Commutator::BCH_Transform( *H_0, Omega.back() );
      }
      else
      {
//        FlowingOps[0] = H_saved.BCH_Transform( Omega.back() );
        FlowingOps[0] = Commutator::BCH_Transform( H_saved, Omega.back() );
      }

      generator.Update(&FlowingOps[0],&Eta);

      // Write details of the flow
      WriteFlowStatus(flowfile);
      WriteFlowStatus(std::cout);

   }

}


#ifndef NO_ODE

// Implement element-wise division and abs and reduce for Operators.
// This is required for adaptive steppers
//vector<Operator> operator/ (const vector<Operator>& num, const vector<Operator>& denom)
std::deque<Operator> operator/ (const std::deque<Operator>& num, const std::deque<Operator>& denom)
{
//   vector<Operator> quotient = num;
   std::deque<Operator> quotient = num;
   for ( size_t i=0;i<num.size();++i )
   {
     quotient[i].ZeroBody /= denom[i].ZeroBody;
     quotient[i].OneBody /= denom[i].OneBody;
     for ( auto& itmat: quotient[i].TwoBody.MatEl )    itmat.second /= denom[i].TwoBody.GetMatrix(itmat.first[0],itmat.first[1]);
   }
   return quotient;
}

//vector<Operator> operator* (const double a, const vector<Operator>& X)
std::deque<Operator> operator* (const double a, const std::deque<Operator>& X)
{
//  vector<Operator> Y = X;
  std::deque<Operator> Y = X;
  for ( auto& y : Y )  y *= a;
  return Y;
}

//vector<Operator> operator+ ( const vector<Operator>& X, const vector<Operator>& Y)
std::deque<Operator> operator+ ( const std::deque<Operator>& X, const std::deque<Operator>& Y)
{
//  vector<Operator> Z = X;
  std::deque<Operator> Z = X;
  for ( size_t i=0;i<Z.size();++i )  Z[i] += Y[i];
  return Z;
}

// Also need the dubious operation of adding a double to an operator.
//vector<Operator> operator+ (const double a, const vector<Operator>& X)
std::deque<Operator> operator+ (const double a, const std::deque<Operator>& X)
{
//   vector<Operator> Y = X;
   std::deque<Operator> Y = X;
   for ( auto& y : Y )
   {
     y.ZeroBody += a;
     y.OneBody += a;
//     for( auto& v : y.OneBody ) v += a;
     for ( auto& itmat: y.TwoBody.MatEl )
      itmat.second += a;
//       for ( auto& v : itmat.second ) v += a;
   }
   return Y;
}

// Return the element-wise absolute value of an operator
// this is needed for ODE adaptive solver
//vector<Operator> abs(const vector<Operator>& OpIn)
std::deque<Operator> abs(const std::deque<Operator>& OpIn)
{
//   vector<Operator> OpOut = OpIn;
   std::deque<Operator> OpOut = OpIn;
   for (auto& opout : OpOut )
   {
     opout.ZeroBody = std::abs(opout.ZeroBody);
     opout.OneBody = arma::abs(opout.OneBody);
     for ( auto& itmat : opout.TwoBody.MatEl )    itmat.second = arma::abs(itmat.second);
   }
   return OpOut;
}

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION < 1.56
#ifdef OLD_BOOST
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_reduce< std::deque<Operator> >
{
   template<class Op>
   double operator()(const std::deque<Operator>& X, Op op, double init)
   {
      for (auto& x : X)
      {
        init = op(init,x.ZeroBody);
        for ( auto& v : x.OneBody )    init = op(init,v);
        for ( auto& itmat : x.TwoBody.MatEl )
        {
            for (auto& v : itmat.second )    init = op(init,v);
        }
      }
      return init;
   }
};
}}}
#endif

// Apply operation to each element of X and return the result
// this is needed for ODE adaptive solver
// USE THIS FOR BOOST VERSION >= 1.56
//struct vector_space_norm_inf< vector<Operator> >
#ifndef OLD_BOOST
namespace boost {namespace numeric {namespace odeint{
template<>
struct vector_space_norm_inf< std::deque<Operator> >
{
   typedef double result_type;
//   double operator()(const vector<Operator>& X)
   double operator()(const std::deque<Operator>& X)
   {
     double norm = 0;
     for ( auto& x : X )
       norm += x.Norm();
     return norm;
   }
};
}}}
#endif

void IMSRGSolver::Solve_ode()
{

   ode_mode = "H";
   WriteFlowStatusHeader(std::cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
//   runge_kutta4< vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   runge_kutta4< std::deque<Operator>, double, std::deque<Operator>, double, vector_space_algebra> stepper;
   auto system = *this;
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   integrate_const(stepper, system, FlowingOps, s, smax, ds, monitor);
   monitor.report();
}

void IMSRGSolver::Solve_ode_adaptive()
{
   ode_mode = "H";
   if (method == "restore_4th_order") ode_mode = "Restored";
   WriteFlowStatusHeader(std::cout);
   WriteFlowStatus(flowfile);
   std::cout << "done writing header and status" << std::endl;
   using namespace boost::numeric::odeint;
   auto system = *this;
//   typedef runge_kutta_dopri5< vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
   typedef runge_kutta_dopri5< std::deque<Operator> , double , std::deque<Operator> ,double , vector_space_algebra > stepper;
//   typedef adams_bashforth_moulton< 4, vector<Operator> , double , vector<Operator> ,double , vector_space_algebra > stepper;
   auto monitor = ode_monitor;
//   size_t steps = integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   integrate_adaptive(make_controlled<stepper>(ode_e_abs,ode_e_rel), system, FlowingOps, s, smax, ds, monitor);
   monitor.report();

}

// Evaluate dx/dt for boost ode
//void IMSRGSolver::operator()( const Operator& x, Operator& dxdt, const double t)
//void IMSRGSolver::operator()( const vector<Operator>& x, vector<Operator>& dxdt, const double t)
void IMSRGSolver::operator()( const std::deque<Operator>& x, std::deque<Operator>& dxdt, const double t)
{
   s = t;
   if (ode_mode == "H")
   {
     FlowingOps[0] = x[0];
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     auto& H_s = FlowingOps[0];
     generator.Update(&H_s,&Eta);
     if (Eta.Norm() < eta_criterion)
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = 0*x[i];
       }
     }
     else
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = Commutator::Commutator(Eta,x[i]);
       }
     }
   }
   else if (ode_mode == "Omega")
   {

     double norm_omega = Omega.back().Norm();
     if (norm_omega > omega_norm_max) // This doesn't seem to works so well just yet...
     {
       NewOmega();
       norm_omega = 0;
     }
     Omega.back() = x.back();
     auto& Omega_s = x.back();
     Operator& H_s = FlowingOps[0];
     if ((Omega.size()+n_omega_written) > 1)
       H_s = Commutator::BCH_Transform(H_saved, Omega_s);
//       H_s = H_saved.BCH_Transform(Omega_s);
     else
       H_s = Commutator::BCH_Transform(*H_0, Omega_s);
//       H_s = H_0->BCH_Transform(Omega_s);
     generator.Update(&H_s,&Eta);
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     dxdt.back() = Eta - 0.5*Commutator::Commutator(Omega_s,Eta);
   }
   else if (ode_mode == "Restored" )
   {
     FlowingOps[0] = x[0];
     FlowingOps[1] = x[1];
     if (dxdt.size() < x.size()) dxdt.resize(x.size());
     dxdt[1] = Operator(x[1]);
     auto& H_s = FlowingOps[0];
     generator.Update(&H_s,&Eta);
     if (Eta.Norm() < eta_criterion)
     {
       for (size_t i=0;i<x.size();++i)
       {
         dxdt[i] = 0*x[i];
       }
     }
     else
     {
       dxdt[0] = Commutator::Commutator(Eta,x[0]+x[1]);
       dxdt[1].Erase();
//       dxdt[1].comm221ss(Eta,x[0]);
       Commutator::comm221ss(Eta,x[0],dxdt[1]);
       // keep only pp and hh parts of d chi/ ds
       for (auto& a : modelspace->holes)
       {
         for (auto& i : modelspace->particles)
         {
           dxdt[1].OneBody(a,i) = 0;
           dxdt[1].OneBody(i,a) = 0;
         }
       }
       for (size_t i=2;i<x.size();++i)
       {
         dxdt[i] = Commutator::Commutator(Eta,x[i]);
       }
     }

   }
   WriteFlowStatus(flowfile);
   WriteFlowStatus(std::cout);
}



void IMSRGSolver::Solve_ode_magnus()
{
   ode_mode = "Omega";
   WriteFlowStatus(std::cout);
   WriteFlowStatus(flowfile);
   using namespace boost::numeric::odeint;
   namespace pl = std::placeholders;
//   runge_kutta4<vector<Operator>, double, vector<Operator>, double, vector_space_algebra> stepper;
   runge_kutta4<std::deque<Operator>, double, std::deque<Operator>, double, vector_space_algebra> stepper;
   auto system = *this;
   auto monitor = ode_monitor;
//   size_t steps = integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   integrate_const(stepper, system, Omega, s, smax, ds, monitor);
   monitor.report();
}



#endif


/// Returns \f$ e^{Omega} \mathcal{O} e^{-Omega} \f$
Operator IMSRGSolver::Transform(Operator& OpIn)
{
  return Transform_Partial(OpIn, 0);
}

Operator IMSRGSolver::Transform(Operator&& OpIn)
{
  return Transform_Partial(OpIn, 0);
}




/// Returns \f$ e^{-Omega} \mathcal{O} e^{Omega} \f$
Operator IMSRGSolver::InverseTransform(Operator& OpIn)
{
//  if (OpIn.GetJRank()+OpIn.GetTRank()+OpIn.GetParity()>0)
//  {
//    OpIn.ResetTensorTransformFirstPass();
//  }
  Operator OpOut = OpIn;
  for (auto omega=Omega.rbegin(); omega !=Omega.rend(); ++omega )
  {
    Operator negomega = -(*omega);
//    OpOut = OpOut.BCH_Transform( negomega );
    OpOut = Commutator::BCH_Transform( OpOut, negomega );
  }
  return OpOut;
}

/// Returns \f$ e^{\Omega} \mathcal{O} e^{-\Omega} \f$
/// for the \f$\Omega_i\f$s with index greater than or equal to n.
Operator IMSRGSolver::Transform_Partial(Operator& OpIn, int n)
{
//  cout << "Begin Transform_Partial" << endl;
  Operator OpOut = OpIn;
  if ((rw != NULL) and rw->GetScratchDir() != "")
  {
    Operator omega(OpIn);
//    char tmp[512];
    for (int i=n;i<n_omega_written;i++)
    {
 //    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
//     std::string fname(tmp);
     std::ostringstream filename;
     filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
     std::ifstream ifs(filename.str(),std::ios::binary);
     omega.ReadBinary(ifs);
//     if (OpIn.GetJRank()>0) cout << "step " << i << endl;
//     OpOut = OpOut.BCH_Transform( omega );
     OpOut = Commutator::BCH_Transform( OpOut, omega );
//     if (OpIn.GetJRank()>0)cout << "done" << endl;
    }
  }

  for (size_t i=std::max(n-n_omega_written,0); i<Omega.size();++i)
  {
//     if (OpIn.GetJRank()>0) cout << "step " << i << endl;
//    OpOut = OpOut.BCH_Transform( Omega[i] );
    OpOut = Commutator::BCH_Transform( OpOut, Omega[i] );
//     if (OpIn.GetJRank()>0)cout << "done" << endl;
  }

  return OpOut;
}


Operator IMSRGSolver::Transform_Partial(Operator&& OpIn, int n)
{
//  cout << "Calling r-value version of Transform_Partial, n = " << n << endl;
  Operator OpOut = OpIn;
  if ((rw != NULL) and rw->GetScratchDir() != "")
  {
    Operator omega(OpIn);
//    char tmp[512];
    for (int i=n;i<n_omega_written;i++)
    {
//     sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
//     std::string fname(tmp);
     std::ostringstream filename;
     filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
     std::ifstream ifs(filename.str(),std::ios::binary);
     omega.ReadBinary(ifs);
//     OpOut = OpOut.BCH_Transform( omega );
     OpOut = Commutator::BCH_Transform( OpOut, omega );
    }
  }

  for (size_t i=std::max(n-n_omega_written,0); i<Omega.size();++i)
  {
//    OpOut = OpOut.BCH_Transform( Omega[i] );
    OpOut = Commutator::BCH_Transform( OpOut, Omega[i] );
  }
  return OpOut;
}

// count number of equations to be solved
int IMSRGSolver::GetSystemDimension()
{
   int dim = 1; // zero-body part

   int N = H_0->OneBody.n_cols;
   dim += N*(N+1)/2;
   dim += H_0->TwoBody.Dimension();
   return dim;
}



void IMSRGSolver::CleanupScratch()
{
  if (n_omega_written<=0) return;
  std::cout << "Cleaning up files written to scratch space" << std::endl;
//  char tmp[512];
  for (int i=0;i<n_omega_written;i++)
  {
    std::ostringstream filename;
    filename << rw->GetScratchDir().c_str() << "/OMEGA_" << std::setw(6) << std::setfill('0') << getpid() << std::setw(3) << std::setfill('0') << i;
//    sprintf(tmp,"%s/OMEGA_%06d_%03d",rw->GetScratchDir().c_str(), getpid(), i);
//    std::string fname(tmp);
//    if ( remove(tmp) !=0 )
    if ( remove(filename.str().c_str()) !=0 )
    {
      std::cout << "Error when attempting to delete " << filename.str() << std::endl;
    }
  }
}



void IMSRGSolver::WriteFlowStatus(std::string fname)
{
   if (fname !="")
   {
     std::ofstream ff(fname,std::ios::app);
     WriteFlowStatus(ff);
   }
}
void IMSRGSolver::WriteFlowStatus(std::ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 18;
      int fprecision = 6;
      auto& H_s = FlowingOps[0];
      f.setf(std::ios::fixed);
      f << std::fixed << std::setw(5) << istep
        << std::setw(10) << std::setprecision(3) << s
        << std::setw(fwidth) << std::setprecision(fprecision) << H_s.ZeroBody
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 0, 0)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 1, 1)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 2, 2)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 3, 3)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 4, 4)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 5, 5)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 6, 6)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 7, 7)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 8, 8)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody( 9, 9)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(10,10)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(11,11)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(12,12)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(13,13)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(14,14)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(15,15)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(16,16)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(17,17)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(18,18)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(19,19)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(20,20)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(21,21)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(22,22)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(23,23)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(24,24)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(25,25)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(26,26)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(27,27)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(28,28)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBody(29,29)
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.Norm()
      //<< std::setw(fwidth) << std::setprecision(fprecision) << H_s.Trace( modelspace->GetAref(), modelspace->GetZref() )
//        << std::setw(fwidth) << std::setprecision(fprecision) << H_s.OneBodyNorm()
//        << std::setw(fwidth) << std::setprecision(fprecision) << H_s.TwoBodyNorm()
//        << std::setw(fwidth) << std::setprecision(fprecision) << Omega.Norm()
        << std::setw(fwidth) << std::setprecision(fprecision) << Omega.back().OneBodyNorm()
        << std::setw(fwidth) << std::setprecision(fprecision) << Omega.back().TwoBodyNorm()
        << std::setw(fwidth) << std::setprecision(fprecision) << Eta.OneBodyNorm()
        << std::setw(fwidth) << std::setprecision(fprecision) << Eta.TwoBodyNorm()
        << std::setw(7)      << std::setprecision(0)          << profiler.counter["N_ScalarCommutators"] + profiler.counter["N_TensorCommutators"]
        << std::setw(fwidth) << std::setprecision(fprecision) << H_s.GetMP2_Energy()
        << std::setw(7)      << std::setprecision(0)          << profiler.counter["N_Operators"]
        << std::setprecision(fprecision)
        << std::setw(12) << std::setprecision(3) << profiler.GetTimes()["real"]
        << std::setw(12) << std::setprecision(3) << profiler.CheckMem()["RSS"]/1024. << " / " << std::skipws << profiler.MaxMemUsage()/1024. << std::fixed
        << std::endl;
   }

}

void IMSRGSolver::WriteFlowStatusHeader(std::string fname)
{
   std::ofstream ff;
   if (fname !="") ff.open(fname,std::ios::app);
   WriteFlowStatusHeader(ff);
}
void IMSRGSolver::WriteFlowStatusHeader(std::ostream& f)
{
   if ( f.good() )
   {
      int fwidth = 16;
      int fprecision = 9;
      f.setf(std::ios::fixed);
      f << std::fixed << std::setw(5) << "i"
        << std::setw(10) << std::setprecision(3) << "s"
        << std::setw(fwidth) << std::setprecision(fprecision) << "E0"
//        << std::setw(fwidth) << std::setprecision(fprecision) << "||H_1||"
//        << std::setw(fwidth) << std::setprecision(fprecision) << "||H_2||"
        //<< std::setw(fwidth) << std::setprecision(fprecision) << "||H||"
        //<< std::setw(fwidth) << std::setprecision(fprecision) << "Tr(H)/Tr(1)"
        << std::setw(fwidth) << std::setprecision(fprecision) << "||Omega_1||"
        << std::setw(fwidth) << std::setprecision(fprecision) << "||Omega_2||"
        << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_1||"
        << std::setw(fwidth) << std::setprecision(fprecision) << "||Eta_2||"
        << std::setw(7)      << std::setprecision(fprecision) << "Ncomm"
        << std::setw(16)     << std::setprecision(fprecision) << "E(MP2)"
        << std::setw(7)      << std::setprecision(fprecision) << "N_Ops"
        << std::setw(16) << std::setprecision(fprecision) << "Walltime (s)"
        << std::setw(19) << std::setprecision(fprecision) << "Memory (MB)"
        << std::endl;
        for (int x=0;x<175;x++) f << "-";
        f << std::endl;
   }

}

// added by T.Miyagi
void IMSRGSolver::SetMiscFile(std::string str)
{
  miscfile = str;
  std::ofstream miscf;
  if (miscfile != "")
  {
    miscf.open(miscfile,std::ofstream::out);
    miscf.close();
  }
}

void IMSRGSolver::WriteStatusMisc(std::string fname)
{
  if (fname !="")
  {
    std::ofstream ff(fname,std::ios::app);
    WriteStatusMisc(ff);
  }
}

void IMSRGSolver::WriteStatusMisc(std::ostream& f)
{
  if ( f.good() )
  {
    auto& H_s = FlowingOps[0];
    f << H_s.GetOrderedTwoBodyMonopoleMatrix(0,0) << std::endl;
    f << Eta.GetOrderedTwoBodyMonopoleMatrix(0,0) << std::endl;
  }
}
