/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                  ____                                         ///
///        _________________           _____________/   /\               _________________        ///
///       /____/_____/_____/|         /____/_____/ /___/  \             /____/_____/_____/|       ///
///      /____/_____/__G_ /||        /____/_____/|/   /\  /\           /____/_____/____ /||       ///
///     /____/_____/__+__/|||       /____/_____/|/ G /  \/  \         /____/_____/_____/|||       ///
///    |     |     |     ||||      |     |     |/___/   /\  /\       |     |     |     ||||       ///
///    |  I  |  M  |     ||/|      |  I  |  M  /   /\  /  \/  \      |  I  |  M  |     ||/|       ///
///    |_____|_____|_____|/||      |_____|____/ + /  \/   /\  /      |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |   /___/   /\  /  \/       |     |     |     ||||       ///
///    |  S  |  R  |     ||/|      |  S  |   \   \  /  \/   /        |  S  |  R  |  G  ||/|       ///
///    |_____|_____|_____|/||      |_____|____\ __\/   /\  /         |_____|_____|_____|/||       ///
///    |     |     |     ||||      |     |     \   \  /  \/          |     |     |     ||||       ///
///    |     |  +  |     ||/       |     |  +  |\ __\/   /           |     |  +  |  +  ||/        ///
///    |_____|_____|_____|/        |_____|_____|/\   \  /            |_____|_____|_____|/         ///
///                                               \___\/                                          ///
///                                                                                               ///
///           imsrg++ : Interface for performing standard IMSRG calculations.                     ///
///                     Usage is imsrg++  option1=value1 option2=value2 ...                       ///
///                     To get a list of options, type imsrg++ help                               ///
///                                                                                               ///
///                                                      - Ragnar Stroberg 2016                   ///
///                                                                                               ///
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//    imsrg++.cc, part of  imsrg++
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


#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"

//using namespace imsrg_util; // lets get rid of this

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
#ifdef BUILDVERSION
  std::cout << "######  imsrg++ build version: " << BUILDVERSION << std::endl;
#endif
  Parameters parameters(argc,argv);
  if (parameters.help_mode) return 0;

  std::string inputtbme = parameters.s("2bme");
  std::string input3bme = parameters.s("3bme");
  std::string reference = parameters.s("reference");
  std::string valence_space = parameters.s("valence_space");
  std::string custom_valence_space = parameters.s("custom_valence_space");
  std::string basis = parameters.s("basis");
  std::string method = parameters.s("method");
  std::string flowfile = parameters.s("flowfile");
  std::string smryfile = parameters.s("summaryfile");
  std::string flow1file = parameters.s("flow1file");
  std::string flow2file = parameters.s("flow2file");
  std::string intfile = parameters.s("intfile");
  std::string core_generator = parameters.s("core_generator");
  std::string valence_generator = parameters.s("valence_generator");
  std::string fmt2 = parameters.s("fmt2");
  std::string fmt3 = parameters.s("fmt3");
  std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
  std::string LECs = parameters.s("LECs");
  std::string scratch = parameters.s("scratch");
  std::string use_brueckner_bch = parameters.s("use_brueckner_bch");
  std::string valence_file_format = parameters.s("valence_file_format");
  std::string occ_file = parameters.s("occ_file");
  std::string goose_tank = parameters.s("goose_tank");
  std::string write_omega = parameters.s("write_omega");
  std::string nucleon_mass_correction = parameters.s("nucleon_mass_correction");
  std::string hunter_gatherer = parameters.s("hunter_gatherer");

  int eMax = parameters.i("emax");
  int E3max = parameters.i("e3max");
  int lmax3 = parameters.i("lmax3");
  int targetMass = parameters.i("A");
  int nsteps = parameters.i("nsteps");
  int file2e1max = parameters.i("file2e1max");
  int file2e2max = parameters.i("file2e2max");
  int file2lmax = parameters.i("file2lmax");
  int file3e1max = parameters.i("file3e1max");
  int file3e2max = parameters.i("file3e2max");
  int file3e3max = parameters.i("file3e3max");

  double hw = parameters.d("hw");
  double smax = parameters.d("smax");
  double ode_tolerance = parameters.d("ode_tolerance");
  double dsmax = parameters.d("dsmax");
  double ds_0 = parameters.d("ds_0");
  double domega = parameters.d("domega");
  double omega_norm_max = parameters.d("omega_norm_max");
  double denominator_delta = parameters.d("denominator_delta");
  double BetaCM = parameters.d("BetaCM");
  double hwBetaCM = parameters.d("hwBetaCM");
  double eta_criterion = parameters.d("eta_criterion");

  std::vector<std::string> opnames = parameters.v("Operators");
  std::vector<std::string> opsfromfile = parameters.v("OperatorsFromFile");
  std::vector<std::string> opnamesPT1 = parameters.v("OperatorsPT1");
  std::vector<std::string> opnamesRPA = parameters.v("OperatorsRPA");
  std::vector<std::string> opnamesTDA = parameters.v("OperatorsTDA");

  std::vector<Operator> ops;
  std::vector<std::string> spwf = parameters.v("SPWF");


  std::ifstream test;
  // test 2bme file
  if (inputtbme != "none")
  {
    test.open(inputtbme);
//    if( not test.good() and fmt2!="oakridge_binary")
    if( not test.good() and  fmt2.find("oakridge")== std::string::npos)
    {
      std::cout << "trouble reading " << inputtbme << " exiting. " << std::endl;
      return 1;
    }
    test.close();
  }
  // test 3bme file
  if (input3bme != "none")
  {
    test.open(input3bme);
    if( not test.good() )
    {
      std::cout << "trouble reading " << input3bme << " exiting. " << std::endl;
      return 1;
    }
    test.close();
  }



  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  rw.Set3NFormat( fmt3 );

  if(smryfile == "default") smryfile = parameters.DefaultSummaryFile();
  std::ofstream summary;
  summary.open(smryfile, std::ofstream::out);

//  ModelSpace modelspace;

  if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
  {
    if (valence_space=="") // if no name is given, then just name it "custom"
    {
      parameters.string_par["valence_space"] = "custom";
      flowfile = parameters.DefaultFlowFile();
      if(flow1file != "") flow1file = parameters.DefaultFlow1File();
      if(flow2file != "") flow2file = parameters.DefaultFlow2File();
      intfile = parameters.DefaultIntFile();
    }
    valence_space = custom_valence_space;
  }


  ModelSpace modelspace = ( reference=="default" ? ModelSpace(eMax,valence_space) : ModelSpace(eMax,reference,valence_space) );

  if (occ_file != "none" and occ_file != "" )
  {
    modelspace.Init_occ_from_file(eMax,valence_space,occ_file);
  }

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  modelspace.SetE3max(E3max);
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);



// For both dagger operators and single particle wave functions, it's convenient to
// just get every orbit in the valence space. So if SPWF="valence" ,  we append all valence orbits
  if ( std::find( spwf.begin(), spwf.end(), "valence" ) != spwf.end() )
  {
    // this erase/remove idiom is needed because remove just shuffles things around rather than actually removing it.
    spwf.erase( std::remove( spwf.begin(), spwf.end(), "valence" ), std::end(spwf) );
    for ( auto v : modelspace.valence )
    {
      spwf.push_back( modelspace.Index2String(v) );
    }
  }

  if ( std::find( opnames.begin(), opnames.end(), "DaggerHF_valence") != opnames.end() )
  {
    opnames.erase( std::remove( opnames.begin(), opnames.end(), "DaggerHF_valence"), std::end(opnames) );
    for ( auto v : modelspace.valence )
    {
      opnames.push_back( "DaggerHF_"+modelspace.Index2String(v) );
    }
    std::cout << "I found DaggerHF_valence, so I'm changing the opnames list to :" << std::endl;
    for ( auto opn : opnames ) std::cout << opn << " ,  ";
    std::cout << std::endl;
  }




//  std::cout << "Making the Hamiltonian..." << std::endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Hbare.SetHermitian();


  if ( goose_tank == "true" or goose_tank == "True")
  {
//    Hbare.SetUseGooseTank(true);
    Commutator::SetUseGooseTank(true);
  }
  std::cout << "Reading interactions..." << std::endl;
  if (inputtbme != "none")
  {
    if (fmt2 == "me2j")
      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
    else if (fmt2 == "navratil" or fmt2 == "Navratil")
      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
    else if (fmt2 == "oslo" )
      rw.ReadTBME_Oslo(inputtbme, Hbare);
    else if (fmt2.find("oakridge") != std::string::npos )
    { // input format should be: singleparticle.dat,vnn.dat
      size_t comma_pos = inputtbme.find_first_of(",");
      if ( fmt2.find("bin") != std::string::npos )
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "binary");
      else
        rw.ReadTBME_OakRidge( inputtbme.substr(0,comma_pos),  inputtbme.substr( comma_pos+1 ), Hbare, "ascii");
    }
    else if (fmt2 == "takayuki" )
      rw.ReadTwoBody_Takayuki( inputtbme, Hbare);
    else if (fmt2 == "nushellx" )
      rw.ReadNuShellX_int( Hbare, inputtbme );
    else if (fmt2 == "tokyo" )
      rw.ReadTokyo(inputtbme, Hbare);

    std::cout << "done reading 2N" << std::endl;
  }

  if (fmt2 != "nushellx")  // Don't need to add kinetic energy if we read a shell model interaction
  {
    Hbare += imsrg_util::Trel_Op(modelspace);
    if (Hbare.OneBody.has_nan())
    {
       std::cout << "  Looks like the Trel op is hosed from the get go." << std::endl;
    }
  }

  if ( nucleon_mass_correction == "true" or nucleon_mass_correction == "True" )
  {  // correction to kinetic energy because M_proton != M_neutron
    Hbare += imsrg_util::Trel_Masscorrection_Op(modelspace);
  }

  if (Hbare.particle_rank >=3)
  {
    rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
    std::cout << "done reading 3N" << std::endl;
  }


  // Add a Lawson term. If hwBetaCM is specified, use that frequency
  if (std::abs(BetaCM)>1e-3)
  {
    if (hwBetaCM < 0) hwBetaCM = modelspace.GetHbarOmega();
    std::ostringstream hcm_opname;
    hcm_opname << "HCM_" << hwBetaCM;
    Hbare += BetaCM * imsrg_util::OperatorFromString( modelspace, hcm_opname.str());
  }

  std::cout << "Creating HF" << std::endl;
  //HartreeFock hf(Hbare);
  HFMBPT hf(Hbare);
  std::cout << "Solving" << std::endl;
  hf.Solve();
  summary << "Ecore HF " << hf.EHF << std::endl;

//  Operator HNO;
  Operator& HNO = Hbare;
  if (basis == "HF" and method !="HF")
    HNO = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    HNO = Hbare.DoNormalOrdering();
  else if (basis == "NAT")
  {
    hf.GetNaturalOrbital();
    HNO = hf.GetNormalOrderedHNAT();
  }

  if ( spwf.size() > 0 )
  {
    imsrg_util::WriteSPWaveFunctions( spwf, hf, intfile);
  }

  HNO -= BetaCM * 1.5*hwBetaCM;
  std::cout << "Hbare 0b = " << HNO.ZeroBody << std::endl;
  // test
  //HNO = HNO.UndoNormalOrdering();
  //modelspace.SetReference(modelspace.core);
  //HNO.SetModelSpace(modelspace);
  //HNO = HNO.DoNormalOrdering();
  //rw.WriteTokyo(HNO,intfile+".snt");
  //exit(0);
  //

  //
  //exit(0); 
  // for check

  
  if (method != "HF")
  {
    std::cout << "Perturbative estimates of gs energy:" << std::endl;
    double EMP2 = HNO.GetMP2_Energy();
    std::cout << "EMP2 = " << EMP2 << std::endl;
//    double EMP3 = HNO.GetMP3_Energy();
    std::array<double,3> Emp_3 = HNO.GetMP3_Energy();
    double EMP3 = Emp_3[0]+Emp_3[1]+Emp_3[2];
    std::cout << "E3_pp = " << Emp_3[0] << "  E3_hh = " << Emp_3[1] << " E3_ph = " << Emp_3[2] << "   EMP3 = " << EMP3 << std::endl;
//    cout << "EMP3 = " << EMP3 << endl;
    std::cout << "To 3rd order, E = " << HNO.ZeroBody+EMP2+EMP3 << std::endl;
  }

  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
    ops.emplace_back( imsrg_util::OperatorFromString(modelspace,opname) );
  }
  // Calculate first order perturbative correction to some operators, if that's what we asked for.
  // Strictly speaking, it doesn't make much sense to do this and then proceed with the IMSRG calculation,
  // but I'm not here to tell people what to do...
  for (auto& opnamept1 : opnamesPT1 )
  {
    ops.emplace_back( imsrg_util::FirstOrderCorr_1b( imsrg_util::OperatorFromString(modelspace,opnamept1)   , HNO ) );
    opnames.push_back( opnamept1+"PT1" );
  }
  for (auto& opnametda : opnamesTDA )
  {  // passing the argument "TDA" just sets the phhp and hpph blocks to zero in the RPA calculation
    ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnametda)   , HNO, "TDA" ) );
    opnames.push_back( opnametda+"TDA" );
  }
  for (auto& opnamerpa : opnamesRPA )
  {
    ops.emplace_back( imsrg_util::RPA_resummed_1b( imsrg_util::OperatorFromString(modelspace,opnamerpa)   , HNO, "RPA" ) );
    opnames.push_back( opnamerpa+"RPA" );
  }


  // the format should look like OpName^j_t_p_r^/path/to/file
  for (auto& tag : opsfromfile)
  {
    std::istringstream ss(tag);
    std::string opname,qnumbers,fname;
    std::vector<int> qn(4);

    getline(ss,opname,'^');
    getline(ss,qnumbers,'^');
    getline(ss,fname,'^');
    ss.str(qnumbers);
    ss.clear();
    for (int i=0;i<4;i++)
    {
      std::string tmp;
      getline(ss,tmp,'_');
      std::istringstream(tmp) >> qn[i];
    }

    int j,t,p,r;
    j = qn[0];
    t = qn[1];
    p = qn[2];
    r = qn[3];
//    std::cout << "Parsed tag. opname = " << opname << "  qnumbers = " << qnumbers << "  " << j << " " << t << " " << p << " " << r << "   file = " << fname << std::endl;
    Operator op(modelspace,j,t,p,r);
    rw.Read2bCurrent_Navratil( fname, op );
    ops.push_back( op );
    opnames.push_back( opname );
  }



//  for (auto& op : ops)
  for (size_t i=0;i<ops.size();++i)
  {
     // We don't transform a DaggerHF, because we want the a^dagger to already refer to the HF basis.
    if ((basis == "HF") and (opnames[i].find("DaggerHF") == std::string::npos)  )
    {
      ops[i] = hf.TransformToHFBasis(ops[i]);
    }
    if (basis == "NAT") ops[i] = hf.TransformHOToNATBasis(ops[i]);
    ops[i] = ops[i].DoNormalOrdering();
    if (method == "MP3")
    {
      double dop = ops[i].MP1_Eval( HNO );
      std::cout << "Operator 1st order correction  " << dop << "  ->  " << ops[i].ZeroBody + dop << std::endl;
    }
  }


  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    std::cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << std::endl;
    std::cout << " HF charge radius = " << ( abs(Rp2.ZeroBody)<1e-6 ? 0.0 : sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) ) << std::endl;
  }
  for (index_t i=0;i<ops.size();++i)
  {
    std::cout << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
    summary << opnames[i] << "core " << basis << " " <<
      std::setprecision(8) << ops[i].ZeroBody << std::endl;
  }


  std::cout << "HF Single particle energies:" << std::endl;
  hf.PrintSPE();
  std::cout << std::endl;

  if ( method == "HF" or method == "MP3")
  {
    HNO.PrintTimes();
    return 0;
  }


  if (method == "FCI")
  {
    HNO = HNO.UndoNormalOrdering();
    rw.WriteNuShellX_int(HNO,intfile+".int");
    rw.WriteNuShellX_sps(HNO,intfile+".sp");

    for (index_t i=0;i<ops.size();++i)
    {
      ops[i] = ops[i].UndoNormalOrdering();
      if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
      {
        rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
      }
      else
      {
        rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
        rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
      }
    }
    HNO.PrintTimes();
    return 0;
  }

//  Operator HlowT = HNO;
//  double Temp = hw;
//  double Efermi = 0;
//  Operator Eye = HNO;
//  Eye.Eye();
//  HlowT.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  Eye.ScaleFermiDirac(HNO, Temp, Efermi);  // 0 is roughly fermi surface? we can do beter...
//  std::cout << "Initial low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) <<"  with normalization  " << Eye.Trace( modelspace.GetAref(),modelspace.GetZref() ) << std::endl;

  IMSRGSolver imsrgsolver(HNO);
  imsrgsolver.SetReadWrite(rw);
  imsrgsolver.SetEtaCriterion(eta_criterion);
  bool brueckner_restart = false;
  if (hunter_gatherer=="true") imsrgsolver.SetHunterGatherer( true);

  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=50000;
    method = "magnus";
  }
  if (method.find("brueckner") != std::string::npos)
  {
    if (method=="brueckner2") brueckner_restart=true;
    if (method=="brueckner1step")
    {
       nsteps = 1;
       core_generator = valence_generator;
    }
    use_brueckner_bch = "true";
    omega_norm_max=500;
    method = "magnus";
  }

  if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
  {
//    Hbare.SetUseBruecknerBCH(true);
//    HNO.SetUseBruecknerBCH(true);
    Commutator::SetUseBruecknerBCH(true);
    std::cout << "Using Brueckner flavor of BCH" << std::endl;
  }

  imsrgsolver.SetMethod(method);
//  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetHin(HNO);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetFlow1File(flow1file);
  imsrgsolver.SetFlow2File(flow2file);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDsmax(dsmax);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  imsrgsolver.SetGenerator(core_generator);
  if (core_generator.find("imaginary")!=std::string::npos)
  {
   if (ds_0>1e-2)
   {
     ds_0 = 1e-4;
     dsmax = 1e-2;
     imsrgsolver.SetDs(ds_0);
     imsrgsolver.SetDsmax(dsmax);
   }
  }
  imsrgsolver.Solve();

//  HlowT = imsrgsolver.Transform(HlowT);
//  std::cout << "After Solve, low temp trace with T = " << Temp << " and Ef = " << Efermi << ":   " << HlowT.Trace(modelspace.GetAref(),modelspace.GetZref()) << std::endl;

  if (method == "magnus")
  {
//    for (size_t i=0;i<ops.size();++i)
//    {
//      Operator tmp = imsrgsolver.Transform(ops[i]);
////      rw.WriteOperatorHuman(tmp,intfile+opnames[i]+"_step1.op");
//    }
//    std::cout << std::endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  if (brueckner_restart)
  {
     arma::mat newC = hf.C * arma::expmat( -imsrgsolver.GetOmega(0).OneBody  );
//     if (input3bme != "none") Hbare.SetParticleRank(3);
     HNO = hf.GetNormalOrderedH(newC);
     imsrgsolver.SetHin(HNO);
     imsrgsolver.s = 0;
     imsrgsolver.Solve();
  }

  if (nsteps > 1) // two-step decoupling, do core first
  {
    if (method == "magnus") smax *= 2;
    if (denominator_delta_orbit != "none")
      imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

    imsrgsolver.SetGenerator(valence_generator);
    modelspace.ResetFirstPass();
    if (valence_generator.find("imaginary")!=std::string::npos)
    {
     if (ds_0>1e-2)
     {
       ds_0 = 1e-4;
       dsmax = 1e-2;
       imsrgsolver.SetDs(ds_0);
       imsrgsolver.SetDsmax(dsmax);
     }
    }
    imsrgsolver.SetSmax(smax);
    imsrgsolver.Solve();
  }



  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) std::cout << "transforming operators" << std::endl;
    for (size_t i=0;i<ops.size();++i)
    {
      std::cout << opnames[i] << " " << std::endl;
      ops[i] = imsrgsolver.Transform(ops[i]);
      std::cout << " (" << ops[i].ZeroBody << " ) " << std::endl;
//      rw.WriteOperatorHuman(ops[i],intfile+opnames[i]+"_step2.op");
    }
    std::cout << std::endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted/ensemble normal ordering
  // we now re-normal order wrt to the core
  // and do any remaining flow.
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
  {
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (not renormal_order)
    {
      for (auto c : modelspace.core)
      {
         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (std::abs(1-modelspace.GetOrbit(c).occ)>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
  if ( renormal_order )
  {

    HNO = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
    std::cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << std::endl;
    HNO = HNO.UndoNormalOrdering();
    //
    imsrgsolver.SetHin(HNO);
    imsrgsolver.WriteStatusFlow1(imsrgsolver.flow1file);
    //

    ms2.SetReference(ms2.core); // change the reference
    HNO.SetModelSpace(ms2);

    std::cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << std::endl;
    HNO = HNO.DoNormalOrdering();

// More flowing is unnecessary, since things should stay decoupled.
    imsrgsolver.SetHin(HNO);
    //
    imsrgsolver.SetHin(HNO);
    imsrgsolver.WriteStatusFlow1(imsrgsolver.flow1file);
    //
//    imsrgsolver.SetEtaCriterion(1e-4);
//    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    std::cout << "Final transformation on the operators..." << std::endl;
    int iop = 0;
    for (auto& op : ops)
    {
      std::cout << opnames[iop++] << std::endl;
      op = op.UndoNormalOrdering();
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    if (valence_file_format == "antoine") // this is still being tested...
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
    }
    std::cout << "Writing files: " << intfile << std::endl;
    rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
    rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
    rw.WriteTokyo(imsrgsolver.GetH_s(),intfile+".snt","");
    //rw.WriteTokyoFull(imsrgsolver.GetH_s(),intfile+".snt");

    if (method == "magnus")
    {
       for (index_t i=0;i<ops.size();++i)
       {
          if ( ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1) and (ops[i].GetNumberLegs()%2==0) )
          {
            rw.WriteNuShellX_op(ops[i],intfile+"_"+opnames[i]+".int");
            rw.WriteTokyo(ops[i],intfile+"_"+opnames[i]+".snt","op");
          }
          else if ( ops[i].GetNumberLegs()%2==1) // odd number of legs -> this is a dagger operator
          {
//            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int"); // do this for now. later make a *.dag format.
            rw.WriteDaggerOperator( ops[i], intfile+opnames[i]+".dag",opnames[i]);
          }
          else
          {
            rw.WriteTensorOneBody(intfile+"_"+opnames[i]+"_1b.op",ops[i],opnames[i]);
            rw.WriteTensorTwoBody(intfile+"_"+opnames[i]+"_2b.op",ops[i],opnames[i]);
            rw.WriteTensorTokyo(intfile+"_"+opnames[i]+".snt",ops[i]);
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    std::cout << "Core Energy = " << std::setprecision(6) << imsrgsolver.GetH_s().ZeroBody << std::endl;
    summary << "Ecore IMSRG " << std::setprecision(8) << imsrgsolver.GetH_s().ZeroBody << std::endl;
    for (index_t i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      std::cout << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
      summary << opnames[i] << " = " << ops[i].ZeroBody << std::endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         std::cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
         std::cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << std::endl;
         summary << "IMSRG point proton radius = " << sqrt( op.ZeroBody ) << std::endl;
         summary << "IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << std::endl;
      }
      if ((op.GetJRank()>0) or (op.GetTRank()>0)) // if it's a tensor, you probably want the full operator
      {
        std::cout << "Writing operator to " << intfile+opnames[i]+".op" << std::endl;
        rw.WriteOperatorHuman(op,intfile+opnames[i]+".op");
      }
    }
  }


//  std::cout << "Made it here and write_omega is " << write_omega << std::endl;
  if (write_omega == "true" or write_omega == "True")
  {
    std::cout << "writing Omega to " << intfile << "_omega.op" << std::endl;
    rw.WriteOperatorHuman(imsrgsolver.Omega.back(),intfile+"_omega.op");
  }

  summary.close();
  Hbare.PrintTimes();

  return 0;
}

