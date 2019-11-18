#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <omp.h>
#include "IMSRG.hh"
#include "PhysicalConstants.hh"
#include "Parameters.hh"
#include "Atom.hh"

namespace Atom
{
  int MainAtom(int argc, char** argv)
  {
    Parameters parameters(argc,argv);
    if (parameters.help_mode) return 0;

    std::string inputtbme = parameters.s("2bme");
    std::string reference = parameters.s("reference");
    std::string valence_space = parameters.s("valence_space");
    std::string custom_valence_space = parameters.s("custom_valence_space");
    std::string basis = parameters.s("basis");
    std::string basis_type = parameters.s("basis_type");
    std::string method = parameters.s("method");
    std::string flowfile = parameters.s("flowfile");
    std::string intfile = parameters.s("intfile");
    std::string core_generator = parameters.s("core_generator");
    std::string valence_generator = parameters.s("valence_generator");
    std::string fmt2 = parameters.s("fmt2");
    std::string denominator_delta_orbit = parameters.s("denominator_delta_orbit");
    std::string scratch = parameters.s("scratch");
    std::string occ_file = parameters.s("occ_file");
    std::string goose_tank = parameters.s("goose_tank");
    std::string IMSRG3 = parameters.s("IMSRG3");
    std::string freeze_occupations = parameters.s("freeze_occupations");
    bool use_NAT_occupations = (parameters.s("use_NAT_occupations")=="true") ? true : false;
    int eMax = parameters.i("emax");
    int lmax = parameters.i("lmax"); // so far I only use this with atomic systems.
    int file2e1max = parameters.i("file2e1max");
    int file2e2max = parameters.i("file2e2max");
    int file2lmax = parameters.i("file2lmax");
    int targetMass = parameters.i("A");
    int nsteps = parameters.i("nsteps");

    double hw = parameters.d("hw");
    double smax = parameters.d("smax");
    double ode_tolerance = parameters.d("ode_tolerance");
    double dsmax = parameters.d("dsmax");
    double ds_0 = parameters.d("ds_0");
    double domega = parameters.d("domega");
    double omega_norm_max = parameters.d("omega_norm_max");
    double denominator_delta = parameters.d("denominator_delta");
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
    if (inputtbme != "none" and fmt2.find("oakridge")==std::string::npos and fmt2 != "schematic" )
    {
      test.open(inputtbme);
      ////    if( not test.good() and fmt2!="oakridge_binary")
      //    if( not test.good() and  fmt2.find("oakridge")== std::string::npos)
      if( not test.good() )
      {
        std::cout << "trouble reading " << inputtbme << "  fmt2 = " << fmt2 << "   exiting. " << std::endl;
        return 1;
      }
      test.close();
    }

    ReadWrite rw;
    rw.SetScratchDir(scratch);
    // Test whether the scratch directory exists and we can write to it.
    // This is necessary because otherwise you get garbage for transformed operators and it's
    // not obvious what went wrong.
    if ( method=="magnus" and  scratch != "" and scratch!= "/dev/null" and scratch != "/dev/null/")
    {
      std::string testfilename = scratch + "/_this_is_a_test_delete_me";
      std::ofstream testout(testfilename);
      testout << "PASSED" << std::endl;
      testout.close();
      std::remove( testfilename.c_str() );
      if ( not testout.good() )
      {
        std::cout << "WARNING in " << __FILE__ <<  " failed test write to scratch directory " << scratch;
        if (opnames.size()>0 )
        {
          std::cout << "   dying now. " << std::endl;
          exit(EXIT_FAILURE);
        }
        else std::cout << std::endl;
      }
    }
    if ( (method=="magnus") and (scratch=="/dev/null" or scratch=="/dev/null/") )
    {
      if ( opnames.size() > 0 )
      {
        std::cout << "WARNING!!! using Magnus with scratch = " << scratch << " but you're also trying to transform some operators: ";
        for (auto opn : opnames ) std::cout << opn << " ";
        std::cout << "   dying now." << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    //  ModelSpace modelspace;

    if (custom_valence_space!="") // if a custom space is defined, the input valence_space is just used as a name
    {
      if (valence_space=="") // if no name is given, then just name it "custom"
      {
        parameters.string_par["valence_space"] = "custom";
        flowfile = parameters.DefaultFlowFile();
        intfile = parameters.DefaultIntFile();
      }
      valence_space = custom_valence_space;
    }

    ModelSpace modelspace;
    modelspace.SetLmax(lmax);
    modelspace.InitAtomicSpace(eMax, basis_type, reference, valence_space);
    if (nsteps < 0) nsteps = modelspace.valence.size()>0 ? 2 : 1;
    Operator Hbare = Operator(modelspace,0,0,0,2);
    Hbare.SetHermitian();
    if ( goose_tank == "true" or goose_tank == "True") Commutator::SetUseGooseTank(true);
    //
    std::cout << "Reading interactions..." << std::endl;
    if (inputtbme != "none")
    {
      if (fmt2 == "tokyo") rw.ReadTokyoAtomic(inputtbme,Hbare);
    }

    if (inputtbme == "none")
    {
      using PhysConst::M_ELECTRON;
      using PhysConst::M_NUCLEON;
      int Z = modelspace.GetTargetZ() ;
      Hbare -= Z*imsrg_util::VCentralCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;
      Hbare += imsrg_util::VCoulomb_Op(modelspace, lmax) * sqrt((M_ELECTRON*1e6)/M_NUCLEON ) ;  // convert oscillator length from fm with nucleon mass to nm with electon mass (in eV).
      Hbare += imsrg_util::KineticEnergy_Op(modelspace); // Don't need to rescale this, because it's related to the oscillator frequency, which we input.
      Hbare /= PhysConst::HARTREE; // Convert to Hartree
    }

    std::cout << "Creating HF" << std::endl;
    HFMBPT hf(Hbare); // HFMBPT inherits from HartreeFock, so no harm done here.
    if (freeze_occupations == "false" )  hf.UnFreezeOccupations();
    std::cout << "Solving" << std::endl;
    hf.Solve();
    return 0;

  }
}
