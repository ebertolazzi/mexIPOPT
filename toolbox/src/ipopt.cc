/*--------------------------------------------------------------------------*\
 |  file: Ipopt.cc                                                          |
 |                                                                          |
 |  IPOPT MATLAB-interface provided by:                                     |
 |      Enrico Bertolazzi (enrico.bertolazzi@unitn.it)                      |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      Via Sommarive 9, I-38123, Trento, Italy                             |
 |                                                                          |
 |  This IPOPT MATLAB-interface is derived from the code by:                |
 |        Peter Carbonetto                                                  |
 |        Dept. of Computer Science                                         |
 |        University of British Columbia, May 19, 2007                      |
\*--------------------------------------------------------------------------*/

/*
 * (C) Copyright 2015 Enrico Bertolazzi
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 */

#include "mex.h"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wignored-attributes"
#endif

#include "IpoptInterfaceCommon.hh"

#include "IpRegOptions.hpp"
#include "IpJournalist.hpp"
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"

#ifdef IPOPT_INTERFACE_USE_MPI
  #include <mpi.h>
#endif

using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::EJournalLevel;
using Ipopt::Journal;
using Ipopt::MatlabJournal;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::TNLP;
using Ipopt::ApplicationReturnStatus;
using Ipopt::SolveStatistics;

/*
// redirect stdout, found at
// https://it.mathworks.com/matlabcentral/answers/132527-in-mex-files-where-does-output-to-stdout-and-stderr-go
*/

#if defined(OS_LINUX) || defined(OS_MAC)
class mystream : public std::streambuf {
protected:
  virtual
  std::streamsize
  xsputn(const char *s, std::streamsize n) override
  { mexPrintf("%.*s", n, s); mexEvalString("drawnow;"); return n; }

  virtual
  int
  overflow(int c=EOF) override
  { if (c != EOF) { mexPrintf("%.1s", &c); mexEvalString("drawnow;"); } return 1; }

};

class scoped_redirect_cout {
public:
  scoped_redirect_cout()
  { old_buf = std::cout.rdbuf(); std::cout.rdbuf(&mout); }
  ~scoped_redirect_cout()
  { std::cout.rdbuf(old_buf); }
private:
  mystream mout;
  std::streambuf *old_buf;
};
static scoped_redirect_cout mycout_redirect;

#endif


namespace IpoptInterface {

  static
  void
  mexFunction_internal( int nlhs, mxArray       *plhs[],
	  	                  int nrhs, mxArray const *prhs[] ) {

    // Check to see if we have the correct number of input and output
    // arguments.
    IPOPT_ASSERT(
      nrhs == 3,
      "Incorrect number of input arguments, expected 3 found " << nrhs
    );
    IPOPT_ASSERT(
      nlhs == 2,
      "Incorrect number of output arguments, expected 2 found " << nlhs
    );

    IPOPT_DEBUG("\nCALL CallbackFunctions\n");
    // Get the second input which specifies the callback functions.
    CallbackFunctions funcs( prhs[0], prhs[1] );

    // Get the third input which specifies the options.
    // Create a new IPOPT application object and process the options.
    IPOPT_DEBUG("\nCALL app\n");
    IpoptApplication app(false);

    IPOPT_DEBUG("\nCALL options\n");
    Options options( funcs.numVariables(), app, prhs[2] ); // app, options

    #if IPOPT_VERSION_MINOR > 11
    app.RethrowNonIpoptException(true);
    #endif

    // The first output argument is the value of the optimization
    // variables obtained at the solution.

    // The second output argument stores other information, such as
    // the exit status, the value of the Lagrange multipliers upon
    // termination, the final state of the auxiliary data, and so on.
    IPOPT_DEBUG("\nCALL info\n");
    MatlabInfo info( plhs[1] );

    // Check to see whether the user provided a callback function for
    // computing the Hessian. This is not needed in the special case
    // when a quasi-Newton approximation to the Hessian is being used.
    IPOPT_ASSERT(
      options.ipoptOptions().useQuasiNewton() ||
      funcs.hessianFuncIsAvailable(),
      "You must supply a callback function for computing the Hessian "
      "unless you decide to use a quasi-Newton approximation to the Hessian"
    );

    // If the user tried to use her own scaling, report an error.
    IPOPT_ASSERT(
      !options.ipoptOptions().userScaling(),
      "The user-defined scaling option does not work in the MATLAB interface for IPOPT"
    );

    // If the user supplied initial values for the Lagrange
    // multipliers, activate the "warm start" option in IPOPT.
    if ( options.multlb().size()>0 &&
         options.multub().size()>0 &&
         (numconstraints(options)==0 || options.multconstr().size()>0) )
      app.Options()->SetStringValue("warm_start_init_point","yes");

    // Set up the IPOPT console.
    IPOPT_DEBUG("\nCALL printLevel\n");
    EJournalLevel printLevel = (EJournalLevel) options.ipoptOptions().printLevel();
    SmartPtr<Journal> console;
    if ( printLevel > 0 ) { //prevents IPOPT display if we don't want it
      IPOPT_DEBUG("\nCALL console\n");
      console = new MatlabJournal(printLevel);
      IPOPT_DEBUG("\nCALL journal\n");
      app.Jnlst()->AddJournal(console);
    }

    // Intialize the IpoptApplication object and process the options.
    IPOPT_DEBUG("\nCALL initialize\n");
    ApplicationReturnStatus exitstatus = app.Initialize();
    IPOPT_ASSERT(
      exitstatus == Ipopt::Solve_Succeeded,
      "IPOPT solver initialization failed"
    );

    // Create a new instance of the constrained, nonlinear program.
    IPOPT_DEBUG("\nCALL matlab\n");
    SmartPtr<TNLP> program = new MatlabProgram( funcs, options, info );

    // Ask Ipopt to solve the problem.
    IPOPT_DEBUG("\nCALL optimize\n");
    exitstatus = app.OptimizeTNLP(program);
    info.setExitStatus(exitstatus);
    plhs[0] = mxDuplicateArray(info.getfield_mx("x"));

    // Collect statistics about Ipopt run
    IPOPT_DEBUG("\nCALL statistic\n");
    if ( IsValid(app.Statistics()) ) {
      SmartPtr<SolveStatistics> stats = app.Statistics();
      info.setIterationCount(stats->IterationCount());
      //Get Function Calls
      Index obj, con, grad, jac, hess;
      stats->NumberOfEvaluations( obj, con, grad, jac, hess );
      info.setFuncEvals( obj, con, grad, jac, hess );
      //CPU Time
      info.setCpuTime( stats->TotalCpuTime() );
    }

    // Free the dynamically allocated memory.
    //mxDestroyArray(x0.mx_ptr());
  }

}

// Function definitions.
// -----------------------------------------------------------------
void
mexFunction(
  int nlhs, mxArray       *plhs[],
  int nrhs, mxArray const *prhs[]
) {

  #ifdef IPOPT_INTERFACE_USE_MPI
  // mumps use mpi
  int flag;
  int ok = MPI_Initialized(&flag);
  if ( !flag ) MPI_Init(nullptr,nullptr);
  #endif

  try {

    IpoptInterface::mexFunction_internal( nlhs, plhs, nrhs, prhs );

  } catch ( std::exception & error ) {

    std::ostringstream ost;
    ost << "\n*** Error using Ipopt Matlab interface: ***\n\n"
        << error.what()
        << "\n*******************************************\n";
    mexErrMsgTxt(ost.str().c_str());

  } catch ( ... ) {

    mexErrMsgTxt("\n*** Error using Ipopt Matlab interface: ***\nUnknown error\n");

  }
  //mexPrintf("\n*** IPOPT DONE ***\n");

}
