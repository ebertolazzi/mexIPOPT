/*--------------------------------------------------------------------------*\
 |  file: IpoptInterfaceCommon.hh                                           |
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

#ifndef IPOPT_INTERFACE_COMMON_HH
#define IPOPT_INTERFACE_COMMON_HH

#pragma once

#include <stdexcept> // for unix system
#include <exception>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>

// STL
#include <algorithm>
#include <vector>
#include <iterator>

#include "mex.h"

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wsuggest-destructor-override"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wsuggest-override"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wcomma"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include "IpUtils.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

#if IPOPT_VERSION_MAJOR < 3
  #error "Ipopt Matlab Interface need IPOPT version >= 3.11"
#endif
#if IPOPT_VERSION_MINOR < 11
  #error "Ipopt Matlab Interface need IPOPT version >= 3.11"
#endif

#ifndef IPOPT_ASSERT
  #define IPOPT_DO_ERROR(MSG) \
    { std::ostringstream ost; ost << MSG; throw std::runtime_error(ost.str()); }
  #define IPOPT_ASSERT(COND,MSG) if ( !(COND) ) IPOPT_DO_ERROR(MSG)
#endif

//#define DEBUG
#ifdef DEBUG
  #define IPOPT_DEBUG(MSG) {             \
    std::ostringstream ost; ost << MSG;  \
    mexPrintf("%s\n",ost.str().c_str()); \
    mexEvalString("drawnow;");           \
  }
#else
  #define IPOPT_DEBUG(MSG)
#endif

// if C++ < C++11 define nullptr
#if __cplusplus <= 199711L
  #ifndef nullptr
    #include <cstddef>
    #define nullptr NULL
  #endif
#endif

#ifdef IPOPT_INTERFACE_MISSING_COPY_N
namespace std {
  template <typename Ta, typename Tb>
  inline
  void
  copy_n( Ta from, int n, Tb to )
  { std::copy( from, from+n, to ); }
}
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace IpoptInterface {

  using Ipopt::ApplicationReturnStatus;
  using Ipopt::TNLP;
  using Ipopt::SolverReturn;
  using Ipopt::AlgorithmMode;
  using Ipopt::IpoptData;
  using Ipopt::IpoptCalculatedQuantities;

  // Type definitions.
  // -----------------------------------------------------------------
  // This line is needed for versions of MATLAB prior to 7.3.
  #ifdef MWINDEXISINT
  typedef int mwIndex;
  #endif

  using Ipopt::Index;
  using Ipopt::Number;

  // Class MatlabFunctionHandle.
  // -----------------------------------------------------------------
  // The purpose of this class is twofold. The first aim is to store
  // information about a MATLAB function handle. (For more information
  // on function handles in MATLAB, type HELP FUNCTION_HANDLE in the
  // MATLAB console). The second purpose is to provide a routine for
  // evaluating the response of the function, provided inputs to the
  // function.
  class MatlabFunctionHandle {
  protected:
    mxArray * m_f = nullptr;  // The MATLAB function handle.
    std::string m_name;
  public:

    // The default constructor creates a null function handle.
    MatlabFunctionHandle( std::string const & name ) : m_name(name) { }

    // The destructor.
    ~MatlabFunctionHandle()
    { if ( m_f != nullptr ) mxDestroyArray( m_f ); m_f = nullptr; }

    std::string const & name() const { return m_name; }

    // This constructor accepts as input a pointer to a MATLAB array.
    // It is up to the user to ensure that the MATLAB array is a valid
    // function handle.
    void
    bind( mxArray const * p, char const error_msg[] );

    // This method is used to call the MATLAB function, provided inputs
    // to the function. It is up to the user to make sure that the
    // outputs array is of the appropriate size. It is up to the user to
    // properly deallocate the outputs.
    void
    eval(
      Index            n_lhs,
      mxArray       ** lhs,
      Index            n_rhs,
      mxArray const ** rhs
    ) const;

    // Returns true if and only if the function handle is not null.
    bool ok() const { return m_f != nullptr; }
  };

  /*
  //   ____
  //  / ___| _ __   __ _ _ __ ___  ___
  //  \___ \| '_ \ / _` | '__/ __|/ _ \
  //   ___) | |_) | (_| | |  \__ \  __/
  //  |____/| .__/ \__,_|_|  |___/\___|
  //        |_|
  */

  typedef struct SparseMatrix {
    Index              m_nnz;     // The height of the matrix.
    Index              m_numRows; // The height of the matrix.
    Index              m_numCols; // The width of the matrix.
    std::vector<Index> m_Jc;      // See mxSetJc in the MATLAB documentation.
    std::vector<Index> m_Ir;      // See mxSetIr in the MATLAB documentation.

    void setup( mxArray * ptr );
    void getStructure( Index rows[], Index cols[] ) const;
    void getValues( std::string const & func, mxArray * ptr, Number values[] ) const;
  } SparseMatrix;

  /*
  //    ____      _ _ _                _
  //   / ___|__ _| | | |__   __ _  ___| | __
  //  | |   / _` | | | '_ \ / _` |/ __| |/ /
  //  | |__| (_| | | | |_) | (_| | (__|   <
  //   \____\__,_|_|_|_.__/ \__,_|\___|_|\_\
  */

  // Class CallbackFunctions.
  // -----------------------------------------------------------------
  // An object of this class does two things. First of all, it stores
  // handles to MATLAB functions (type HELP FUNCTION_HANDLE in the
  // MATLAB console) for all the necessary and optional callback
  // functions for IPOPT. Secondly, this class actually provides the
  // routines for calling these functions with the necessary inputs and
  // outputs.
  class CallbackFunctions {
    MatlabFunctionHandle m_obj{"obj"};                     // Objective callback function.
    MatlabFunctionHandle m_grad{"grad"};                   // Gradient callback function.
    MatlabFunctionHandle m_constraint{"constraint"};       // Constraint callback function.
    MatlabFunctionHandle m_jacobian{"jacobian"};           // Jacobian callback function.
    MatlabFunctionHandle m_jacstruc{"jacobian_structure"}; // Jacobian structure function.
    MatlabFunctionHandle m_hessian{"hessian"};             // Hessian callback function.
    MatlabFunctionHandle m_hesstruc{"hessian_structure"};  // Hessian structure function.
    MatlabFunctionHandle m_iter{"iter"};                   // Iterative callback function.

    Index     mx_x_nc;
    Index     mx_x_nv;
    mxArray * mx_x;
    bool      m_x_is_cell_array; // true if x is a cell array

    std::vector<Number> m_x0;

    mutable SparseMatrix m_Jacobian;
    mutable SparseMatrix m_Hessian;

  public:

    // The constructor must be provided with a single MATLAB array, and
    // this MATLAB array must be a structure array in which each field
    // is a function handle.
    explicit
    CallbackFunctions( mxArray const * mx_x0, mxArray const * ptr );

    // The destructor.
    ~CallbackFunctions();

    void fillx( Index n, Number const x[] ) const;
    mxArray const * mx_getx() const { return mx_x; }

    std::vector<Number> const & getx0() const { return m_x0; }

    bool
    from_cell_array( mxArray const * ptr, Index n, Number * x ) const;

    Index numVariables() const { return mx_x_nv; }

    // These functions return true if the respective callback functions
    // are available.
    bool constraintFuncIsAvailable () const { return m_constraint.ok(); }
    bool jacobianFuncIsAvailable   () const { return m_jacobian.ok(); }
    bool hessianFuncIsAvailable    () const { return m_hessian.ok(); }
    bool iterFuncIsAvailable       () const { return m_iter.ok(); }

    // These functions execute the various callback functions with the
    // appropriate inputs and outputs. Here, m is the number of constraints.
    // The first function returns the value of the objective at x.
    Number computeObjective( Index m, Number const x[] ) const;

    // This function computes the value of the gradient at x, and
    // returns the gradient entries in the array g, which must be of
    // length equal to the number of optimization variables.
    void computeGradient( Index m, Number const x[], Number g[] ) const;

    // This function computes the response of the vector-valued
    // constraint function at x, and stores the result in the array c
    // which must be of length m.
    void computeConstraints( Index n, Number const x[], Index m, Number c[] ) const;

    // This function gets the structure of the sparse m x n Jacobian matrix.
    Index
    getJacobianNnz( ) const
    { return m_Jacobian.m_nnz; }

    void
    loadJacobianStructure( Index n, Index m ) const;

    void
    getJacobianStructure(  Index rows[], Index cols[] ) const
    { m_Jacobian.getStructure( rows, cols ); }

    // This function gets the structure of the sparse n x n Hessian matrix.
    Index
    getHessianNnz( ) const
    { return m_Hessian.m_nnz; }

    void
    loadHessianStructure( Index n ) const;

    void
    getHessianStructure( Index rows[], Index cols[] ) const
    { m_Hessian.getStructure( rows, cols ); }

    // This function computes the Jacobian of the constraints at x.
    void
    computeJacobian(
      Index        m,
      Index        n,
      Number const x[],
      Number       values[]
    ) const;

    // This function computes the Hessian of the Lagrangian at x.
    void
    computeHessian(
      Index        n,
      Number const x[],
      Number       sigma,
      Index        m,
      Number const lambda[],
      Number       values[]
    ) const;

    // Call the intermediate callback function. A return value of false
    // tells IPOPT to terminate.
    bool
    iterCallback(
      Index  iter,
      Number f,
      Number inf_pr,
      Number inf_du,
      Number mu,
      Number d_norm,
      Number regularization_size,
      Number alpha_du,
      Number alpha_pr,
      Index  ls_trials,
      Ipopt::IpoptData const           * ip_data,
      Ipopt::IpoptCalculatedQuantities * ip_cq
    ) const;
  };

  /*
  //   ___        __
  //  |_ _|_ __  / _| ___
  //   | || '_ \| |_ / _ \
  //   | || | | |  _| (_) |
  //  |___|_| |_|_|  \___/
  */

  // -----------------------------------------------------------------
  // An object of this class stores all the information we will pass
  // back to MATLAB upon termination of IPOPT.
  class MatlabInfo {
  protected:
    mxArray * m_info_ptr;  // All the information is stored in a MATLAB array.

  public:

    // Create a new info object and store the information in a MATLAB
    // array. The input pointer will point to the the newly created
    // MATLAB array. Since the user has an opportunity to modify the
    // MATLAB array pointed to by "ptr", we do not destroy the array
    // when the MatlabInfo object is destroyed. It is up to the user to
    // do that.
    explicit MatlabInfo( mxArray *& ptr );

    // The destructor.
    ~MatlabInfo();

    void setfield( mxArray const * ptr, char const * field );
    void setfield( size_t n, Number const * x, char const * field );
    void setfield( Number x, char const * field ) { setfield( 1, &x, field ); }

    mxArray *
    getfield_mx( char const * field ) const;

    Number const *
    getfield( char const * field ) const;

    // Access and modify the exit status and solution statistics.
    ApplicationReturnStatus getExitStatus() const;
    void setExitStatus( ApplicationReturnStatus status ) { setfield( (Number)status, "status"); }
    void setFuncEvals( Index obj, Index con, Index grad, Index jac, Index hess );
    void setIterationCount( Index iter ) { setfield( iter, "iter"); }
    void setCpuTime( Number cpu ) { setfield( cpu, "cpu"); }
  };

  /*
  //    ___        _   _
  //   / _ \ _ __ | |_(_) ___  _ __  ___
  //  | | | | '_ \| __| |/ _ \| '_ \/ __|
  //  | |_| | |_) | |_| | (_) | | | \__ \
  //   \___/| .__/ \__|_|\___/|_| |_|___/
  //        |_|
  */
  /*
  // Class IpoptOptions.
  // -----------------------------------------------------------------
  // This class processes the IPOPT options as specified by a user in the
  // MATLAB environment.
  */
  class IpoptOptions {
  protected:
    Ipopt::IpoptApplication & m_app;  // The IPOPT application object.

    // These three functions are used by the class constructor.
    void setOption        ( char const * label, mxArray const * ptr );
    void setStringOption  ( char const * label, mxArray const * ptr );
    void setIntegerOption ( char const * label, mxArray const * ptr );
    void setNumberOption  ( char const * label, mxArray const * ptr );

  public:

    // The constructor accepts as input an IPOPT application object and
    // a MATLAB array. The latter input must be a structure array, with
    // field names corresponding to the names of options in IPOPT.
    IpoptOptions( Ipopt::IpoptApplication & app, mxArray const * ptr );

    // The destructor.
    ~IpoptOptions();

    // The first function returns true if and only if the user has
    // specified a quasi-Newton approximation to the Hessian instead of
    // the exact Hessian. The second function returns true if and only
    // if the user has activated the derivative checker. The third
    // function returns true if and only if a user-specified scaling of
    // the problem is activated. The fourth function returns the print
    // level for the IPOPT console. The remaining two functions return
    // the floating-point value for positive and negative infinity,
    // respectively.
    bool useQuasiNewton () const;
    bool useDerivChecker() const;
    bool userScaling    () const;
    Index printLevel() const;
  };

  // Class Options.
  // -----------------------------------------------------------------
  // This class processes the options input from MATLAB.
  class Options {
  protected:
    Index               m_n;       // The number of optimization variables.
    Index               m_m;       // The number of constraints.
    std::vector<Number> m_lb;      // Lower bounds on the variables.
    std::vector<Number> m_ub;      // Upper bounds on the variables.
    std::vector<Number> m_cl;      // Lower bounds on constraints.
    std::vector<Number> m_cu;      // Upper bounds on constraints.
    std::vector<Number> m_zl;      // Lagrange multipliers for lower bounds.
    std::vector<Number> m_zu;      // Lagrange multipliers for upper bounds.
    std::vector<Number> m_lambda;  // Lagrange multipliers for constraints.

    Number m_neginfty, m_posinfty;

    IpoptOptions m_ipopt;   // The IPOPT options.

    // These are helper functions used by the class constructor.
    void loadLowerBounds( mxArray const * ptr );
    void loadUpperBounds( mxArray const * ptr );
    void loadConstraintBounds( mxArray const * ptr );
    void loadMultipliers( mxArray const * pt );

  public:

    // The constructor expects as input a point to a MATLAB array, in
    // particular a structure array with the appropriate fields. Note
    // that the Options object does *not* possess an independent copy of
    // some of the MATLAB data (such as the auxiliary data).
    Options( Index n, Ipopt::IpoptApplication & app, mxArray const * ptr );

    // The destructor.
    ~Options();

    // Get the number of variables and the number of constraints.
    friend Index numvars        ( Options const & options ) { return options.m_n; }
    friend Index numconstraints ( Options const & options ) { return options.m_m; }

    // Access the lower and upper bounds on the variables and constraints.
    std::vector<Number> const & lowerbounds () const { return m_lb; }
    std::vector<Number> const & upperbounds () const { return m_ub; }
    std::vector<Number> const & constraintlb() const { return m_cl; }
    std::vector<Number> const & constraintub() const { return m_cu; }

    // Access the IPOPT options object.
    IpoptOptions const & ipoptOptions() const { return m_ipopt; }

    // Access the Lagrange multpliers.
    std::vector<Number> const & multlb    () const { return m_zl; }
    std::vector<Number> const & multub    () const { return m_zu; }
    std::vector<Number> const & multconstr() const { return m_lambda; }

  };

  /*
  //   ____
  //  |  _ \ _ __ ___   __ _ _ __ __ _ _ __ ___
  //  | |_) | '__/ _ \ / _` | '__/ _` | '_ ` _ \
  //  |  __/| | | (_) | (_| | | | (_| | | | | | |
  //  |_|   |_|  \___/ \__, |_|  \__,_|_| |_| |_|
  //                   |___/
  */
  // Class MatlabProgram
  // -----------------------------------------------------------------
  class MatlabProgram : public TNLP {
  protected:
    CallbackFunctions const & m_funcs;    // Callback routines.
    Options const &           m_options;  // Further program info.
    MatlabInfo &              m_info;     // Info passed back to MATLAB.

  public:

    // The constructor.
    MatlabProgram(
      CallbackFunctions const & funcs,
      Options           const & options,
      MatlabInfo              & info
    );

    // The destructor.
    virtual ~MatlabProgram() override;

    // Method to return some info about the nonlinear program.
    virtual
    bool
    get_nlp_info(
      Index          & n,
      Index          & m,
      Index          & sizeOfJ,
      Index          & sizeOfH,
      IndexStyleEnum & indexStyle
    ) override;

    // Return the bounds for the problem.
    virtual
    bool
    get_bounds_info(
      Index    n,
      Number * lb,
      Number * ub,
      Index    m,
      Number * cl,
      Number * cu
    ) override;

    // Return the starting point for the algorithm.
    virtual
    bool
    get_starting_point(
      Index    n,
      bool     initializeVars,
      Number * vars,
      bool     initializez,
      Number * zl,
      Number * zu,
      Index    m,
      bool     initializeLambda,
      Number * lambda
    ) override;

    // Compute the value of the objective.
    virtual
    bool
    eval_f(
      Index          n,
      Number const * vars,
      bool           ignore,
      Number &       f
    ) override;

    // Compute the gradient of the objective.
    virtual
    bool
    eval_grad_f(
      Index          n,
      Number const * vars,
      bool           ignore,
      Number *       grad
    ) override;

    // Evaluate the constraint residuals.
    virtual
    bool
    eval_g(
      Index          n,
      Number const * vars,
      bool           ignore,
      Index          m,
      Number *       g
    ) override;

    // This method either returns: 1.) The structure of the Jacobian
    // (if "Jacobian" is zero), or 2.) The values of the Jacobian (if
    // "Jacobian" is not zero).
    virtual
    bool
    eval_jac_g(
      Index          numVariables,
      Number const * variables,
      bool           ignoreThis,
      Index          numConstraints,
      Index          Jx_nnz,
      Index        * rows,
      Index        * cols,
      Number       * Jx
    ) override;

    // This method either returns: 1.) the structure of the Hessian of
    // the Lagrangian (if "Hessian" is zero), or 2.) the values of the
    // Hessian of the Lagrangian (if "Hesson" is not zero).
    virtual
    bool
    eval_h(
      Index          n,
      Number const * vars,
      bool           ignore,
      Number         sigma,
      Index          m,
      Number const * lambda,
      bool           ignoretoo,
      Index          Hx_nnz,
      Index        * rows,
      Index        * cols,
      Number       * Hx
    ) override;

    // This method is called when the algorithm is complete.
    virtual
    void
    finalize_solution(
      SolverReturn                status,
      Index                       numVariables,
      Number const              * variables,
      Number const              * zl,
      Number const              * zu,
      Index                       numConstraints,
      Number const              * constraints,
      Number const              * lambda,
      Number                      objective,
      IpoptData const           * ip_data,
      IpoptCalculatedQuantities * ip_cq
    ) override;

    // Intermediate callback method. It is called once per iteration
    // of the IPOPT algorithm.
    virtual
    bool
    intermediate_callback(
      AlgorithmMode               mode,
      Index                       t,
      Number                      f,
      Number                      inf_pr,
      Number                      inf_du,
      Number                      mu,
      Number                      d_norm,
      Number                      regularization_size,
      Number                      alpha_du,
      Number                      alpha_pr,
      Index                       ls_trials,
      IpoptData const           * ip_data,
      IpoptCalculatedQuantities * ip_cq
    ) override;

  };

}

/*
//       _                              _
//      | | ___  _   _ _ __ _ __   __ _| |
//   _  | |/ _ \| | | | '__| '_ \ / _` | |
//  | |_| | (_) | |_| | |  | | | | (_| | |
//   \___/ \___/ \__,_|_|  |_| |_|\__,_|_|
*/

#include "IpJournalist.hpp"

namespace Ipopt {

  // Class MatlabJournal.
  // ---------------------------------------------------------------
  // This class encapsulates journal output to the MATLAB console.
  class MatlabJournal : public Journal {
  public:

    // The constructor.
    MatlabJournal( EJournalLevel default_level );

    // The destructor.
    virtual ~MatlabJournal() override;

  protected:

    // These functions override the functions in the Journal class.
    virtual
    void
    PrintImpl(
      EJournalCategory category,
      EJournalLevel    level,
      char const *     str
    ) override;

    virtual
    std::string
    Name() override
    { return "MatlabJournal"; }

    virtual
    void
    PrintfImpl(
      EJournalCategory category,
      EJournalLevel    level,
      char const *     pformat,
      va_list          ap
    ) override;

    virtual
    void
    FlushBufferImpl() override
    {}

  };
}

#endif
