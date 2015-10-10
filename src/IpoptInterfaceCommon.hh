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
 |        Original code is published under the Eclipse Public License.      |
\*--------------------------------------------------------------------------*/

/*
//  Copyright (c) 2015, Enrico Bertolazzi
//  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//
//      * Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//      * Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in
//        the documentation and/or other materials provided with the distribution
//      * Neither the name of the  nor the names
//        of its contributors may be used to endorse or promote products derived
//        from this software without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
//  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef IPOPT_INTERFACE_COMMON_HH
#define IPOPT_INTERFACE_COMMON_HH

#include "mex.h"

#include "IpUtils.hpp"
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"

#include <stdexcept> // for unix system
#include <exception>
#include <cstdio>
#include <sstream>

// STL
#include <algorithm>
#include <vector>

#if IPOPT_VERSION_MAJOR < 3
  #error "Ipopt Matlab Interface need IPOPT version >= 3.11"
#endif
#if IPOPT_VERSION_MINOR < 11
  #error "Ipopt Matlab Interface need IPOPT version >= 3.11"
#endif

#ifndef IPOPT_ASSERT
  #define IPOPT_DO_ERROR(MSG) \
    { std::ostringstream ost ; ost << MSG ; throw std::runtime_error(ost.str()) ; }
  #define IPOPT_ASSERT(COND,MSG) if ( !(COND) ) IPOPT_DO_ERROR(MSG)
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
  { std::copy( from, from+n, to ) ; }
}
#endif

namespace IpoptInterface {

  using Ipopt::ApplicationReturnStatus ;
  using Ipopt::TNLP ;
  using Ipopt::SolverReturn ;
  using Ipopt::AlgorithmMode ;
  using Ipopt::IpoptData ;
  using Ipopt::IpoptCalculatedQuantities ;

  // Type definitions.
  // -----------------------------------------------------------------
  // This line is needed for versions of MATLAB prior to 7.3.
  #ifdef MWINDEXISINT
  typedef int mwIndex;
  #endif

  using Ipopt::Index ;
  using Ipopt::Number ;

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
    mxArray     * f;  // The MATLAB function handle.
    std::string f_name ;
  public:

    // The default constructor creates a null function handle.
    MatlabFunctionHandle() : f(nullptr) { }

    // The destructor.
    ~MatlabFunctionHandle() { if ( f != nullptr ) mxDestroyArray(f) ; }

    std::string const & name() const { return f_name ; }

    // This constructor accepts as input a pointer to a MATLAB array.
    // It is up to the user to ensure that the MATLAB array is a valid
    // function handle.
    void
    bind( mxArray const * p, char const * error_msg ) ;

    // This method is used to call the MATLAB function, provided inputs
    // to the function. It is up to the user to make sure that the
    // outputs array is of the appropriate size. It is up to the user to
    // properly deallocate the outputs.
    void
    eval( Index           n_lhs,
          mxArray       * lhs[],
          Index           n_rhs,
          mxArray const * rhs[] ) const ;

    // Returns true if and only if the function handle is not null.
    operator bool() const { return f != nullptr ; };
  };

  /*
  //   ____
  //  / ___| _ __   __ _ _ __ ___  ___
  //  \___ \| '_ \ / _` | '__/ __|/ _ \
  //   ___) | |_) | (_| | |  \__ \  __/
  //  |____/| .__/ \__,_|_|  |___/\___|
  //        |_|
  */

  typedef struct {
    Index              nnz     ; // The height of the matrix.
    Index              numRows ; // The height of the matrix.
    Index              numCols ; // The width of the matrix.
    std::vector<Index> Jc      ; // See mxSetJc in the MATLAB documentation.
    std::vector<Index> Ir      ; // See mxSetIr in the MATLAB documentation.

    // get the pattern if the sparse matrix in ptr and store in the structure.
    void
    setup( mxArray * ptr ) {
      // Get the height, width and number of non-zeros.
      numRows = (Index) mxGetM(ptr);
      numCols = (Index) mxGetN(ptr);
      // The precise number of non-zero elements is contained in the
      // last entry of the jc array. (There is one jc entry for each
      // column in the matrix, plus an extra one.)
      Jc.resize(numCols+1);
      std::copy_n(mxGetJc(ptr),numCols+1,Jc.begin());
      // Copy the row and column indices, and the values of the nonzero entries.
      nnz = Jc.back() ;
      Ir.resize(nnz);
      std::copy_n(mxGetIr(ptr),nnz,Ir.begin());
    }

    // extract the pattern of the sparse matrix in the ipopt format
    void
    getStructure( Index rows[], Index cols[] ) const {
      for ( Index c = 0, i = 0 ; c < numCols ; ++c ) {
        for ( ; i < Jc[c+1] ; ++i ) { cols[i] = c ; rows[i] = Ir[i] ; }
      }
    }

    // extyract values from the sparse matrix in ptr and store in values
    // which correspond to nonzeros of the sparse matrix described by
    // Ir and Jc. The nonzeros of sparse matrix in ptr must be a subset of the
    // nonzeros of the pattern described by Ir and Jc.
    void
    getValues( char const * func, mxArray * ptr, Number values[] ) const {
      // il patterm può essere un sottoinsieme
      mwIndex const * mxJc = mxGetJc(ptr) ;
      mwIndex const * mxIr = mxGetIr(ptr) ;
      double  const * v    = mxGetPr(ptr) ;

      std::fill_n(values,nnz,0) ;
      Index k = 0 ;
      for ( Index c = 0, i = 0 ; c < numCols ; ++c ) {
        for ( i = mxJc[c], k = Jc[c] ; i < mxJc[c+1] ; ++i, ++k ) {
          while ( Ir[k] < mxIr[i] && k < Jc[c+1] ) ++k ; // skip not set elements
          IPOPT_ASSERT( Ir[k] == mxIr[i],
                        "In MATLAB function " << func << "\nbad pattern" ) ;
          values[k] = v[i] ;
        }
      }
    }
  } SparseMatrix ;

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
    MatlabFunctionHandle obj_func ;        // Objective callback function.
    MatlabFunctionHandle grad_func ;       // Gradient callback function.
    MatlabFunctionHandle constraint_func ; // Constraint callback function.
    MatlabFunctionHandle jacobian_func ;   // Jacobian callback function.
    MatlabFunctionHandle jacstruc_func ;   // Jacobian structure function.
    MatlabFunctionHandle hessian_func ;    // Hessian callback function.
    MatlabFunctionHandle hesstruc_func ;   // Hessian structure function.
    MatlabFunctionHandle iter_func ;       // Iterative callback function.

    Index nc, nv ;
    mxArray * mx_x ;
    bool x_is_cell_array ; // true if x is a cell array

    std::vector<Number> x0 ;

    mutable SparseMatrix Jacobian ;
    mutable SparseMatrix Hessian ;

  public:

    // The constructor must be provided with a single MATLAB array, and
    // this MATLAB array must be a structure array in which each field
    // is a function handle.
    explicit CallbackFunctions( mxArray const * mx_x0, mxArray const * ptr );

    // The destructor.
    ~CallbackFunctions();

    void fillx( Index n, Number const x[] ) const ;
    mxArray const * mx_getx() const { return mx_x ; }

    std::vector<Number> const & getx0() const { return x0 ; }

    bool
    from_cell_array( mxArray const * ptr, Index n, Number * x ) const ;

    Index numVariables() const { return nv ; }

    // These functions return true if the respective callback functions
    // are available.
    bool constraintFuncIsAvailable () const { return constraint_func ; }
    bool jacobianFuncIsAvailable   () const { return jacobian_func ; }
    bool hessianFuncIsAvailable    () const { return hessian_func ; }
    bool iterFuncIsAvailable       () const { return iter_func ; }

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
    Index getJacobianNnz( ) const { return Jacobian.nnz ; }

    void loadJacobianStructure( Index n, Index m ) const;

    void
    getJacobianStructure(  Index rows[], Index cols[] ) const
    { Jacobian.getStructure(rows,cols) ; }

    // This function gets the structure of the sparse n x n Hessian matrix.
    Index getHessianNnz( ) const { return Hessian.nnz ; }

    void loadHessianStructure( Index n ) const;

    void
    getHessianStructure( Index rows[], Index cols[] ) const
    { Hessian.getStructure(rows,cols) ; }

    // This function computes the Jacobian of the constraints at x.
    void
    computeJacobian( Index        m,
                     Index        n,
                     Number const x[],
                     Number       values[] ) const;

    // This function computes the Hessian of the Lagrangian at x.
    void
    computeHessian( Index        n,
                    Number const x[],
                    Number       sigma,
                    Index        m,
                    Number const lambda[],
                    Number       values[] ) const;

    // Call the intermediate callback function. A return value of false
    // tells IPOPT to terminate.
    bool
    iterCallback ( Index  iter,
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
			             Ipopt::IpoptCalculatedQuantities * ip_cq ) const;
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
    mxArray * ptr;  // All the information is stored in a MATLAB array.

  public:

    // Create a new info object and store the information in a MATLAB
    // array. The input pointer will point to the the newly created
    // MATLAB array. Since the user has an opportunity to modify the
    // MATLAB array pointed to by "ptr", we do not destroy the array
    // when the MatlabInfo object is destroyed. It is up to the user to
    // do that.
    explicit MatlabInfo( mxArray *& ptr );

    // The destructor.
    ~MatlabInfo() { };

    void setfield( mxArray const *  ptr, char const * field ) ;
    void setfield( size_t n, Number const * x, char const * field ) ;
    void setfield( Number x, char const * field ) { setfield( 1, &x, field ) ; }

    mxArray *
    getfield_mx( char const * field ) const ;

    Number const *
    getfield( char const * field ) const ;

    // Access and modify the exit status and solution statistics.
    ApplicationReturnStatus getExitStatus() const;
    void setExitStatus( ApplicationReturnStatus status ) { setfield( (Number)status, "status") ; }
    void setFuncEvals(Index obj, Index con, Index grad, Index jac, Index hess);
    void setIterationCount( Index iter ) { setfield( iter, "iter") ; }
    void setCpuTime( Number cpu ) { setfield( cpu, "cpu") ; }
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
    Ipopt::IpoptApplication & app;  // The IPOPT application object.

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
    ~IpoptOptions() { };

    // The first function returns true if and only if the user has
    // specified a quasi-Newton approximation to the Hessian instead of
    // the exact Hessian. The second function returns true if and only
    // if the user has activated the derivative checker. The third
    // function returns true if and only if a user-specified scaling of
    // the problem is activated. The fourth function returns the print
    // level for the IPOPT console. The remaining two functions return
    // the floating-point value for positive and negative infinity,
    // respectively.
    bool useQuasiNewton () const ;
    bool useDerivChecker() const ;
    bool userScaling    () const ;
    Index printLevel() const ;
  };

  // Class Options.
  // -----------------------------------------------------------------
  // This class processes the options input from MATLAB.
  class Options {
  protected:
    Index               n ;       // The number of optimization variables.
    Index               m ;       // The number of constraints.
    std::vector<Number> lb ;      // Lower bounds on the variables.
    std::vector<Number> ub ;      // Upper bounds on the variables.
    std::vector<Number> cl ;      // Lower bounds on constraints.
    std::vector<Number> cu ;      // Upper bounds on constraints.
    std::vector<Number> zl ;      // Lagrange multipliers for lower bounds.
    std::vector<Number> zu ;      // Lagrange multipliers for upper bounds.
    std::vector<Number> lambda ;  // Lagrange multipliers for constraints.

    Number neginfty, posinfty ;

    IpoptOptions        ipopt ;   // The IPOPT options.

    // These are helper functions used by the class constructor.
    void loadLowerBounds( mxArray const * ptr ) ;
    void loadUpperBounds( mxArray const * ptr ) ;
    void loadConstraintBounds( mxArray const * ptr );
    void loadMultipliers( mxArray const * pt ) ;

  public:

    // The constructor expects as input a point to a MATLAB array, in
    // particular a structure array with the appropriate fields. Note
    // that the Options object does *not* possess an independent copy of
    // some of the MATLAB data (such as the auxiliary data).
    Options( Index n, Ipopt::IpoptApplication & app, mxArray const * ptr );

    // The destructor.
    ~Options();

    // Get the number of variables and the number of constraints.
    friend Index numvars        ( Options const & options ) { return options.n; }
    friend Index numconstraints ( Options const & options ) { return options.m; }

    // Access the lower and upper bounds on the variables and constraints.
    std::vector<Number> const & lowerbounds () const { return lb ; }
    std::vector<Number> const & upperbounds () const { return ub ; }
    std::vector<Number> const & constraintlb() const { return cl ; }
    std::vector<Number> const & constraintub() const { return cu ; }

    // Access the IPOPT options object.
    IpoptOptions const ipoptOptions() const { return ipopt; }

    // Access the Lagrange multpliers.
    std::vector<Number> const & multlb    () const { return zl ; }
    std::vector<Number> const & multub    () const { return zu ; }
    std::vector<Number> const & multconstr() const { return lambda ; }

  } ;

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
    CallbackFunctions const & funcs;    // Callback routines.
    Options const &           options;  // Further program info.
    MatlabInfo &              info;     // Info passed back to MATLAB.

  public:

    // The constructor.
    MatlabProgram( CallbackFunctions const & funcs,
		               Options           const & options,
                   MatlabInfo              & info );

    // The destructor.
    virtual ~MatlabProgram();

    // Method to return some info about the nonlinear program.
    virtual
    bool
    get_nlp_info ( Index & n,
                   Index & m,
                   Index & sizeOfJ,
                   Index & sizeOfH,
		               IndexStyleEnum & indexStyle ) ;

    // Return the bounds for the problem.
    virtual
    bool
    get_bounds_info( Index    n,
                     Number * lb,
                     Number * ub,
                     Index    m,
				             Number * cl,
                     Number * cu );

    // Return the starting point for the algorithm.
    virtual
    bool
    get_starting_point( Index    n,
                        bool     initializeVars,
                        Number * vars,
		                    bool     initializez,
                        Number * zl,
                        Number * zu,
		                    Index    m,
                        bool     initializeLambda,
                        Number * lambda ) ;

    // Compute the value of the objective.
    virtual
    bool
    eval_f( Index n, Number const * vars, bool ignore, Number & f ) ;

    // Compute the gradient of the objective.
    virtual
    bool
    eval_grad_f( Index n, Number const * vars, bool ignore, Number * grad );

    // Evaluate the constraint residuals.
    virtual
    bool
    eval_g( Index n, Number const * vars, bool ignore, Index m, Number * g );

    // This method either returns: 1.) The structure of the Jacobian
    // (if "Jacobian" is zero), or 2.) The values of the Jacobian (if
    // "Jacobian" is not zero).
    virtual
    bool
    eval_jac_g( Index          numVariables,
                Number const * variables,
			          bool           ignoreThis,
                Index          numConstraints,
			          Index          Jx_nnz,
                Index        * rows,
                Index        * cols,
                Number       * Jx ) ;

    // This method either returns: 1.) the structure of the Hessian of
    // the Lagrangian (if "Hessian" is zero), or 2.) the values of the
    // Hessian of the Lagrangian (if "Hesson" is not zero).
    virtual
    bool
    eval_h( Index          n,
            Number const * vars,
            bool           ignore,
            Number         sigma,
		        Index          m,
            Number const * lambda,
            bool           ignoretoo,
		        Index          Hx_nnz,
            Index        * rows,
            Index        * cols,
            Number       * Hx ) ;

    // This method is called when the algorithm is complete.
    virtual
    void
    finalize_solution( SolverReturn      status,
                       Index             numVariables,
                       Number const    * variables,
                       Number const    * zl,
                       Number const    * zu,
                       Index             numConstraints,
                       Number const    * constraints,
		                   Number const    * lambda,
                       Number            objective,
                       IpoptData const * ip_data,
                       IpoptCalculatedQuantities* ip_cq );

    // Intermediate callback method. It is called once per iteration
    // of the IPOPT algorithm.
    virtual
    bool
    intermediate_callback( AlgorithmMode mode,
                           Index         t,
                           Number        f,
				                   Number        inf_pr,
                           Number        inf_du,
				                   Number        mu,
                           Number        d_norm,
				                   Number        regularization_size,
				                   Number        alpha_du,
                           Number        alpha_pr,
				                   Index         ls_trials,
				                   IpoptData const * ip_data,
				                   IpoptCalculatedQuantities* ip_cq );

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
    MatlabJournal(EJournalLevel default_level) : Journal("matlab", default_level) {}

    // The destructor.
    virtual ~MatlabJournal() { };

  protected:

    // These functions override the functions in the Journal class.
    virtual
    void
    PrintImpl( EJournalCategory category,
               EJournalLevel    level,
               char const *     str) {
      mexPrintf(str);
    }

    virtual
    void
    PrintfImpl( EJournalCategory category,
                EJournalLevel    level,
                char const *     pformat,
                va_list ap ) ;

    virtual
    void FlushBufferImpl()
    {}

  };
}

#endif
