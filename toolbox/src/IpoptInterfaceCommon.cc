/*--------------------------------------------------------------------------*\
 |  file: IpoptInterfaceCommon.cc                                           |
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

#include "IpoptInterfaceCommon.hh"

#include "IpTNLP.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpOrigIpoptNLP.hpp"

using Ipopt::ApplicationReturnStatus;
using Ipopt::IpoptApplication;
using Ipopt::SmartPtr;
using Ipopt::IsValid;
using Ipopt::RegisteredOption;
using Ipopt::SmartPtr;
using Ipopt::Snprintf;

namespace IpoptInterface {

  /*
  //   ____
  //  / ___| _ __   __ _ _ __ ___  ___
  //  \___ \| '_ \ / _` | '__/ __|/ _ \
  //   ___) | |_) | (_| | |  \__ \  __/
  //  |____/| .__/ \__,_|_|  |___/\___|
  //        |_|
  */

  void
  SparseMatrix::setup( mxArray * ptr ) {
    // Get the height, width and number of non-zeros.
    numRows = Index( mxGetM(ptr) );
    numCols = Index( mxGetN(ptr) );
    // The precise number of non-zero elements is contained in the
    // last entry of the jc array. (There is one jc entry for each
    // column in the matrix, plus an extra one.)
    Jc.resize(numCols+1);
    mwIndex const * p = mxGetJc(ptr);
    for ( Index i=0; i <= numCols; ++i, ++p ) Jc[i] = Index(*p);
    //std::copy_n(mxGetJc(ptr),numCols+1,Jc.begin());
    // Copy the row and column indices, and the values of the nonzero entries.
    nnz = Jc.back();
    Ir.resize(nnz);
    p = mxGetIr(ptr);
    for ( Index i=0; i < nnz; ++i, ++p ) Ir[i] = Index(*p);
    // std::copy_n(mxGetIr(ptr),nnz,Ir.begin());
  }

  // extract the pattern of the sparse matrix in the ipopt format
  void
  SparseMatrix::getStructure( Index rows[], Index cols[] ) const {
    for ( Index c = 0, i = 0; c < numCols; ++c ) {
      for (; i < Jc[c+1]; ++i ) { cols[i] = c; rows[i] = Ir[i]; }
    }
  }

  // extyract values from the sparse matrix in ptr and store in values
  // which correspond to nonzeros of the sparse matrix described by
  // Ir and Jc. The nonzeros of sparse matrix in ptr must be a subset of the
  // nonzeros of the pattern described by Ir and Jc.
  void
  SparseMatrix::getValues( char const * func, mxArray * ptr, Number values[] ) const {
    // il patterm può essere un sottoinsieme
    mwIndex const * mxJc = mxGetJc(ptr);
    mwIndex const * mxIr = mxGetIr(ptr);
    double  const * v    = mxGetPr(ptr);

    std::fill_n( values, nnz, 0 );
    for ( Index c = 0; c < numCols; ++c ) {
      mwIndex i  = mxJc[c];
      mwIndex k  = Jc[c];
      mwIndex i1 = mxJc[c+1];
      mwIndex k1 = Jc[c+1];
      for (; i < i1; ++i, ++k ) {
        while ( Ir[k] < mxIr[i] && k < k1 ) ++k; // skip not set elements
        IPOPT_ASSERT(
          Ir[k] == mxIr[i],
          "In MATLAB function " << func <<
          "\nelement (" << mxIr[i]+1 << "," << c+1 <<
          ") not found in pattern"
        );
        values[k] = v[i];
      }
    }
  }

  static
  size_t
  from_matlab(
    mxArray const *       ptr,
    std::vector<Number> & data,
    char const            msg[]
  ) {

    size_t n = mxGetNumberOfElements(ptr);
    if ( mxIsCell(ptr) ) {

      // Compute the number of optimization variables.
      size_t nv = 0;
      for ( size_t i = 0; i < n; ++i ) {
        mxArray const * p = mxGetCell(ptr,i);  // Get the ith cell.
        IPOPT_ASSERT(
          p != nullptr,
          msg << " in from_matlab failed to access cell n. " << i
        );
        IPOPT_ASSERT(
          mxIsDouble(p) && !mxIsComplex(p) && !mxIsSparse(p),
          msg << " in from_matlab needs a cell array in which each cell is a"
                 " full real array in DOUBLE precision"
        );
        nv += mxGetNumberOfElements(p);
      }

      data.clear();
      data.reserve(nv);

      for ( size_t i = 0; i < n; ++i ) {
        mxArray const * p = mxGetCell(ptr,i);
        IPOPT_ASSERT(
          p != nullptr,
          msg << " in from_matlab failed to access cell n. " << i
        );
        nv = mxGetNumberOfElements(p);
        Number const * pp = mxGetPr(p);
        IPOPT_ASSERT(
          pp != nullptr,
          msg << " in from_matlab failed to access data for cell n. " << i
        );
        std::copy( pp, pp+nv, std::back_inserter(data) );
      }
    } else {
      IPOPT_ASSERT(
        mxIsDouble(ptr) && !mxIsComplex(ptr) && !mxIsSparse(ptr),
        msg << " in from_matlab data must be a real double non sparse vector, found: "
            << mxGetClassName(ptr)
      );
      Number const * p = mxGetPr(ptr);
      IPOPT_ASSERT(
        p != nullptr,
        msg << " in from_matlab failed to access data from vector"
      );

      data.clear();
      data.reserve(n);

      std::copy( p, p+n, std::back_inserter(data) );
    }
    return size_t(data.size());
  }

  /*
  //   ____
  //  / ___| _ __   __ _ _ __ ___  ___
  //  \___ \| '_ \ / _` | '__/ __|/ _ \
  //   ___) | |_) | (_| | |  \__ \  __/
  //  |____/| .__/ \__,_|_|  |___/\___|
  //        |_|
  */

  static
  bool
  isLowerTri( mxArray const * ptr ) {

    // Get the height and width of the matrix.
    mwIndex nrow = mxGetM(ptr);
    mwIndex ncol = mxGetN(ptr);

    // Check whether the sparse matrix is symmetric.
    if ( nrow != ncol ) return false;

    // A sparse lower triangular matrix has the property that the
    // column indices are always less than or equal to the row indices.
    bool      ok = true;  // The return value.
    mwIndex * jc = mxGetJc(ptr);
    mwIndex * ir = mxGetIr(ptr);
    mwIndex i = 0;
    for ( mwIndex j = 0; j < ncol && ok; ++j )
      for (; i < jc[j+1] && ok; ++i )
        ok = (j <= ir[i]);
    return ok;
  }

  static
  bool
  inIncOrder( mxArray const * ptr ) {
    Index    ncol = (Index) mxGetN(ptr);
    mwIndex* jc   = mxGetJc(ptr);
    mwIndex* ir   = mxGetIr(ptr);
    bool     ok   = true;
    for ( Index j = 0; j < ncol && ok; ++j )
      for ( mwIndex i = jc[j]+1; i < jc[j+1] && ok; ++i )
        ok = ( ir[i-1] < ir[i] );
    return ok;
  }

  static
  void
  mxToIpoptInf( std::vector<Number> & v, Number ipopt_inf ) {
    for ( std::vector<Number>::iterator iv = v.begin(); iv != v.end(); ++iv )
      if ( mxIsInf(*iv) )
        *iv = ipopt_inf;
  }

  // This constructor accepts as input a pointer to a MATLAB array. It
  // is up to the user to ensure that the MATLAB array is a valid
  // function handle.
  void
  MatlabFunctionHandle::bind( mxArray const * p, char const * error_msg ) {
    IPOPT_ASSERT(
      p != nullptr,
      "You must specify a callback routine for computing " << error_msg
    );
    IPOPT_ASSERT(
      !mxIsEmpty(p) && mxIsClass(p,"function_handle"),
      "You did not provide a valid function handle for computing " << error_msg
    );

    if ( f != nullptr ) mxDestroyArray(f);
    f = mxDuplicateArray(p);

    mxArray *outputs[1];
    mxArray *exception = mexCallMATLABWithTrap( 1, outputs, 1, &f, "func2str" );

    if ( exception != nullptr  ) {
      mxArray *msg;
      mexCallMATLAB(1, &msg, 1, &exception, "getReport" );
      char * c_msg = mxArrayToString(msg);
      std::string cpp_msg = c_msg;
      mxFree(c_msg);
      IPOPT_DO_ERROR("in function: func2str\n" << cpp_msg);
    }

    char * c_msg = mxArrayToString(outputs[0]);
    f_name = c_msg;
    mxFree(c_msg);
  }

  void
  MatlabFunctionHandle::eval(
    Index           n_lhs,
    mxArray       * lhs[],
    Index           n_rhs,
    mxArray const * rhs[]
  ) const {

    // Construct the inputs to "feval".
    mxArray** finputs = new mxArray*[n_rhs+1];
    finputs[0] = f;
    for ( Index i = 0; i < n_rhs; ++i )
      finputs[i+1] = mxDuplicateArray(rhs[i]);

    // Call "feval".
    mxArray *exception = mexCallMATLABWithTrap(
      n_lhs, lhs, n_rhs+1, finputs, "feval"
    );

    // Free the dynamically allocated memory.
    for ( Index i = 1; i <= n_rhs; ++i )
      mxDestroyArray(finputs[i]);
    delete[] finputs;

    if ( exception != nullptr  ) {
      mxArray *msg;
      mexCallMATLAB(1, &msg, 1, &exception, "getReport" );
      char * c_msg = mxArrayToString(msg);
      std::string cpp_msg = c_msg;
      mxFree(c_msg);
      IPOPT_DO_ERROR("in function: " << f_name << "\n" << cpp_msg);
    }
  }

  /*
  //   ___        __
  //  |_ _|_ __  / _| ___
  //   | || '_ \| |_ / _ \
  //   | || | | |  _| (_) |
  //  |___|_| |_|_|  \___/
  */
  // Function definitions for class MatlabInfo.
  // ------------------------------------------------------------------
  MatlabInfo::MatlabInfo( mxArray *& ptr) : ptr(nullptr) {
    // Create the structure array.

    char const * fieldnames[9] = {
      "x",
      "status",
      "zl",
      "zu",
      "lambda",
      "iter",
      "cpu",
      "objective",
      "eval"
    };
    this->ptr = ptr = mxCreateStructMatrix(1,1,9,fieldnames);

    // Initialize fields.
    for ( Index i = 0; i < 8; ++i ) // do not initialize "eval"
      mxSetField(ptr,0,fieldnames[i],mxCreateDoubleScalar(0));

    // Build Eval Structure
    static char const * evalfields[5] = {
      "objective",
      "constraints",
      "gradient",
      "jacobian",
      "hessian"
    };

    mxArray *evalStruct = mxCreateStructMatrix(1,1,5,evalfields);
    for ( Index i = 0; i < 5; ++i )
      mxSetField( evalStruct, 0, evalfields[i], mxCreateDoubleScalar(0) );

    mxSetField( ptr, 0, "eval", evalStruct );
  }

  void
  MatlabInfo::setfield( mxArray const * pin, char const * field ) {
    mxArray * p = mxGetField( ptr, 0, field );
    IPOPT_ASSERT(
      p != nullptr, "MatlabInfo::setfield missing field `" << field << "'"
    );
    // First destroy any previous values.
    mxDestroyArray( p );
    mxSetField( ptr, 0, field, mxDuplicateArray( pin  ));
  }

  void
  MatlabInfo::setfield( size_t n, Number const * x, char const * field ) {
    mxArray * p = mxGetField(ptr,0,field);
    IPOPT_ASSERT(
      p != nullptr, "MatlabInfo::setfield missing field `" << field << "'"
    );
    // First destroy any previous values.
    mxDestroyArray(p);
    p = mxCreateDoubleMatrix( n, 1, mxREAL );
    std::copy_n( x, n, mxGetPr(p) );
    mxSetField( ptr, 0, field, p );
  }

  mxArray *
  MatlabInfo::getfield_mx( char const * field ) const {
    mxArray * p = mxGetField( ptr, 0, field );
    IPOPT_ASSERT(
      p != nullptr,
      "MatlabInfo::getfield_mx missing field `" << field << "'"
    );
    return p;
  }

  Number const *
  MatlabInfo::getfield( char const * field ) const {
    Number const * p = mxGetPr(getfield_mx(field));
    IPOPT_ASSERT(
      p != nullptr,
      "MatlabInfo::getfield('" << field << "') cant access values"
    );
    return p;
  }

  ApplicationReturnStatus
  MatlabInfo::getExitStatus() const {
    mxArray const * p = getfield_mx("status");
    return (ApplicationReturnStatus) (int) *mxGetPr(p);
  }

  void
  MatlabInfo::setFuncEvals(
    Index obj,
    Index con,
    Index grad,
    Index jac,
    Index hess
  ) {
    mxArray* p = getfield_mx("eval");
    *mxGetPr( mxGetField( p, 0, "objective") )   = obj;
    *mxGetPr( mxGetField( p, 0, "constraints") ) = con;
    *mxGetPr( mxGetField( p, 0, "gradient") )    = grad;
    *mxGetPr( mxGetField( p, 0, "jacobian") )    = jac;
    *mxGetPr( mxGetField( p, 0, "hessian") )     = hess;
  }

  /*
  //    ___        _   _
  //   / _ \ _ __ | |_(_) ___  _ __  ___
  //  | | | | '_ \| __| |/ _ \| '_ \/ __|
  //  | |_| | |_) | |_| | (_) | | | \__ \
  //   \___/| .__/ \__|_|\___/|_| |_|___/
  //        |_|
  */

  // Function definitions for class IpoptOptions.
  // -----------------------------------------------------------------
  IpoptOptions::IpoptOptions( IpoptApplication & app, mxArray const * ptr ) : app(app) {

    if ( ptr != nullptr ) { // avoid segfault if no ipopt field in options
      // Check to make sure the MATLAB array is a structure array.
      IPOPT_ASSERT(
        mxIsStruct(ptr),
        "The OPTIONS input must be a structure array; found " << mxGetClassName(ptr) <<
        "\ntype HELP STRUCT in the MATLAB console for more information"
      );

      // Each field in the structure array should correspond to an option
      // in IPOPT. Repeat for each field.
      Index n = mxGetNumberOfFields(ptr);
      for ( Index i = 0; i < n; ++i ) {
        char const * label = mxGetFieldNameByNumber(ptr,i);
        mxArray *    p     = mxGetFieldByNumber(ptr,0,i);
        setOption(label,p);
      }
    }
  }

  bool
  IpoptOptions::useQuasiNewton() const {
    std::string value;
    app.Options()->GetStringValue("hessian_approximation",value,"");
    bool b = !value.compare("limited-memory");
    return b;
  }

  bool
  IpoptOptions::useDerivChecker() const {
    std::string value;
    app.Options()->GetStringValue("derivative_test",value,"");
    bool b = value.compare("none");
    return b;
  }

  bool
  IpoptOptions::userScaling() const {
    std::string value;
    app.Options()->GetStringValue("nlp_scaling_method",value,"");
    bool b = !value.compare("user-scaling");
    return b;
  }

  Index
  IpoptOptions::printLevel() const {
    Index value;  // The return value.
    app.Options()->GetIntegerValue("print_level",value,"");
    return value;
  }

  void
  IpoptOptions::setOption( char const * label, mxArray const * ptr ) {
    // Check to make sure we have a valid option.
    SmartPtr<const RegisteredOption> option = app.RegOptions()->GetOption(label);
    IPOPT_ASSERT(
      IsValid(option),
      "You have specified a nonexistent IPOPT option (\"" << label << "\")"
    );
    Ipopt::RegisteredOptionType type = option->Type();
    if      (type == Ipopt::OT_String)  setStringOption(label,ptr);
    else if (type == Ipopt::OT_Integer) setIntegerOption(label,ptr);
    else                                setNumberOption(label,ptr);
  }

  void
  IpoptOptions::setStringOption( char const * label, mxArray const * ptr ) {
    // Check whether the option value is a string.
    IPOPT_ASSERT(
      mxIsChar(ptr),
      "IPOPT option value for option \"" << label << "\" should be a string"
    );

    // Get the option value.
    char* value = mxArrayToString(ptr);

    // Set the option.
    bool success = app.Options()->SetStringValue(label,value);
    IPOPT_ASSERT(
      success, "Invalid value for IPOPT option \"" << label << "\""
    );

    // Free the dynamically allocated memory.
    mxFree(value);
  }

  void
  IpoptOptions::setIntegerOption( char const * label, mxArray const * ptr) {

    // Check whether the option value is a number.
    IPOPT_ASSERT(
      mxIsDouble(ptr),
      "IPOPT option value for option \"" << label << "\" should be an integer"
    );

    // Set either the integer option.
    Number value   = mxGetScalar(ptr);
    bool   success = app.Options()->SetIntegerValue(label,(Index) value);
    IPOPT_ASSERT(
      success, "Invalid value for integer IPOPT option \"" << label << "\""
    );
  }

  void
  IpoptOptions::setNumberOption( char const * label, mxArray const * ptr ) {

    // Check whether the option value is a number.
    IPOPT_ASSERT(
      mxIsDouble(ptr),
      "IPOPT option value for option \"" << label << "\" should be a number"
    );

    // Set either the numeric option.
    Number value   = mxGetScalar(ptr);
    bool   success = app.Options()->SetNumericValue(label,value);
    IPOPT_ASSERT(
      success, "Invalid value for numeric IPOPT option \"" << label << "\""
    );
  }

  // Function definitions for class Options.
  // -----------------------------------------------------------------
  Options::Options( Index n_in, Ipopt::IpoptApplication& app, mxArray const * ptr)
  : n(n_in)
  , m(0) // will be set by loadConstraintBounds
  , ipopt(app,mxGetField(ptr,0,"ipopt")) // Process the IPOPT options.
  {
    app.Options()->GetNumericValue("nlp_upper_bound_inf",posinfty,"");
    app.Options()->GetNumericValue("nlp_lower_bound_inf",neginfty,"");

    // Load the bounds on the variables.
    loadLowerBounds(ptr);
    loadUpperBounds(ptr);

    // Load the bounds on the constraints.
    loadConstraintBounds(ptr);

    // Load the Lagrange multipliers.
    loadMultipliers(ptr);
  }

  Options::~Options() {
  }

  // Function definitions for static members of class Options.
  // -----------------------------------------------------------------
  void
  Options::loadLowerBounds( mxArray const * ptr ) {
    mxArray const * p = mxGetField(ptr,0,"lb");
    if ( p != nullptr ) {
      // Load the upper bounds and check to make sure they are valid.
      Index N = Index( from_matlab(p,lb,"Options::loadLowerBounds") );
      IPOPT_ASSERT(
        N == n,
        "Lower bounds array must have one element for each optimization"
        " variable, found a vector of N = " << N <<
        " elements, expected n = " << n
      );
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      mxToIpoptInf( lb, neginfty );
    } else {
      // If the lower bounds have not been specified, set them to negative infinity.
      lb.resize(n);
      std::fill( lb.begin(), lb.end(), neginfty );
    }
  }

  void
  Options::loadUpperBounds( mxArray const * ptr ) {
    mxArray const * p = mxGetField(ptr,0,"ub");
    if ( p != nullptr ) {
      // Load the upper bounds and check to make sure they are valid.
      Index N = Index( from_matlab( p, ub, "Options::loadUpperBounds" ) );
      IPOPT_ASSERT(
        N == n,
        "Upper bounds array must have one element for each optimization variable"
      );
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      mxToIpoptInf( ub, posinfty );
    } else {
      // If the upper bounds have not been specified, set them to positive infinity.
      ub.resize(n);
      std::fill( ub.begin(), ub.end(), posinfty );
    }
  }

  void
  Options::loadConstraintBounds( mxArray const * ptr ) {
    m = 0;  // The return value is the number of constraints.
    // LOAD CONSTRAINT BOUNDS.
    // If the user has specified constraints bounds, then she must
    // specify *both* the lower and upper bounds.
    mxArray const * pl = mxGetField(ptr,0,"cl");
    mxArray const * pu = mxGetField(ptr,0,"cu");
    if ( pl != nullptr || pu != nullptr ) {
      // Check to make sure the constraint bounds are valid.
      IPOPT_ASSERT(
        pl != nullptr && pu != nullptr,
        "You must specify both lower and upper bounds on the constraints"
      );
      IPOPT_ASSERT(
        mxIsDouble(pl) && mxIsDouble(pu) &&
        (mxGetNumberOfElements(pl) == mxGetNumberOfElements(pu)),
        "The lower and upper bounds must both be double-precision arrays "
        "with the same number of elements"
      );

      // Get the number of constraints.
      m = (Index) mxGetNumberOfElements(pl);

      // Load the lower bounds on the constraints and convert MATLAB's
      // convention of infinity to IPOPT's convention of infinity.
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      cl.resize(m); std::copy_n(mxGetPr(pl),m,cl.begin()); mxToIpoptInf( cl, neginfty );
      cu.resize(m); std::copy_n(mxGetPr(pu),m,cu.begin()); mxToIpoptInf( cu, posinfty );
    }
  }

  void
  Options::loadMultipliers( mxArray const * ptr ) {
    // Load the Lagrange multipliers associated with the lower bounds.
    mxArray const * p = mxGetField(ptr,0,"zl");
    static char const * msg1 = "must be a double-precision array with one element\nfor each optimization variable";
    if ( p != nullptr ) {
      static char const * msg = "The initial point for the Lagrange multipliers associated with the lower bounds\n";
      IPOPT_ASSERT( mxIsDouble(p) && !mxIsComplex(p), msg << msg1 );
      IPOPT_ASSERT(
        mxGetNumberOfElements(p) == n,
        msg << "must be a vector of length " << n <<
        " found of lenght " << mxGetNumberOfElements(p)
      );
      zl.resize(n);
      std::copy_n(mxGetPr(p),n,zl.begin());
    } else {
      zl.clear();
    }
    // Load the Lagrange multipliers associated with the upper bounds.
    p = mxGetField(ptr,0,"zu");
    if ( p != nullptr ) {
      static char const * msg = "The initial point for the Lagrange multipliers associated with the upper bounds\n";
      IPOPT_ASSERT( mxIsDouble(p) && !mxIsComplex(p), msg << msg1 );
      IPOPT_ASSERT(
        mxGetNumberOfElements(p) == n,
        msg << "must be a vector of length " << n <<
        " found of lenght " << mxGetNumberOfElements(p)
      );
      zu.resize(n);
      std::copy_n(mxGetPr(p),n,zu.begin());
    } else {
      zu.clear();
    }
    // Load the Lagrange multipliers associated with the equality and inequality constraints.
    p = mxGetField(ptr,0,"lambda");
    if ( p != nullptr ) {
      static char const * msg = "The initial point for the Lagrange multipliers associated with the constraints\n";
      IPOPT_ASSERT(
        mxIsDouble(p) && !mxIsComplex(p),
        msg << "must be a double-precision array with one element for each constraint"
      );
      IPOPT_ASSERT(
        mxGetNumberOfElements(p) == m,
        msg << "must be a vector of length " << n <<
        " found of length " << mxGetNumberOfElements(p)
      );
      lambda.resize(m);
      std::copy_n(mxGetPr(p),m,lambda.begin());
    } else {
      lambda.clear();
    }
  }

  /*
  //    ____      _ _ _                _
  //   / ___|__ _| | | |__   __ _  ___| | __
  //  | |   / _` | | | '_ \ / _` |/ __| |/ /
  //  | |__| (_| | | | |_) | (_| | (__|   <
  //   \____\__,_|_|_|_.__/ \__,_|\___|_|\_\
  */
  // Functions for class CallbackFunctions.
  // -----------------------------------------------------------------
  CallbackFunctions::CallbackFunctions(
    mxArray const * mx_x0,
    mxArray const * ptr
  ) {

    from_matlab( mx_x0, x0, "in CallbackFunctions" );

    mxArray const * p;  // A pointer to a MATLAB array.
    this -> x_is_cell_array = mxIsCell( mx_x0 );
    if ( this -> x_is_cell_array ) {
      this -> mx_x_nv = 0;
      this -> mx_x_nc = Index( mxGetNumberOfElements( mx_x0 ));
      this -> mx_x    = mxCreateCellMatrix( this -> mx_x_nc, 1 );
      for ( Index i = 0; i < this -> mx_x_nc; ++i ) {
        mxArray const * p = mxGetCell( mx_x0, i );
        this -> mx_x_nv += (Index) mxGetNumberOfElements(p);
        mxSetCell( this -> mx_x, i, mxDuplicateArray(p) );
      }
    } else {
      this -> mx_x    = mxDuplicateArray( mx_x0 );
      this -> mx_x_nv = Index( mxGetNumberOfElements( mx_x ) );
      this -> mx_x_nc = 0;
    }

    // Check whether we are provided with a structure array.
    IPOPT_ASSERT(
      mxIsStruct(ptr),
      "in CallbackFunctions the second input must be a STRUCT, found: " <<
      mxGetClassName(ptr)
    );

    // Get the function handle for computing the objective.
    p = mxGetField( ptr, 0, "objective" );
    obj_func.bind( p, "the objective function" );

    // Get the function handle for computing the gradient.
    p = mxGetField( ptr, 0, "gradient" );
    grad_func.bind( p, "the gradient of the objective" );

    // Get the function handle for computing the constraints, if such a function was specified.
    p = mxGetField( ptr, 0, "constraints" );
    if ( p != nullptr )
      constraint_func.bind( p, "the response of the constraints" );

    // Get the function handle for computing the Jacobian.
    // This function is necessary if there are constraints.
    p = mxGetField( ptr, 0, "jacobian" );
    if ( p != nullptr ) {
      jacobian_func.bind( p, "the first derivatives (Jacobian) of the constraints" );
    } else {
      IPOPT_ASSERT(
        !constraint_func,
        "You must provide a function that returns the first derivatives"
        " (Jacobian) of the constraints"
      );
    }

    // Get the function handle for computing the sparsity structure of the Jacobian.
    // This function is necessary if the Jacobian is being computed.
    p = mxGetField( ptr, 0, "jacobianstructure" );
    if ( p != nullptr ) {
      jacstruc_func.bind(p,"the sparsity structure of the Jacobian");
    } else {
      IPOPT_ASSERT(
        !jacobian_func,
        "You must provide a function that returns the sparsity structure of the Jacobian"
      );
    }

    // Get the function handle for computing the Hessian. This function is always optional.
    p = mxGetField( ptr, 0, "hessian" );
    if ( p != nullptr ) hessian_func.bind(p,"the Hessian of the Lagrangian");

    // Get the function handle for computing the sparsity structure of the Hessian of the Lagrangian.
    // This function is necessary if the Hessian is being computed.
    p = mxGetField( ptr, 0, "hessianstructure" );
    if ( p != nullptr ) {
      hesstruc_func.bind( p,"the sparsity structure of the Hessian" );
    } else {
      IPOPT_ASSERT(
        !hessian_func,
        "You must provide a function that returns the sparsity structure of the Hessian"
      );
    }

    // Get the iterative callback function handle. This function is always optional.
    p = mxGetField( ptr, 0, "iterfunc" );
    if ( p != nullptr )
      iter_func.bind( p, "the iterative callback" );
  }

  CallbackFunctions::~CallbackFunctions() {
  }

  bool
  CallbackFunctions::from_cell_array(
    mxArray const * ptr,
    Index           n,
    Number        * x
  ) const {
    if ( !x_is_cell_array ) return false;
    Index ntot = 0;
    for ( Index i = 0; i < this -> mx_x_nc; ++i ) {
      mxArray const * p = mxGetCell(ptr,i);
      Index nn = (Index) mxGetNumberOfElements(p);
      ntot += nn;
      if ( ntot > n ) return false;
      std::copy_n(mxGetPr(p),nn,x);
      x += nn;
    }
    return ntot == n;
  }

  void
  CallbackFunctions::fillx( Index n, Number const * x ) const {
    if ( x_is_cell_array ) {
      for ( Index i = 0; i < this -> mx_x_nc; ++i ) {
        mxArray const * p = mxGetCell(mx_x,i);
        Index nn = (Index) mxGetNumberOfElements(p);
        std::copy_n(x,nn,mxGetPr(p));
        x += nn;
      }
    } else {
      std::copy_n( x, n, mxGetPr(mx_x) );
    }
  }

  Number
  CallbackFunctions::computeObjective( Index n, Number const x[] ) const {

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    obj_func.eval(1,&ptr,1,(mxArray const **)&mx_x);

    // Get the output from the MATLAB callback function, which is the
    // value of the objective function at x.
    IPOPT_ASSERT(
      mxIsDouble(ptr) && !mxIsComplex(ptr),
      "In MATLAB function " << obj_func.name() << "\n"
      "The first return value of the objective callback function"
      " must be a real double scalar, found: " << mxGetClassName(ptr)
    );
    IPOPT_ASSERT(
      mxGetNumberOfElements(ptr) == 1,
      "In MATLAB function " << obj_func.name() << "\n"
      "The first return value of the objective callback function has " <<
      mxGetNumberOfElements(ptr) << " elements, expected 1"
    );

    if ( mxIsSparse(ptr) ) {
      // convert sparse objective (unlikely but possible) to full
      mxArray * ptr1;
      mexCallMATLAB(1, &ptr1, 1, &ptr, "full");
      std::swap(ptr,ptr1);
      mxDestroyArray(ptr1); // destroy old sparse vector
    }

    Number f = *mxGetPr(ptr);

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
    return f;
  }

  void
  CallbackFunctions::computeGradient(
    Index        n,
    Number const x[],
    Number       g[]
  ) const {

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    grad_func.eval(1,&ptr,1,(mxArray const **)&mx_x);

    // se cell array converto
    if ( mxIsCell(ptr) ) {
      IPOPT_ASSERT(
        from_cell_array( ptr, n, g ),
        "In MATLAB function " << grad_func.name() << "\n"
        "The gradient callback return a cell array not with the same structure of x0"
      );
    } else {
      // Get the output from the MATLAB callback function, which is the
      // value of the gradient of the objective function at x.
      IPOPT_ASSERT(
        mxIsDouble(ptr) && !mxIsComplex(ptr),
        "In MATLAB function " << grad_func.name() << "\n"
        "The gradient callback must return a real double vector, found: " <<
        mxGetClassName(ptr)
      );
      if (mxIsSparse(ptr)) {
        mxArray * ptr1;
        // convert sparse gradient to full (simplest method, not fastest)
        mexCallMATLAB(1, &ptr1, 1, &ptr, "full");
        std::swap(ptr,ptr1);
        mxDestroyArray(ptr1); // destroy old sparse vector
      }
      IPOPT_ASSERT(
        mxGetNumberOfElements(ptr) == n,
        "In MATLAB function " << grad_func.name() << "\n"
        "The gradient callback must return a real double vector of size = " << n <<
        " while it return a vector of size = " << mxGetNumberOfElements(ptr)
      );

      std::copy_n( mxGetPr(ptr), n, g );
    }

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
  }

  void
  CallbackFunctions::computeConstraints(
    Index        n,
    Number const x[],
    Index        m,
    Number       c[]
  ) const {

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    constraint_func.eval(1,&ptr,1,(mxArray const **)&mx_x);

    // Get the output from the MATLAB callback function, which is the
    // value of vector-valued constraint function at x.
    IPOPT_ASSERT(
      mxIsDouble(ptr) && !mxIsComplex(ptr),
      "In MATLAB function " << constraint_func.name() << "\n"
      "The constraint callback must return a real double vector, found: " <<
      mxGetClassName(ptr)
    );

    IPOPT_ASSERT(
      m == mxGetNumberOfElements(ptr),
      "In MATLAB function " << constraint_func.name() << "\n"
      "The contraints callback must return a real double vector of size = " << m <<
      " while it return a vector of size = " << mxGetNumberOfElements(ptr)
    );

    if (mxIsSparse(ptr)) {
      mxArray * ptr1;
      // convert sparse constraint vector (unlikely but possible) to full
      mexCallMATLAB(1, &ptr1, 1, &ptr, "full");
      std::swap(ptr,ptr1);
      mxDestroyArray(ptr1); // destroy old sparse vector
    }
    std::copy_n(mxGetPr(ptr),m,c);

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
  }

  static
  void
  checkJaconbian( char const * name, Index n, Index m, mxArray * ptr ) {
    // Get the output from the MATLAB callback function, which is the
    // sparse matrix specifying the value the Jacobian.
    IPOPT_ASSERT(
      mxIsSparse(ptr),
      "In MATLAB function " << name << "\n"
      "Jacobian must be a real sparse matrix, found full matrix " <<
      mxGetM(ptr) << " x " << mxGetN(ptr)
    );
    IPOPT_ASSERT(
      !mxIsComplex(ptr),
      "In MATLAB function " << name << "\n"
      "Jacobian must be a real sparse matrix, found complex sparse matrix"
    );
    IPOPT_ASSERT(
      mxGetM(ptr) == m && mxGetN(ptr) == n,
      "In MATLAB function " << name << "\n"
      "Jacobian must be an (m=" << m << ") x (n=" << n << ") sparse matrix\n"
      "where m is the number of constraints and n is the number of variables,\n"
      "but found m=" << mxGetM(ptr) << " and n=" << mxGetN(ptr)
    );
    IPOPT_ASSERT(
      inIncOrder(ptr),
      "In MATLAB function " << name << "\n"
      "Jacobian must be a sparse matrix with row indices in increasing order"
    );
  }

  void
  CallbackFunctions::loadJacobianStructure( Index n, Index m ) const {

    // Call the MATLAB callback function.
    mxArray * ptr;
    jacstruc_func.eval(1,&ptr,0,(mxArray const **)nullptr);

    checkJaconbian( jacstruc_func.name().c_str(), n, m, ptr );
    Jacobian.setup( ptr );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
  }

  static
  void
  checkHessian( char const * name, Index n, mxArray * ptr ) {
    // Get the output from the MATLAB callback function, which is the
    // sparse matrix specifying the structure of the Hessian.
    IPOPT_ASSERT(
      !mxIsComplex(ptr),
      "In MATLAB function " << name << "\n"
      "Hessian must be a real matrix, found: " << mxGetClassName(ptr)
    );
    IPOPT_ASSERT(
      mxIsSparse(ptr),
      "In MATLAB function " << name << " Hessian must be a sparse matrix"
    );
    IPOPT_ASSERT(
      mxGetM(ptr) == n && mxGetN(ptr) == n,
      "In MATLAB function " << name << "\n"
      "Hessian must be a " << n << " x " << n << " matrix\n"
      "found an " << mxGetM(ptr)  << " x " <<  mxGetN(ptr) << " matrix"
    );
    IPOPT_ASSERT(
      isLowerTri(ptr),
      "In MATLAB function " << name << "\n"
      "Hessian must be a real sparse lower triangular matrix.\n"
      "Matrix is not lower triangular"
    );
    IPOPT_ASSERT(
      inIncOrder(ptr),
      "In MATLAB function " << name << "\n"
      "Hessian must be an n x n sparse, symmetric and lower triangular matrix "
      "with row indices in increasing order.\n"
      "Row indices are not in increasing order"
    );
  }

  void
  CallbackFunctions::loadHessianStructure( Index n ) const {

    // Call the MATLAB callback function.
    mxArray * ptr;
    hesstruc_func.eval(1,&ptr,0,(mxArray const **)nullptr);

    checkHessian( hesstruc_func.name().c_str(), n, ptr );
    Hessian.setup( ptr );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
  }

  void
  CallbackFunctions::computeJacobian(
    Index        m,
    Index        n,
    Number const x[],
    Number       values[]
  ) const {

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    jacobian_func.eval(1,&ptr,1,(mxArray const **)&mx_x);

    checkJaconbian( jacobian_func.name().c_str(), n, m, ptr );

    if (!mxIsDouble(ptr)) {
      mxArray * ptr1;
      // convert non-double Jacobian to double
      mexCallMATLAB(1, &ptr1, 1, &ptr, "double");
      std::swap(ptr,ptr1);
      mxDestroyArray(ptr1);
    }

    Jacobian.getValues( jacobian_func.name().c_str(), ptr, values );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
  }

  void
  CallbackFunctions::computeHessian(
    Index        n,
    Number const x[],
    Number       sigma,
    Index        m,
    Number const lambda[],
    Number       values[]
  ) const {
 
    fillx( n, x );

    // Create the input arguments to the MATLAB routine, sigma and lambda.
    mxArray * psigma  = mxCreateDoubleScalar(sigma);
    mxArray * plambda = mxCreateDoubleMatrix(m,1,mxREAL);
    std::copy_n( lambda, m, mxGetPr(plambda) );

    // Call the MATLAB callback function.
    mxArray const * inputs[3] = { mx_x, psigma, plambda };
    mxArray       * ptr;
    hessian_func.eval( 1, &ptr, 3, inputs );

    checkHessian( hesstruc_func.name().c_str(), n, ptr );
    /*
    if (!mxIsDouble(ptr)) {
      mxArray * ptr1;
      // convert non-double Hessian to double
      mexCallMATLAB( 1, &ptr1, 1, &ptr, "double" );
      std::swap( ptr, ptr1 );
      mxDestroyArray( ptr1 );
    }
    */

    IPOPT_ASSERT(
      mxIsDouble(ptr),
      "in computeHessian output data must be of type double, found: " <<
      mxGetClassName(ptr)
    );
    IPOPT_ASSERT(
      !mxIsComplex(ptr),
      "in computeHessian output data cannot be complex!"
    );
    IPOPT_ASSERT(
      mxIsSparse(ptr),
      "in computeHessian output data must be a sparse matrix"
    );
    Number const * v = mxGetPr(ptr);
    IPOPT_ASSERT(
      v != nullptr, "in computeHessian failed to access data from vector"
    );

    Hessian.getValues( jacobian_func.name().c_str(), ptr, values );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
    mxDestroyArray(psigma);
    mxDestroyArray(plambda);
  }

  bool
  CallbackFunctions::iterCallback(
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
    Ipopt::IpoptData const *           ip_data,
    Ipopt::IpoptCalculatedQuantities * ip_cq
  ) const {

    if ( !iter_func ) return true;

    // Create the input arguments to the MATLAB routine.
    mxArray* pt = mxCreateDoubleScalar(iter);
    mxArray* pf = mxCreateDoubleScalar(f);

    // Create structure to hold extra IPOPT variables
    static char const * varfields[8] = {
      "inf_pr", "inf_du", "mu",
      "d_norm", "regularization_size",
      "alpha_du", "alpha_pr", "ls_trials"
    };

    mxArray *varStruct = mxCreateStructMatrix(1,1,8,varfields);
    mxSetField(varStruct,0,"inf_pr",mxCreateDoubleScalar(inf_pr));
    mxSetField(varStruct,0,"inf_du",mxCreateDoubleScalar(inf_du));
    mxSetField(varStruct,0,"mu",mxCreateDoubleScalar(mu));
    mxSetField(varStruct,0,"d_norm",mxCreateDoubleScalar(d_norm));
    mxSetField(varStruct,0,"regularization_size",mxCreateDoubleScalar(regularization_size));
    mxSetField(varStruct,0,"alpha_du",mxCreateDoubleScalar(alpha_du));
    mxSetField(varStruct,0,"alpha_pr",mxCreateDoubleScalar(alpha_pr));
    mxSetField(varStruct,0,"ls_trials",mxCreateDoubleScalar(ls_trials));

    // Call the MATLAB callback function.
    mxArray const * inputs[3] = { pt, pf, varStruct };
    mxArray       * outputs[1];

    iter_func.eval(1,outputs,3,inputs);

    // Get the output from the MATLAB callback function, which is a
    // boolean value telling whether or not IPOPT should continue.
    IPOPT_ASSERT(
      mxIsLogicalScalar(outputs[0]),
      "The return value for the iterative callback must either be TRUE or FALSE"
    );
    bool b = mxIsLogicalScalarTrue(outputs[0]);

    // Free the dynamically allocated memory.
    mxDestroyArray(outputs[0]);
    mxDestroyArray(pt);
    mxDestroyArray(pf);
    mxDestroyArray(varStruct);

    return b;
  }

  /*
  //   ____
  //  |  _ \ _ __ ___   __ _ _ __ __ _ _ __ ___
  //  | |_) | '__/ _ \ / _` | '__/ _` | '_ ` _ \
  //  |  __/| | | (_) | (_| | | | (_| | | | | | |
  //  |_|   |_|  \___/ \__, |_|  \__,_|_| |_| |_|
  //                   |___/
  */

  // Function definitions for class MatlabProgram.
  // ---------------------------------------------------------------
  MatlabProgram::MatlabProgram(
    CallbackFunctions const & funcs,
    Options           const & options,
    MatlabInfo              & info
  )
  : funcs(funcs)
  , options(options)
  , info(info) { }

  MatlabProgram::~MatlabProgram() {
  }

  bool
  MatlabProgram::get_nlp_info(
    Index          & n,
    Index          & m,
    Index          & J_nnz,
    Index          & H_nnz,
    IndexStyleEnum & indexStyle 
  ) {
    // Get the number of variables and constraints.
    n = numvars(options);
    m = numconstraints(options);

    // Get the size of the Jacobian.
    if ( m > 0 ) {
      IPOPT_ASSERT(
        funcs.jacobianFuncIsAvailable(),
        "You need to specify the callback functions for computing the Jacobian "
        "and the sparsity structure of the Jacobian"
      );
      funcs.loadJacobianStructure(n,m);
      J_nnz = funcs.getJacobianNnz();
    } else {
      J_nnz = 0;
    }
    // Get the size of the symmetric Hessian matrix. We don't need to
    // store the actual result, we just need to look at the number of
    // non-zero entries in the lower triangular part of the matrix.
    if (!options.ipoptOptions().useQuasiNewton()) {
      IPOPT_ASSERT(
        funcs.hessianFuncIsAvailable(),
        "You need to specify the callback functions for computing the Hessian "
        "and the sparsity structure of the Hessian"
      );
      funcs.loadHessianStructure(n);
      H_nnz = funcs.getHessianNnz();
    } else {
      H_nnz = 0;
    }
    // Use C-style indexing.
    indexStyle = C_STYLE;
    return true;
  }

  bool
  MatlabProgram::get_bounds_info(
    Index    n,
    Number * lb,
    Number * ub,
    Index    m,
    Number * cl,
    Number * cu
  ) {
    // Fill in the structures with the bounds information.
    std::copy_n(options.lowerbounds().begin(),n,lb);
    std::copy_n(options.upperbounds().begin(),n,ub);
    std::copy_n(options.constraintlb().begin(),m,cl);
    std::copy_n(options.constraintub().begin(),m,cu);
    return true;
  }

  bool
  MatlabProgram::get_starting_point(
    Index    n,
    bool     initializeVars,
    Number * vars,
    bool     initializez,
    Number * zl,
    Number * zu,
    Index    m,
    bool     initializeLambda,
    Number * lambda
  ) {

    IPOPT_ASSERT(
      n == numvars(options),
      "Bad number of variables required in get_starting_point, expected = " <<
      numvars(options) << " required = " << n
    );
    IPOPT_ASSERT(
      m == numconstraints(options),
      "Bad number of constraints required in get_starting_point, expected = " <<
      numconstraints(options) << " required = " << m
    );

    // Check to see whether IPOPT is requesting the initial point for the primal variables.
    if ( initializeVars ) std::copy_n(funcs.getx0().begin(),n,vars);

    // Check to see whether IPOPT is requesting the initial point for the Lagrange
    // multipliers associated with the bounds on the optimization variables.
    if ( initializez ) {
      IPOPT_ASSERT(
        options.multlb().size() > 0 && options.multub().size() > 0,
        "Initialization of Lagrange multipliers for lower and upper bounds "
        "requested but initial values are not supplied"
      );
      std::copy_n(options.multlb().begin(),n,zl);
      std::copy_n(options.multub().begin(),n,zu);
    }

    // Check to see whether IPOPT is requesting the initial point for the Lagrange
    // multipliers corresponding to the equality and inequality constraints.
    if ( initializeLambda && m>0 ) {
      IPOPT_ASSERT(
        options.multconstr().size() > 0,
        "Initialization of Lagrange multipliers for constraints are requested "
        "but initial values are not provided"
      );
      std::copy_n(options.multconstr().begin(),m,lambda);
    }

    return true;
  }

  bool
  MatlabProgram::eval_f(
    Index          n,
    Number const * vars,
    bool           ignore,
    Number       & f
  ) {
    f = funcs.computeObjective(n,vars);
    return true;
  }

  bool
  MatlabProgram::eval_grad_f(
    Index          n,
    Number const * vars,
    bool           ignore,
    Number       * grad
  ) {
    funcs.computeGradient( n, vars, grad );
    return true;
  }

  bool
  MatlabProgram::eval_g(
    Index          n,
    Number const * vars,
    bool           ignore,
    Index          m,
    Number       * g
  ) {
    if ( m > 0 ) {
      IPOPT_ASSERT(
        funcs.constraintFuncIsAvailable(),
        "You need to specify the callback function for computing the constraints"
      );
      funcs.computeConstraints(n,vars,m,g);
    }
    return true;
  }

  bool
  MatlabProgram::eval_jac_g(
    Index          numVariables,
    Number const * variables,
    bool           ignoreThis,
    Index          numConstraints,
    Index          Jx_nnz,
    Index        * rows,
    Index        * cols,
    Number       * Jx
  ) {
  
    if ( numConstraints > 0 ) {
      IPOPT_ASSERT(
        funcs.jacobianFuncIsAvailable(),
        "You need to specify the callback functions "
        "for computing the Jacobian and the sparsity structure of the Jacobian"
      );

      // If the input Jx is 0, then return the sparsity structure of
      // the Jacobian. Otherwise, return the values of the nonzero
      // entries.
      if ( Jx == nullptr ) {
        // Delete any previous structure information concerning the Jacobian matrix.
        // Get the sparse matrix structure of the Jacobian.
        // funcs.loadJacobianStructure(numVariables,numConstraints);
        IPOPT_ASSERT(
          funcs.getJacobianNnz() == Jx_nnz,
          "The constraint Jacobian passed back from the MATLAB "
          "routine has an incorrect number of nonzero entries"
        );
        // Copy the sparse matrix structure to the IPOPT sparse matrix format.
        funcs.getJacobianStructure(rows,cols);
      } else {
        // Return the value of the Jacobian.
        funcs.computeJacobian(numConstraints,numVariables,variables,Jx);
      }
    }
    return true;
  }

  bool
  MatlabProgram::eval_h(
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
  ) {
    // If the input Hx is 0, then return the sparsity structure of the
    // Hessian. Otherwise, return the values of the nonzero entries.
    if (Hx == nullptr) {
      IPOPT_ASSERT(
        funcs.hessianFuncIsAvailable(),
        "You need to specify the callback functions for computing "
        "the Hessian and the sparsity structure of the Hessian"
      );
      // Delete any previous structure information concerning the Hessian matrix.
      // Return the sparse matrix structure of the symmetric Hessian.
      //funcs.getHessianStructure(n,H);
      IPOPT_ASSERT(
        funcs.getHessianNnz() == Hx_nnz,
        "The Hessian passed back from the MATLAB routine "
        "has an incorrect number of nonzero entries"
      );
      // Copy the sparse matrix structure to the IPOPT sparse matrix format.
      funcs.getHessianStructure(rows,cols);
    } else {
      // Return the value of the lower triangular portion of the Hessian.
      funcs.computeHessian(n,vars,sigma,m,lambda,Hx);
    }
    return true;
  }

  void
  MatlabProgram::finalize_solution(
    SolverReturn   status,
    Index          numVariables,
    Number const * variables,
    Number const * zl,
    Number const * zu,
    Index          numConstraints,
    Number const * constraints,
    Number const * lambda,
    Number         objective,
    IpoptData const * ip_data,
    IpoptCalculatedQuantities* ip_cq
  ) {

    // Store the current solution the value of the Lagrange multipliers at the solution.
    info.setfield(numVariables,zl,"zl");
    info.setfield(numVariables,zu,"zu");
    info.setfield(numConstraints,lambda,"lambda");
    info.setfield(objective,"objective");
    // se necessario converto in cell array
    funcs.fillx( numVariables, variables );
    info.setfield(funcs.mx_getx(),"x");
  }

  bool
  MatlabProgram::intermediate_callback(
    AlgorithmMode mode,
    Index        iter,
    Number       f,
    Number       inf_pr,
    Number       inf_du,
    Number       mu,
    Number       d_norm,
    Number       regularization_size,
    Number       alpha_du,
    Number       alpha_pr,
    Index        ls_trials,
    IpoptData const * ip_data,
    IpoptCalculatedQuantities* ip_cq
  ) {

    if (funcs.iterFuncIsAvailable())
      return funcs.iterCallback(
        iter, f, inf_pr, inf_du, mu, d_norm,
        regularization_size, alpha_du, alpha_pr,
        ls_trials, ip_data, ip_cq
      );
    else
      return true;
  }
}

/*
//       _                              _
//      | | ___  _   _ _ __ _ __   __ _| |
//   _  | |/ _ \| | | | '__| '_ \ / _` | |
//  | |_| | (_) | |_| | |  | | | | (_| | |
//   \___/ \___/ \__,_|_|  |_| |_|\__,_|_|
*/

namespace Ipopt {

  #ifdef HAVE_VSNPRINTF
    #define IPOPT_VSNPRINTF vsnprintf
  #elif defined(HAVE__VSNPRINTF)
    #define IPOPT_VSNPRINTF _vsnprintf
  #endif

  #ifdef HAVE_VA_COPY
    #define VA_LIST_AP apcopy
  #else
    #define VA_LIST_AP ap
  #endif

  void
  MatlabJournal::PrintfImpl(
    EJournalCategory category,
    EJournalLevel    level,
    char const *     pformat,
    va_list          ap
  ) {

    Index const maxStrLen = 1024;
    char        s[maxStrLen];
    int         nchar;
    #ifdef IPOPT_VSNPRINTF
      #ifdef HAVE_VA_COPY
      va_list apcopy;
      va_copy(apcopy, ap);
      #endif
      nchar = IPOPT_VSNPRINTF(s,maxStrLen,pformat,VA_LIST_AP);
      #ifdef HAVE_VA_COPY
      va_end(apcopy);
      #endif
    #else
      nchar = vsprintf(s,pformat,ap);
    #endif
    mexPrintf("%s",s);
    //mexEvalString("drawnow;"); // to dump string.
    mexEvalString("pause(.001);"); // to dump string.

    IPOPT_ASSERT(
      nchar < maxStrLen,
      "String buffer it too short for all the characters to be printed to MATLAB console"
    );
  }
}
