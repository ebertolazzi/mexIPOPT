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

#ifdef __clang__
#pragma clang diagnostic ignored "-Wignored-attributes"
#endif

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

    IPOPT_DEBUG("In SparseMatrix::setup");

    IPOPT_ASSERT(
      mxIsSparse( ptr ),
      "Error in SparseMatrix::setup, expected sparse matrix"
    );

    // Get the height, width and number of non-zeros.
    m_numRows = Index( mxGetM(ptr) );
    m_numCols = Index( mxGetN(ptr) );
    // The precise number of non-zero elements is contained in the
    // last entry of the jc array. (There is one jc entry for each
    // column in the matrix, plus an extra one.)
    m_Jc.resize(m_numCols+1);
    mwIndex const * p = mxGetJc(ptr);
    for ( Index i=0; i <= m_numCols; ++i, ++p ) m_Jc[i] = Index(*p);
    //std::copy_n(mxGetJc(ptr),numCols+1,Jc.begin());
    // Copy the row and column indices, and the values of the nonzero entries.
    m_nnz = m_Jc.back();
    m_Ir.resize(m_nnz);
    p = mxGetIr(ptr);
    for ( Index i=0; i < m_nnz; ++i, ++p ) m_Ir[i] = Index(*p);
    // std::copy_n(mxGetIr(ptr),nnz,Ir.begin());
    IPOPT_DEBUG("Exit SparseMatrix::setup");
  }

  // extract the pattern of the sparse matrix in the ipopt format
  void
  SparseMatrix::getStructure( Index rows[], Index cols[] ) const {
    IPOPT_DEBUG("In SparseMatrix::getStructure");
    for ( Index c = 0, i = 0; c < m_numCols; ++c ) {
      Index iend = m_Jc[c+1];
      for (; i < iend; ++i ) { cols[i] = c; rows[i] = m_Ir[i]; }
    }
    IPOPT_DEBUG("Out SparseMatrix::getStructure");
  }

  // extyract values from the sparse matrix in ptr and store in values
  // which correspond to nonzeros of the sparse matrix described by
  // Ir and Jc. The nonzeros of sparse matrix in ptr must be a subset of the
  // nonzeros of the pattern described by Ir and Jc.
  void
  SparseMatrix::getValues( std::string const & func, mxArray * ptr, Number values[] ) const {

    IPOPT_DEBUG("In SparseMatrix::getValues");

    IPOPT_ASSERT(
      mxIsSparse( ptr ),
      "Error in SparseMatrix::getValues, expected sparse matrix"
    );

    // il patterm può essere un sottoinsieme
    mwIndex const * mxJc = mxGetJc(ptr);
    mwIndex const * mxIr = mxGetIr(ptr);
    double  const * v    = mxGetPr(ptr);

    std::fill_n( values, m_nnz, 0 );
    mwIndex i, k, i1, k1;
    for ( mwIndex c = 0; c < mwIndex(m_numCols); ++c ) {
      i = mxJc[c], i1 = mxJc[c+1];
      k = m_Jc[c]; k1 = m_Jc[c+1];
      for (; i < i1; ++i, ++k ) {
        mwIndex mxi = mxIr[i];
        while ( k < k1 && mwIndex(m_Ir[k]) < mxi ) ++k; // skip not set elements
        if ( k < k1 && mwIndex(m_Ir[k]) == mxi ) {
          IPOPT_ASSERT(
            std::isfinite(v[i]),
            "In MATLAB function " << func <<
            "\nelement (" << mxi+1 << "," << c+1 << ") is NaN\n"
          );
          values[k] = v[i];
        } else {
          IPOPT_DO_ERROR(
            "In MATLAB function " << func <<
            "\nelement (" << mxi+1 << "," << c+1 << ") not found in pattern"
          );
        }
      }
    }
    IPOPT_DEBUG("Out SparseMatrix::getValues");
  }

  static
  size_t
  from_matlab(
    mxArray const *       ptr,
    std::vector<Number> & data,
    char const            msg[]
  ) {

    IPOPT_DEBUG("In from_matlab");

    mwIndex n = mxGetNumberOfElements(ptr);
    if ( mxIsCell(ptr) ) {

      // Compute the number of optimization variables.
      mwIndex nv = 0;
      for ( mwIndex i = 0; i < n; ++i ) {
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

      for ( mwIndex i = 0; i < n; ++i ) {
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

    IPOPT_DEBUG("Out from_matlab");

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

    IPOPT_DEBUG("In isLowerTri");

    IPOPT_ASSERT(
      mxIsSparse( ptr ),
      "Error in isLowerTri, expected sparse matrix"
    );

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
    for ( mwIndex j = 0; j < ncol && ok; ++j ) {
      mwIndex jcend = jc[j+1];
      for (; i < jcend && ok; ++i )
        ok = (j <= ir[i]);
    }

    IPOPT_DEBUG("Out isLowerTri");

    return ok;
  }

  static
  bool
  IsInIncOrder( mxArray const * ptr ) {

    IPOPT_DEBUG("In IsInIncOrder");

    IPOPT_ASSERT(
      mxIsSparse( ptr ),
      "Error in IsInIncOrder, expected sparse matrix"
    );

    Index    ncol = (Index) mxGetN(ptr);
    mwIndex* jc   = mxGetJc(ptr);
    mwIndex* ir   = mxGetIr(ptr);
    bool     ok   = true;
    for ( Index j = 0; j < ncol && ok; ++j ) {
      mwIndex jend = jc[j+1];
      mwIndex i    = jc[j];
      for ( ++i; i < jend && ok; ++i ) ok = ( ir[i-1] < ir[i] );
    }

    IPOPT_DEBUG("Out IsInIncOrder");

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
  MatlabFunctionHandle::bind( mxArray const * p, char const error_msg[] ) {

    IPOPT_DEBUG("In MatlabFunctionHandle::bind");

    IPOPT_ASSERT(
      p != nullptr,
      "You must specify a callback routine for computing " << error_msg
    );
    IPOPT_ASSERT(
      !mxIsEmpty(p) && mxIsClass(p,"function_handle"),
      "You did not provide a valid function handle for computing " << error_msg
    );

    if ( m_f != nullptr ) mxDestroyArray( m_f );
    m_f = mxDuplicateArray( p );

    mxArray *outputs[1];
    mxArray *exception = mexCallMATLABWithTrap( 1, outputs, 1, &m_f, "func2str" );

    if ( exception != nullptr  ) {
      mxArray *msg;
      mexCallMATLAB(1, &msg, 1, &exception, "getReport" );
      char * c_msg = mxArrayToString(msg);
      std::string cpp_msg = c_msg;
      mxFree(c_msg);
      IPOPT_DO_ERROR( "in function: func2str\n" << cpp_msg );
    }

    char * c_msg = mxArrayToString( outputs[0] );
    m_name = c_msg;
    mxFree(c_msg);

    IPOPT_DEBUG("Out MatlabFunctionHandle::bind");
  }

  void
  MatlabFunctionHandle::eval(
    Index            n_lhs,
    mxArray       ** lhs,
    Index            n_rhs,
    mxArray const ** rhs
  ) const {

    IPOPT_DEBUG("In MatlabFunctionHandle::eval");

    // Construct the inputs to "feval".
    std::vector<mxArray*> finputs(n_rhs+1);
    finputs[0] = m_f;
    for ( Index i = 0; i < n_rhs; ++i )
      finputs[i+1] = mxDuplicateArray(rhs[i]);

    // Call "feval".
    mxArray *exception = mexCallMATLABWithTrap(
      n_lhs, lhs, n_rhs+1, &finputs.front(), "feval"
    );

    // Free the dynamically allocated memory.
    for ( Index i = 1; i <= n_rhs; ++i )
      mxDestroyArray(finputs[i]);

    if ( exception != nullptr  ) {
      mxArray *msg;
      mexCallMATLAB(1, &msg, 1, &exception, "getReport" );
      char * c_msg = mxArrayToString(msg);
      std::string cpp_msg = c_msg;
      mxFree(c_msg);
      IPOPT_DO_ERROR( "in function: " << m_name << "\n" << cpp_msg );
    }

    IPOPT_DEBUG("Out MatlabFunctionHandle::eval");
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
  MatlabInfo::MatlabInfo( mxArray *& ptr ) : m_info_ptr(nullptr) {
    // Create the structure array.
    IPOPT_DEBUG("In MatlabInfo::MatlabInfo");

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
    m_info_ptr = ptr = mxCreateStructMatrix(1,1,9,fieldnames);

    // Initialize fields.
    for ( Index i = 0; i < 8; ++i ) // do not initialize "eval"
      mxSetField( m_info_ptr, 0, fieldnames[i], mxCreateDoubleScalar(0) );

    // Build Eval Structure
    static char const * evalfields[5] = {
      "objective",
      "constraints",
      "gradient",
      "jacobian",
      "hessian"
    };

    mxArray *evalStruct = mxCreateStructMatrix( 1, 1, 5, evalfields );
    for ( Index i = 0; i < 5; ++i )
      mxSetField( evalStruct, 0, evalfields[i], mxCreateDoubleScalar(0) );

    mxSetField( m_info_ptr, 0, "eval", evalStruct );

    IPOPT_DEBUG("Out MatlabInfo::MatlabInfo");
  }

  // The destructor.
  MatlabInfo::~MatlabInfo() {
    IPOPT_DEBUG("In MatlabInfo::~MatlabInfo");
  }

  void
  MatlabInfo::setfield( mxArray const * pin, char const * field ) {

    IPOPT_DEBUG("In MatlabInfo::setfield");

    mxArray * p = mxGetField( m_info_ptr, 0, field );
    IPOPT_ASSERT(
      p != nullptr, "MatlabInfo::setfield missing field `" << field << "'"
    );
    // First destroy any previous values.
    mxDestroyArray( p );
    mxSetField( m_info_ptr, 0, field, mxDuplicateArray( pin ));

    IPOPT_DEBUG("Out MatlabInfo::setfield");
  }

  void
  MatlabInfo::setfield( size_t n, Number const * x, char const * field ) {

    IPOPT_DEBUG("In MatlabInfo::setfield2");

    mxArray * p = mxGetField( m_info_ptr, 0, field );
    IPOPT_ASSERT(
      p != nullptr, "MatlabInfo::setfield missing field `" << field << "'"
    );
    // First destroy any previous values.
    mxDestroyArray(p);
    p = mxCreateDoubleMatrix( n, 1, mxREAL );
    std::copy_n( x, n, mxGetPr(p) );
    mxSetField( m_info_ptr, 0, field, p );

    IPOPT_DEBUG("Out MatlabInfo::setfield2");
  }

  mxArray *
  MatlabInfo::getfield_mx( char const * field ) const {

    IPOPT_DEBUG("In MatlabInfo::getfield_mx");

    mxArray * p = mxGetField( m_info_ptr, 0, field );
    IPOPT_ASSERT(
      p != nullptr,
      "MatlabInfo::getfield_mx missing field `" << field << "'"
    );

    IPOPT_DEBUG("Out MatlabInfo::getfield_mx");

    return p;
  }

  Number const *
  MatlabInfo::getfield( char const * field ) const {

    IPOPT_DEBUG("In MatlabInfo::getfield");

    Number const * p = mxGetPr(getfield_mx(field));
    IPOPT_ASSERT(
      p != nullptr,
      "MatlabInfo::getfield('" << field << "') cant access values"
    );

    IPOPT_DEBUG("Out MatlabInfo::getfield");

    return p;
  }

  ApplicationReturnStatus
  MatlabInfo::getExitStatus() const {

    IPOPT_DEBUG("In MatlabInfo::getExitStatus");

    mxArray const * p = getfield_mx("status");

    IPOPT_DEBUG("Out MatlabInfo::getExitStatus");
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

    IPOPT_DEBUG("In MatlabInfo::setFuncEvals");

    mxArray* p = getfield_mx("eval");
    *mxGetPr( mxGetField( p, 0, "objective") )   = obj;
    *mxGetPr( mxGetField( p, 0, "constraints") ) = con;
    *mxGetPr( mxGetField( p, 0, "gradient") )    = grad;
    *mxGetPr( mxGetField( p, 0, "jacobian") )    = jac;
    *mxGetPr( mxGetField( p, 0, "hessian") )     = hess;

    IPOPT_DEBUG("Out MatlabInfo::setFuncEvals");

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
  IpoptOptions::IpoptOptions( IpoptApplication & app, mxArray const * ptr ) : m_app(app) {

    IPOPT_DEBUG("In IpoptOptions::IpoptOptions");

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

    IPOPT_DEBUG("Out IpoptOptions::IpoptOptions");
  }

  // -----------------------------------------------------------------
  IpoptOptions::~IpoptOptions( ) {
    IPOPT_DEBUG("In IpoptOptions::~IpoptOptions");
  }

  bool
  IpoptOptions::useQuasiNewton() const {

    IPOPT_DEBUG("In IpoptOptions::useQuasiNewton");

    std::string value;
    m_app.Options()->GetStringValue("hessian_approximation",value,"");
    bool b = !value.compare("limited-memory");

    IPOPT_DEBUG("Out IpoptOptions::useQuasiNewton");

    return b;
  }

  bool
  IpoptOptions::useDerivChecker() const {

    IPOPT_DEBUG("In IpoptOptions::useDerivChecker");

    std::string value;
    m_app.Options()->GetStringValue("derivative_test",value,"");
    bool b = value.compare("none");

    IPOPT_DEBUG("Out IpoptOptions::useDerivChecker");

    return b;
  }

  bool
  IpoptOptions::userScaling() const {

    IPOPT_DEBUG("In IpoptOptions::userScaling");

    std::string value;
    m_app.Options()->GetStringValue("nlp_scaling_method",value,"");
    bool b = !value.compare("user-scaling");

    IPOPT_DEBUG("Out IpoptOptions::userScaling");

    return b;
  }

  Index
  IpoptOptions::printLevel() const {

    IPOPT_DEBUG("In IpoptOptions::printLevel");

    Index value;  // The return value.
    m_app.Options()->GetIntegerValue("print_level",value,"");

    IPOPT_DEBUG("Out IpoptOptions::printLevel");

    return value;
  }

  void
  IpoptOptions::setOption( char const * label, mxArray const * ptr ) {

    IPOPT_DEBUG("In IpoptOptions::setOption");

    // Check to make sure we have a valid option.
    SmartPtr<RegisteredOption const> option = m_app.RegOptions()->GetOption(label);
    IPOPT_ASSERT(
      IsValid(option),
      "You have specified a nonexistent IPOPT option (\"" << label << "\")"
    );
    Ipopt::RegisteredOptionType type = option->Type();
    if      (type == Ipopt::OT_String)  setStringOption(label,ptr);
    else if (type == Ipopt::OT_Integer) setIntegerOption(label,ptr);
    else                                setNumberOption(label,ptr);

    IPOPT_DEBUG("Out IpoptOptions::setOption");
  }

  void
  IpoptOptions::setStringOption( char const * label, mxArray const * ptr ) {

    IPOPT_DEBUG("In IpoptOptions::setStringOption");

    // Check whether the option value is a string.
    IPOPT_ASSERT(
      mxIsChar(ptr),
      "IPOPT option value for option \"" << label << "\" should be a string"
    );

    // Get the option value.
    char* value = mxArrayToString(ptr);

    // Set the option.
    bool success = m_app.Options()->SetStringValue(label,value);
    IPOPT_ASSERT(
      success, "Invalid value for IPOPT option \"" << label << "\""
    );

    // Free the dynamically allocated memory.
    mxFree(value);

    IPOPT_DEBUG("Out IpoptOptions::setStringOption");
  }

  void
  IpoptOptions::setIntegerOption( char const * label, mxArray const * ptr) {

    IPOPT_DEBUG("In IpoptOptions::setIntegerOption");

    // Check whether the option value is a number.
    IPOPT_ASSERT(
      mxIsDouble(ptr),
      "IPOPT option value for option \"" << label << "\" should be an integer"
    );

    // Set either the integer option.
    Number value   = mxGetScalar(ptr);
    bool   success = m_app.Options()->SetIntegerValue(label,(Index) value);
    IPOPT_ASSERT(
      success, "Invalid value for integer IPOPT option \"" << label << "\""
    );

    IPOPT_DEBUG("Out IpoptOptions::setIntegerOption");
  }

  void
  IpoptOptions::setNumberOption( char const * label, mxArray const * ptr ) {

    IPOPT_DEBUG("In IpoptOptions::setNumberOption");

    // Check whether the option value is a number.
    IPOPT_ASSERT(
      mxIsDouble(ptr),
      "IPOPT option value for option \"" << label << "\" should be a number"
    );

    // Set either the numeric option.
    Number value   = mxGetScalar(ptr);
    bool   success = m_app.Options()->SetNumericValue(label,value);
    IPOPT_ASSERT(
      success, "Invalid value for numeric IPOPT option \"" << label << "\""
    );

    IPOPT_DEBUG("Out IpoptOptions::setNumberOption");
  }

  // Function definitions for class Options.
  // -----------------------------------------------------------------
  Options::Options( Index n_in, Ipopt::IpoptApplication& app, mxArray const * ptr)
  : m_n(n_in)
  , m_m(0) // will be set by loadConstraintBounds
  , m_ipopt( app, mxGetField(ptr,0,"ipopt") ) // Process the IPOPT options.
  {

    IPOPT_DEBUG("In Options::Options");

    app.Options()->GetNumericValue("nlp_upper_bound_inf",m_posinfty,"");
    app.Options()->GetNumericValue("nlp_lower_bound_inf",m_neginfty,"");

    // Load the bounds on the variables.
    loadLowerBounds(ptr);
    loadUpperBounds(ptr);

    // Load the bounds on the constraints.
    loadConstraintBounds(ptr);

    // Load the Lagrange multipliers.
    loadMultipliers(ptr);

    IPOPT_DEBUG("Out Options::Options");
  }

  Options::~Options() {
    IPOPT_DEBUG("In Options:~Options");
  }

  // Function definitions for static members of class Options.
  // -----------------------------------------------------------------
  void
  Options::loadLowerBounds( mxArray const * ptr ) {

    IPOPT_DEBUG("In Options:loadLowerBounds");

    mxArray const * p = mxGetField( ptr, 0, "lb" );
    if ( p != nullptr ) {
      // Load the upper bounds and check to make sure they are valid.
      Index N = Index( from_matlab( p, m_lb, "Options::loadLowerBounds" ) );
      IPOPT_ASSERT(
        N == m_n,
        "Lower bounds array must have one element for each optimization"
        " variable, found a vector of N = " << N <<
        " elements, expected n = " << m_n
      );
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      mxToIpoptInf( m_lb, m_neginfty );
    } else {
      // If the lower bounds have not been specified, set them to negative infinity.
      m_lb.resize( m_n );
      std::fill( m_lb.begin(), m_lb.end(), m_neginfty );
    }

    IPOPT_DEBUG("Out Options:loadLowerBounds");
  }

  void
  Options::loadUpperBounds( mxArray const * ptr ) {

    IPOPT_DEBUG("In Options:loadUpperBounds");

    mxArray const * p = mxGetField( ptr, 0, "ub" );
    if ( p != nullptr ) {
      // Load the upper bounds and check to make sure they are valid.
      Index N = Index( from_matlab( p, m_ub, "Options::loadUpperBounds" ) );
      IPOPT_ASSERT(
        N == m_n,
        "Upper bounds array must have one element for each optimization variable"
      );
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      mxToIpoptInf( m_ub, m_posinfty );
    } else {
      // If the upper bounds have not been specified, set them to positive infinity.
      m_ub.resize( m_n );
      std::fill( m_ub.begin(), m_ub.end(), m_posinfty );
    }

    IPOPT_DEBUG("Out Options:loadUpperBounds");
  }

  void
  Options::loadConstraintBounds( mxArray const * ptr ) {

    IPOPT_DEBUG("In Options:loadConstraintBounds");

    m_m = 0;  // The return value is the number of constraints.
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
      m_m = (Index) mxGetNumberOfElements(pl);

      // Load the lower bounds on the constraints and convert MATLAB's
      // convention of infinity to IPOPT's convention of infinity.
      // Convert MATLAB's convention of infinity to IPOPT's convention of infinity.
      m_cl.resize(m_m);
      std::copy_n(mxGetPr(pl),m_m,m_cl.begin());
      mxToIpoptInf( m_cl, m_neginfty );

      m_cu.resize(m_m);
      std::copy_n(mxGetPr(pu),m_m,m_cu.begin());
      mxToIpoptInf( m_cu, m_posinfty );
    }

    IPOPT_DEBUG("Out Options:loadConstraintBounds");
  }

  void
  Options::loadMultipliers( mxArray const * ptr ) {

    IPOPT_DEBUG("In Options:loadMultipliers");

    // Load the Lagrange multipliers associated with the lower bounds.
    mxArray const * p = mxGetField(ptr,0,"zl");
    static char const * msg1 = "must be a double-precision array with one element\nfor each optimization variable";
    if ( p != nullptr ) {
      static char const * msg = "The initial point for the Lagrange multipliers associated with the lower bounds\n";
      IPOPT_ASSERT( mxIsDouble(p) && !mxIsComplex(p), msg << msg1 );
      IPOPT_ASSERT(
        Index(mxGetNumberOfElements(p)) == Index(m_n),
        msg << "must be a vector of length " << m_n
            << " found of lenght " << mxGetNumberOfElements(p)
      );
      m_zl.resize( m_n );
      std::copy_n( mxGetPr(p), m_n, m_zl.begin() );
    } else {
      m_zl.clear();
    }
    // Load the Lagrange multipliers associated with the upper bounds.
    p = mxGetField(ptr,0,"zu");
    if ( p != nullptr ) {
      static char const * msg = "The initial point for the Lagrange multipliers associated with the upper bounds\n";
      IPOPT_ASSERT( mxIsDouble(p) && !mxIsComplex(p), msg << msg1 );
      IPOPT_ASSERT(
        Index(mxGetNumberOfElements(p)) == Index( m_n ),
        msg << "must be a vector of length " << m_n
            << " found of lenght " << mxGetNumberOfElements(p)
      );
      m_zu.resize( m_n );
      std::copy_n( mxGetPr(p), m_n, m_zu.begin() );
    } else {
      m_zu.clear();
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
        Index(mxGetNumberOfElements(p)) == Index( m_m ),
        msg << "must be a vector of length " << m_n
            << " found of length " << mxGetNumberOfElements(p)
      );
      m_lambda.resize( m_m );
      std::copy_n( mxGetPr(p), m_m, m_lambda.begin() );
    } else {
      m_lambda.clear();
    }

    IPOPT_DEBUG("Out Options:loadMultipliers");
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

    IPOPT_DEBUG("In CallbackFunctions::CallbackFunctions");

    from_matlab( mx_x0, m_x0, "in CallbackFunctions" );

    mxArray const * p;  // A pointer to a MATLAB array.
    m_x_is_cell_array = mxIsCell( mx_x0 );
    if ( m_x_is_cell_array ) {
      mx_x_nv = 0;
      mx_x_nc = Index( mxGetNumberOfElements( mx_x0 ));
      mx_x    = mxCreateCellMatrix( mx_x_nc, 1 );
      for ( Index i = 0; i < mx_x_nc; ++i ) {
        p = mxGetCell( mx_x0, i );
        mx_x_nv += (Index) mxGetNumberOfElements( p );
        mxSetCell( mx_x, i, mxDuplicateArray( p ) );
      }
    } else {
      mx_x    = mxDuplicateArray( mx_x0 );
      mx_x_nv = Index( mxGetNumberOfElements( mx_x ) );
      mx_x_nc = 0;
    }

    // Check whether we are provided with a structure array.
    IPOPT_ASSERT(
      mxIsStruct(ptr),
      "in CallbackFunctions the second input must be a STRUCT, found: " <<
      mxGetClassName(ptr)
    );

    // Get the function handle for computing the objective.
    p = mxGetField( ptr, 0, "objective" );
    m_obj.bind( p, "the objective function" );

    // Get the function handle for computing the gradient.
    p = mxGetField( ptr, 0, "gradient" );
    m_grad.bind( p, "the gradient of the objective" );

    // Get the function handle for computing the constraints, if such a function was specified.
    p = mxGetField( ptr, 0, "constraints" );
    if ( p != nullptr ) m_constraint.bind( p, "the response of the constraints" );

    // Get the function handle for computing the Jacobian.
    // This function is necessary if there are constraints.
    p = mxGetField( ptr, 0, "jacobian" );
    if ( p != nullptr ) m_jacobian.bind( p, "the first derivatives (Jacobian) of the constraints" );

    // Get the function handle for computing the sparsity structure of the Jacobian.
    // This function is necessary if the Jacobian is being computed.
    p = mxGetField( ptr, 0, "jacobianstructure" );
    if ( p != nullptr ) m_jacstruc.bind(p,"the sparsity structure of the Jacobian");

    // Get the function handle for computing the Hessian. This function is always optional.
    p = mxGetField( ptr, 0, "hessian" );
    if ( p != nullptr ) m_hessian.bind(p,"the Hessian of the Lagrangian");

    // Get the function handle for computing the sparsity structure of the Hessian of the Lagrangian.
    // This function is necessary if the Hessian is being computed.
    p = mxGetField( ptr, 0, "hessianstructure" );
    if ( p != nullptr ) m_hesstruc.bind( p,"the sparsity structure of the Hessian" );

    // Get the iterative callback function handle. This function is always optional.
    p = mxGetField( ptr, 0, "iterfunc" );
    if ( p != nullptr ) m_iter.bind( p, "the iterative callback" );

    IPOPT_DEBUG("Out CallbackFunctions::CallbackFunctions");
  }

  CallbackFunctions::~CallbackFunctions() {
    IPOPT_DEBUG("In CallbackFunctions::~CallbackFunctions");
  }

  bool
  CallbackFunctions::from_cell_array(
    mxArray const * ptr,
    Index           n,
    Number        * x
  ) const {
    IPOPT_DEBUG("In CallbackFunctions::from_cell_array");
    if ( !m_x_is_cell_array ) return false;
    Index ntot = 0;
    for ( Index i = 0; i < this -> mx_x_nc; ++i ) {
      mxArray const * p = mxGetCell(ptr,i);
      Index nn = (Index) mxGetNumberOfElements(p);
      ntot += nn;
      if ( ntot > n ) return false;
      std::copy_n(mxGetPr(p),nn,x);
      x += nn;
    }
    IPOPT_DEBUG("Out CallbackFunctions::from_cell_array");
    return ntot == n;
  }

  void
  CallbackFunctions::fillx( Index n, Number const * x ) const {
    IPOPT_DEBUG("In CallbackFunctions::fillx");
    if ( m_x_is_cell_array ) {
      for ( Index i = 0; i < mx_x_nc; ++i ) {
        mxArray const * p = mxGetCell(mx_x,i);
        Index nn = (Index) mxGetNumberOfElements(p);
        std::copy_n( x, nn, mxGetPr(p) );
        x += nn;
      }
    } else {
      std::copy_n( x, n, mxGetPr(mx_x) );
    }
    IPOPT_DEBUG("Out CallbackFunctions::fillx");
  }

  Number
  CallbackFunctions::computeObjective( Index n, Number const x[] ) const {

    IPOPT_DEBUG("In CallbackFunctions::computeObjective");

    IPOPT_ASSERT(
      m_obj.ok(),
      "You must provide a function that compute the objective function"
    );

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_obj.eval( 1, &ptr, 1, (mxArray const **)&mx_x);

    // Get the output from the MATLAB callback function, which is the
    // value of the objective function at x.
    IPOPT_ASSERT(
      mxIsDouble(ptr) && !mxIsComplex(ptr),
      "In MATLAB function " << m_obj.name() << "\n"
      "The first return value of the objective callback function"
      " must be a real double scalar, found: " << mxGetClassName(ptr)
    );
    IPOPT_ASSERT(
      Index(mxGetNumberOfElements(ptr)) == Index(1),
      "In MATLAB function " << m_obj.name() << "\n"
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

    IPOPT_DEBUG("Out CallbackFunctions::computeObjective");

    return f;
  }

  void
  CallbackFunctions::computeGradient(
    Index        n,
    Number const x[],
    Number       g[]
  ) const {

    IPOPT_DEBUG("In CallbackFunctions::computeGradient");

    IPOPT_ASSERT(
      m_grad.ok(),
      "You must provide a function that compute the gradient of the objective function"
    );

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_grad.eval( 1, &ptr, 1, (mxArray const **)&mx_x );

    // se cell array converto
    if ( mxIsCell(ptr) ) {
      IPOPT_ASSERT(
        from_cell_array( ptr, n, g ),
        "In MATLAB function " << m_grad.name() << "\n"
        "The gradient callback return a cell array not with the same structure of x0"
      );
    } else {
      // Get the output from the MATLAB callback function, which is the
      // value of the gradient of the objective function at x.
      IPOPT_ASSERT(
        mxIsDouble(ptr) && !mxIsComplex(ptr),
        "In MATLAB function " << m_grad.name() << "\n"
        "The gradient callback must return a real double vector, found: " <<
        mxGetClassName(ptr)
      );
      if ( mxIsSparse(ptr) ) {
        mxArray * ptr1;
        // convert sparse gradient to full (simplest method, not fastest)
        mexCallMATLAB( 1, &ptr1, 1, &ptr, "full" );
        std::swap( ptr, ptr1 );
        mxDestroyArray( ptr1 ); // destroy old sparse vector
      }
      IPOPT_ASSERT(
        Index(mxGetNumberOfElements(ptr)) == Index(n),
        "In MATLAB function " << m_grad.name() << "\n"
        "The gradient callback must return a real double vector of size = " << n <<
        " while it return a vector of size = " << mxGetNumberOfElements(ptr)
      );

      std::copy_n( mxGetPr(ptr), n, g );
    }

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);

    IPOPT_DEBUG("Out CallbackFunctions::computeGradient");
  }

  void
  CallbackFunctions::computeConstraints(
    Index        n,
    Number const x[],
    Index        m,
    Number       c[]
  ) const {

    IPOPT_DEBUG("In CallbackFunctions::computeConstraints");

    IPOPT_ASSERT(
      m_grad.ok(),
      "You must provide a function that compute the constraints"
    );

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_constraint.eval( 1, &ptr, 1, (mxArray const **)&mx_x );

    // Get the output from the MATLAB callback function, which is the
    // value of vector-valued constraint function at x.
    IPOPT_ASSERT(
      mxIsDouble(ptr) && !mxIsComplex(ptr),
      "In MATLAB function " << m_constraint.name() << "\n"
      "The constraint callback must return a real double vector, found: " <<
      mxGetClassName(ptr)
    );

    IPOPT_ASSERT(
      Index(mxGetNumberOfElements(ptr)) == Index(m),
      "In MATLAB function " << m_constraint.name() << "\n"
      "The constraints callback must return a real double vector of size = " << m <<
      " while it return a vector of size = " << mxGetNumberOfElements(ptr)
    );

    if (mxIsSparse(ptr)) {
      mxArray * ptr1;
      // convert sparse constraint vector (unlikely but possible) to full
      mexCallMATLAB( 1, &ptr1, 1, &ptr, "full" );
      std::swap( ptr, ptr1 );
      mxDestroyArray( ptr1 ); // destroy old sparse vector
    }
    std::copy_n( mxGetPr(ptr), m, c );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);

    IPOPT_DEBUG("Out CallbackFunctions::computeConstraints");
  }

  static
  void
  checkJacobian( std::string const & name, Index n, Index m, mxArray * ptr ) {

    IPOPT_DEBUG("In checkJacobian");

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
      Index(mxGetM(ptr)) == Index(m) && Index(mxGetN(ptr)) == Index(n),
      "In MATLAB function " << name << "\n"
      "Jacobian must be an (m=" << m << ") x (n=" << n << ") sparse matrix\n"
      "where m is the number of constraints and n is the number of variables,\n"
      "but found m=" << mxGetM(ptr) << " and n=" << mxGetN(ptr)
    );
    IPOPT_ASSERT(
      IsInIncOrder(ptr),
      "In MATLAB function " << name << "\n"
      "Jacobian must be a sparse matrix with row indices in increasing order"
    );

    IPOPT_DEBUG("Out checkJacobian");
  }

  void
  CallbackFunctions::loadJacobianStructure( Index n, Index m ) const {

    IPOPT_DEBUG("In CallbackFunctions::loadJacobianStructure");

    IPOPT_ASSERT(
      m_jacstruc.ok(),
      "You must provide a function that returns the sparsity structure of the Jacobian"
    );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_jacstruc.eval( 1, &ptr, 0, (mxArray const **)nullptr );

    checkJacobian( m_jacstruc.name(), n, m, ptr );
    m_Jacobian.setup( ptr );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);

    IPOPT_DEBUG("Out CallbackFunctions::loadJacobianStructure");
  }

  static
  void
  checkHessian( std::string const & name, Index n, mxArray * ptr ) {

    IPOPT_DEBUG("In checkHessian");

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
      Index(mxGetM(ptr)) == Index(n) && Index(mxGetN(ptr)) == Index(n),
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
      IsInIncOrder(ptr),
      "In MATLAB function " << name << "\n"
      "Hessian must be an n x n sparse, symmetric and lower triangular matrix "
      "with row indices in increasing order.\n"
      "Row indices are not in increasing order"
    );

    IPOPT_DEBUG("Out checkHessian");
  }

  void
  CallbackFunctions::loadHessianStructure( Index n ) const {

    IPOPT_DEBUG("In CallbackFunctions::loadHessianStructure");

    IPOPT_ASSERT(
      m_hesstruc.ok(),
      "You must provide a function that returns the sparsity structure of the Hessian"
    );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_hesstruc.eval( 1, &ptr, 0, (mxArray const **)nullptr );

    checkHessian( m_hesstruc.name(), n, ptr );
    m_Hessian.setup( ptr );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);

    IPOPT_DEBUG("Out CallbackFunctions::loadHessianStructure");
  }

  void
  CallbackFunctions::computeJacobian(
    Index        m,
    Index        n,
    Number const x[],
    Number       values[]
  ) const {

    IPOPT_DEBUG("In CallbackFunctions::computeJacobian");

    IPOPT_ASSERT(
      m_jacobian.ok(),
      "You must provide a function that returns the first derivatives"
      " (Jacobian) of the constraints"
    );

    fillx( n, x );

    // Call the MATLAB callback function.
    mxArray * ptr;
    m_jacobian.eval( 1, &ptr, 1, (mxArray const **)&mx_x );

    checkJacobian( m_jacobian.name(), n, m, ptr );

    if (!mxIsDouble(ptr)) {
      mxArray * ptr1;
      // convert non-double Jacobian to double
      mexCallMATLAB(1, &ptr1, 1, &ptr, "double");
      std::swap(ptr,ptr1);
      mxDestroyArray(ptr1);
    }

    m_Jacobian.getValues( m_jacobian.name(), ptr, values );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);

    IPOPT_DEBUG("Out CallbackFunctions::computeJacobian");
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

    IPOPT_DEBUG("In CallbackFunctions::computeHessian");

    IPOPT_ASSERT(
      m_hessian.ok(),
      "You must provide a function that returns the Hessian"
    );

    fillx( n, x );

    // Create the input arguments to the MATLAB routine, sigma and lambda.
    mxArray * psigma  = mxCreateDoubleScalar(sigma);
    mxArray * plambda = mxCreateDoubleMatrix(m,1,mxREAL);
    std::copy_n( lambda, m, mxGetPr(plambda) );

    // Call the MATLAB callback function.
    mxArray const * inputs[3] = { mx_x, psigma, plambda };
    mxArray       * ptr;
    m_hessian.eval( 1, &ptr, 3, inputs );

    checkHessian( m_hessian.name(), n, ptr );
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

    m_Hessian.getValues( m_hessian.name(), ptr, values );

    // Free the dynamically allocated memory.
    mxDestroyArray(ptr);
    mxDestroyArray(psigma);
    mxDestroyArray(plambda);

    IPOPT_DEBUG("Out CallbackFunctions::computeHessian");
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

    IPOPT_DEBUG("In CallbackFunctions::iterCallback");

    if ( !m_iter.ok() ) return true;

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

    m_iter.eval(1,outputs,3,inputs);

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

    IPOPT_DEBUG("Out CallbackFunctions::iterCallback");

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
  : m_funcs(funcs)
  , m_options(options)
  , m_info(info) {
    IPOPT_DEBUG("In MatlabProgram::MatlabProgram");
  }

  MatlabProgram::~MatlabProgram() {
    IPOPT_DEBUG("In MatlabProgram::~MatlabProgram");
  }

  bool
  MatlabProgram::get_nlp_info(
    Index          & n,
    Index          & m,
    Index          & J_nnz,
    Index          & H_nnz,
    IndexStyleEnum & indexStyle
  ) {

    IPOPT_DEBUG("In MatlabProgram::get_nlp_info");

    // Get the number of variables and constraints.
    n = numvars(m_options);
    m = numconstraints(m_options);

    // Get the size of the Jacobian.
    if ( m > 0 ) {
      IPOPT_ASSERT(
        m_funcs.jacobianFuncIsAvailable(),
        "You need to specify the callback functions for computing the Jacobian "
        "and the sparsity structure of the Jacobian"
      );
      m_funcs.loadJacobianStructure(n,m);
      J_nnz = m_funcs.getJacobianNnz();
    } else {
      J_nnz = 0;
    }
    // Get the size of the symmetric Hessian matrix. We don't need to
    // store the actual result, we just need to look at the number of
    // non-zero entries in the lower triangular part of the matrix.
    if ( !m_options.ipoptOptions().useQuasiNewton() ) {
      IPOPT_ASSERT(
        m_funcs.hessianFuncIsAvailable(),
        "You need to specify the callback functions for computing the Hessian "
        "and the sparsity structure of the Hessian"
      );
      m_funcs.loadHessianStructure(n);
      H_nnz = m_funcs.getHessianNnz();
    } else {
      H_nnz = 0;
    }
    // Use C-style indexing.
    indexStyle = C_STYLE;

    IPOPT_DEBUG("Out MatlabProgram::get_nlp_info");

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

    IPOPT_DEBUG("In MatlabProgram::get_bounds_info");

    // Fill in the structures with the bounds information.
    std::copy_n( m_options.lowerbounds().begin(),  n, lb );
    std::copy_n( m_options.upperbounds().begin(),  n, ub );
    std::copy_n( m_options.constraintlb().begin(), m, cl );
    std::copy_n( m_options.constraintub().begin(), m, cu );

    IPOPT_DEBUG("Out MatlabProgram::get_bounds_info");

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

    IPOPT_DEBUG("In MatlabProgram::get_starting_point");

    IPOPT_ASSERT(
      n == numvars(m_options),
      "Bad number of variables required in get_starting_point, expected = " <<
      numvars(m_options) << " required = " << n
    );
    IPOPT_ASSERT(
      m == numconstraints(m_options),
      "Bad number of constraints required in get_starting_point, expected = " <<
      numconstraints(m_options) << " required = " << m
    );

    // Check to see whether IPOPT is requesting the initial point for the primal variables.
    if ( initializeVars ) std::copy_n( m_funcs.getx0().begin(), n, vars );

    // Check to see whether IPOPT is requesting the initial point for the Lagrange
    // multipliers associated with the bounds on the optimization variables.
    if ( initializez ) {
      IPOPT_ASSERT(
        m_options.multlb().size() > 0 && m_options.multub().size() > 0,
        "Initialization of Lagrange multipliers for lower and upper bounds "
        "requested but initial values are not supplied"
      );
      std::copy_n( m_options.multlb().begin(), n, zl );
      std::copy_n( m_options.multub().begin(), n, zu );
    }

    // Check to see whether IPOPT is requesting the initial point for the Lagrange
    // multipliers corresponding to the equality and inequality constraints.
    if ( initializeLambda && m>0 ) {
      IPOPT_ASSERT(
        m_options.multconstr().size() > 0,
        "Initialization of Lagrange multipliers for constraints are requested "
        "but initial values are not provided"
      );
      std::copy_n( m_options.multconstr().begin(), m, lambda );
    }

    IPOPT_DEBUG("Out MatlabProgram::get_starting_point");

    return true;
  }

  bool
  MatlabProgram::eval_f(
    Index          n,
    Number const * vars,
    bool           ignore,
    Number       & f
  ) {
    IPOPT_DEBUG("In MatlabProgram::eval_f");
    f = m_funcs.computeObjective( n, vars );
    IPOPT_DEBUG("Out MatlabProgram::eval_f");
    return true;
  }

  bool
  MatlabProgram::eval_grad_f(
    Index          n,
    Number const * vars,
    bool           ignore,
    Number       * grad
  ) {
    IPOPT_DEBUG("In MatlabProgram::eval_grad_f");
    m_funcs.computeGradient( n, vars, grad );
    IPOPT_DEBUG("Out MatlabProgram::eval_grad_f");
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
    IPOPT_DEBUG("In MatlabProgram::eval_g");
    if ( m > 0 ) {
      IPOPT_ASSERT(
        m_funcs.constraintFuncIsAvailable(),
        "You need to specify the callback function for computing the constraints"
      );
      m_funcs.computeConstraints( n, vars, m, g );
    }
    IPOPT_DEBUG("Out MatlabProgram::eval_g");
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

    IPOPT_DEBUG("In MatlabProgram::eval_jac_g");

    if ( numConstraints > 0 ) {
      IPOPT_ASSERT(
        m_funcs.jacobianFuncIsAvailable(),
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
          m_funcs.getJacobianNnz() == Jx_nnz,
          "The constraint Jacobian passed back from the MATLAB "
          "routine has an incorrect number of nonzero entries"
        );
        // Copy the sparse matrix structure to the IPOPT sparse matrix format.
        m_funcs.getJacobianStructure( rows, cols );
      } else {
        // Return the value of the Jacobian.
        m_funcs.computeJacobian( numConstraints, numVariables, variables, Jx );
      }
    }

    IPOPT_DEBUG("Out MatlabProgram::eval_jac_g");

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

    IPOPT_DEBUG("In MatlabProgram::eval_h");

    // If the input Hx is 0, then return the sparsity structure of the
    // Hessian. Otherwise, return the values of the nonzero entries.
    if ( Hx == nullptr ) {
      IPOPT_ASSERT(
        m_funcs.hessianFuncIsAvailable(),
        "You need to specify the callback functions for computing "
        "the Hessian and the sparsity structure of the Hessian"
      );
      // Delete any previous structure information concerning the Hessian matrix.
      // Return the sparse matrix structure of the symmetric Hessian.
      //funcs.getHessianStructure(n,H);
      IPOPT_ASSERT(
        m_funcs.getHessianNnz() == Hx_nnz,
        "The Hessian passed back from the MATLAB routine "
        "has an incorrect number of nonzero entries"
      );
      // Copy the sparse matrix structure to the IPOPT sparse matrix format.
      m_funcs.getHessianStructure( rows, cols );
    } else {
      // Return the value of the lower triangular portion of the Hessian.
      m_funcs.computeHessian( n, vars, sigma, m, lambda, Hx );
    }

    IPOPT_DEBUG("Out MatlabProgram::eval_h");

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

    IPOPT_DEBUG("In MatlabProgram::finalize_solution");

    // Store the current solution the value of the Lagrange multipliers at the solution.
    m_info.setfield( numVariables,   zl,     "zl"     );
    m_info.setfield( numVariables,   zu,     "zu"     );
    m_info.setfield( numConstraints, lambda, "lambda" );
    m_info.setfield( objective, "objective" );
    // se necessario converto in cell array
    m_funcs.fillx( numVariables, variables );
    m_info.setfield( m_funcs.mx_getx(),"x");

    IPOPT_DEBUG("Out MatlabProgram::finalize_solution");
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

    IPOPT_DEBUG("In MatlabProgram::intermediate_callback");

    if ( m_funcs.iterFuncIsAvailable() )
      return m_funcs.iterCallback(
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

  MatlabJournal::MatlabJournal( EJournalLevel default_level )
  : Journal( "MatlabJournal", default_level ) {
    IPOPT_DEBUG("MatlabJournal::MatlabJournal");
  }

  MatlabJournal::~MatlabJournal() {
    IPOPT_DEBUG("MatlabJournal::MatlabJournal");
  };

  void
  MatlabJournal::PrintImpl(
    EJournalCategory category,
    EJournalLevel    level,
    char const *     str
  ) {
    IPOPT_DEBUG("In MatlabJournal::PrintImpl");
    if ( str == nullptr ) {
      mexPrintf("MatlabJournal::PrintImpl, passed null string!!!!\n");
    } else {
      mexPrintf("%s",str);
      mexEvalString("drawnow;");
    }
    IPOPT_DEBUG("Out MatlabJournal::PrintImpl");
  }

  void
  MatlabJournal::PrintfImpl(
    EJournalCategory category,
    EJournalLevel    level,
    char const *     pformat,
    va_list          ap
  ) {

    IPOPT_DEBUG("In MatlabJournal::PrintfImpl");

    Index const maxStrLen = 4096;
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

    IPOPT_DEBUG("In MatlabJournal::PrintfImpl 2");
    mexPrintf("%s",s);
    mexEvalString("drawnow;");

    //mexEvalStringWithTrap("drawnow;"); // to dump string.

    //IPOPT_DEBUG("In MatlabJournal::PrintfImpl 3");
    //mexEvalStringWithTrap("pause(.001);"); // to dump string.

    IPOPT_ASSERT(
      nchar < maxStrLen,
      "String buffer it too short for all the characters to be printed to MATLAB console"
    );

    IPOPT_DEBUG("Out MatlabJournal::PrintfImpl");
  }
}
