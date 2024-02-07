unit uCvode;


// Dephi interface to CVODE DLL  (CVODE Version 2.6.2 - Sept 2016)

interface

uses Windows, SysUtils, nVector_h;


const
  // lmm
  CV_ADAMS = 1;
  CV_BDF =  2;

  // iter
  CV_FUNCTIONAL = 1;
  CV_NEWTON  =   2;

  CV_NORMAL = 1;
  CV_ONE_STEP = 2;

  CV_DEFAULT_RTOL   = 1E-6;
  CV_DEFAULT_ATOL   = 1E-12;
  CV_DEFAULT_MAXNUMSTEPS = 8000;
  CV_DEFAULT_BDFOrder    = 5;
  CV_DEFAULT_AdamsOrder  = 12;

  // cvode main solver module  (including CVodeinit)
  CV_SUCCESS   =  0;
  CV_TSTOP_RETURN =  1;  // CVode succeeded by reaching the specified stopping point
  CV_ROOT_RETURN  =  2;  // CVode succeeded and found one or more roots.
  CV_WARNING      =  99; // CVode succeeded but an unusual situation occurred.

  CV_TOO_MUCH_WORK= -1;    // The solver took mxstep internal steps but could not reach tout.
  CV_TOO_MUCH_ACC = -2;    // The solver could not satisfy the accuracy demanded by the user for some internal step
  CV_ERR_FAILURE  = -3;    // Error test failures occurred too many times during one internal time step or minimum step size was reached.
  CV_CONV_FAILURE = -4;    // Convergence test failures occurred too many times during one internal time step or minimum step size was reached
  CV_LINIT_FAIL =   -5;    // The linear solver’s initialization function failed.
  CV_LSETUP_FAIL =  -6;    // The linear solver’s setup function failed in an unrecoverable manner.
  CV_LSOLVE_FAIL =  -7;    // The linear solver’s solve function failed in an unrecoverable manner.
  CV_RHSFUNC_FAIL = -8;    // The right-hand side function failed in an unrecoverable manner.
  CV_FIRST_RHSFUNC_ERR =   -9;  // The right-hand side function failed at the first call
  CV_REPTD_RHSFUNC_ERR =   -10; // The right-hand side function had repetead recoverable errors
  CV_UNREC_RHSFUNC_ERR =   -11; // The right-hand side function had a recoverable error, but no recovery is possible.
   CV_RTFUNC_FAIL      =   -12; // The rootfinding function failed in an unrecoverable manner.

  CV_MEM_FAIL  = -20;  // A memory allocation failed.
  CV_MEM_NULL  = -21;  // The cvode mem argument was NULL.
  CV_ILL_INPUT = -22;  // One of the function inputs is illegal
  CV_NO_MALLOC = -23;  // The cvode memory block was not allocated by a call to CVodeMalloc.
  CV_BAD_K =     -24;  // The derivative order k is larger than the order used
  CV_BAD_T =     -25;  // The time t is outside the last step taken.
  CV_BAD_DKY =   -26;  // The output derivative vector is NULL.
  CV_TOO_CLOSE = -27;  // The output and initial times are too close to each other.

  // cvdls linear solver modules
  CVDLS_SUCCESS = 0;
  CVDLS_MEM_NULL = -1;
  CVDLS_LMEM_NULL = -2;
  CVDLS_ILL_INPUT = -3;
  CVDLS_MEM_FAIL = -4;
  CVDLS_JACFUNC_UNRECVR = -5;
  CVDLS_JACFUNC_RECVR = -6;


type
  // These will be used to declare dynamic arrays
  TDArray = array of double;
  TCallBackFcn = function (time : double; y, ydot : N_Vector; useData : pointer) : integer; cdecl;
  TModelfcn = procedure (time : Double; y, dydt : TDArray) of object;

  TUserData = class (TObject)
     n : integer;
     Model : TModelFcn;
     yv, ydotv : TDArray;
     constructor Create (n : integer; model : TModelFcn);
     destructor  destroy; override;
  end;

  TCvode = class (TObject)
     private
      nSize : integer;
      cv_y : N_Vector;
      cvodeMem : Pointer;
      Active : Boolean;
      timeReached : Double;
      errorFlag : integer;

      userData : TUserData;
      procedure   checkForError (ierr : integer);
      procedure   freeMemory;
     public
      constructor Create (n : integer; model : TModelFcn; eps : Double); overload;
      constructor Create (n : Integer); overload;
      destructor  Destroy; override;
      procedure   reset;
      procedure   initialise (y : TDArray; timeStart : double);
      procedure   setModel (model : TModelfcn);
      function    execute (y : TDArray; StartTime, HStep : Double) : Double;

      procedure   setRTol (i : integer; Value : Double);
      function    getRtol (i : integer) : Double;

      procedure   setATol (i : integer; Value : Double);
      function    getAtol (i : integer) : Double;

      procedure   setMaxSteps (Value : integer);
      function    getMaxSteps : integer;

      procedure   setMaxAdamsOrder (value : integer);
      function    getMaxAdamsOrder : integer;

      procedure   setMaxBDFOrder (value : integer);
      function    getMaxBDFOrder : integer;

      property    rtol[Index: Integer] : double read GetRtol write SetRtol;
      property    atol[Index: Integer] : double read GetAtol write SetAtol;
      property    maxSteps : integer read getMaxSteps write setMaxSteps;
      property    AdamsOrder : integer read getMaxAdamsOrder write setMaxAdamsOrder;
      property    BDFOrder : integer read getMaxBDFOrder write setMaxBDFOrder;
  end;

function LoadDLL : Boolean;
procedure Release;

var
    persistent_cvode_relativeTolerance : double = CV_DEFAULT_RTOL;
    persistent_cvode_absoluteTolerance : double = CV_DEFAULT_ATOL;
    persistent_cvode_maximumSteps      : integer = CV_DEFAULT_MAXNUMSTEPS;
    persistent_cvode_BDFOrder          : integer = CV_DEFAULT_BDFOrder;
    persistent_cvode_maxAdamsOrder     : integer = CV_DEFAULT_AdamsOrder;

implementation

type
  TCVodeErrorFunction = procedure (errorCode : integer; module, functionName, msg : PAnsiChar; errorData : pointer); cdecl;

  TNewCvode_Vector = function (n : integer) : N_Vector; cdecl;
  TFreeCvode_Vector = procedure (v : N_Vector); cdecl;
  TCVode_GetArrayPointer = function (v : N_Vector):  pDArray; cdecl;
  TCvode_SetVector = procedure (abstol : N_Vector; Index : integer; Value : double); cdecl;
  TCvode_GetVector = function (v : N_Vector; Index : integer) : double; cdecl;

  TCVodeSetMaxOrder = function (p : pointer; maxOrder : integer) : integer; cdecl;
  TSetMaxErrTestFails = function (p : pointer; maxnef : integer) : integer; cdecl;
  TSetMaxConvFails = function (p : pointer; maxconv : integer) : integer; cdecl;

  TCVodeCreate = function (lmm, iter : integer) : pointer; cdecl;
  TCVodeFreeMem = procedure (var p : Pointer); cdecl;
  TCVodeInit = function (cvodeMem : pointer; f : TCallBackFcn; initialTime : double; initialY : N_Vector) : integer; cdecl;
  TCVodeReInit = function (cvodeMem : pointer; initialTime : double; initialY : pTGeneric_N_Vector) : integer; cdecl;
  TCvodeSetUserData = function (cvodeMem : pointer; userData : pointer) : integer; cdecl;
  TCvDense = function (p : Pointer; size : integer) : integer; cdecl;
  TCvodeRun = function (cvodeMem : pointer; tout : double; y : N_Vector; var t : double; itask : integer) : integer; cdecl;
  TCVodeSStolerances = function (cvodeMem : pointer; reltol, abstol : double) : integer; cdecl;
  TCVodeSetMaxNumSteps = function (cvodeMem : Pointer; mxsteps : integer) : integer; cdecl;
  TCVodeSetErrHandlerF = function (cvodeMem : Pointer; errorFunction : TCVodeErrorFunction; errorData : pointer) : integer; cdecl;

var
    Handle : THandle;

    newCvode_Vector : TNewCvode_Vector = nil;
    freeCvode_Vector : TFreeCvode_Vector = nil;
    CVode_GetArrayPointer : TCVode_GetArrayPointer = nil;

    CVodeCreate : TCVodeCreate = nil;
    CVodeFreeMem : TCVodeFreeMem = nil;
    CVodeInit : TCVodeInit = nil;
    CVodeReInit : TCVodeReInit = nil;
    CVodeSetUserData : TCvodeSetUserData = nil;
    CvDense : TCvDense = nil;
    CVodeRun : TCVodeRun;
    CVodeSStolerances : TCVodeSStolerances;
    CVodeSetMaxNumSteps : TCVodeSetMaxNumSteps = nil;

    CVodeSetMaxOrder : TCVodeSetMaxOrder = nil;
    SetMaxErrTestFails : TSetMaxErrTestFails = nil;
    SetMaxConvFails : TSetMaxConvFails = nil;
    CVodeSetErrHandlerF : TCVodeSetErrHandlerF = nil;


function loadNVectorDLL : Boolean;
var libPath : string;
    ret : cardinal;
    err :string;

begin
  libPath := ExtractFilePath(ParamStr(0)) + 'nvecserial.dll';
  if not FileExists (libPath) then
     raise Exception.Create ('Unable to locate nvecserial.dll at:' + sLineBreak + sLineBreak + libPath);
  Handle := LoadLibrary(PWideChar (libPath));

  if Handle = 0 then
     begin
     ret := GetLastError();
     err := SysErrorMessage(ret);
     Writeln(err);
     end;

  if Handle >= 32 then
     begin
     @newCvode_Vector := GetProcAddress (Handle, 'N_VNew_Serial');
     @freeCvode_Vector := GetProcAddress (Handle, 'N_VDestroy_Serial');
     @CVode_GetArrayPointer := GetProcAddress (Handle, 'N_VGetArrayPointer_Serial')
     end
  else
     exit (false);

  result := true;
end;


function LoadDLL : Boolean;
var libPath : string;
begin
  if not loadNVectorDLL then
     raise Exception.Create ('Failed to correctly load nvecserial dll');

  libPath := ExtractFilePath(ParamStr(0)) + 'cvode_262.dll';
  if not FileExists (libPath) then
     raise Exception.Create ('Unable to locate cvode_262.dll at:' + sLineBreak + sLineBreak + libPath);

  Handle := LoadLibrary(PWideChar (libPath));
  if Handle >= 32 then
     begin
     @CVodeCreate := GetProcAddress (Handle, 'CVodeCreate'); if not Assigned (CVodeCreate) then exit (false);
     @CvodeInit := GetProcAddress (Handle, 'CVodeInit');     if not Assigned (CvodeInit) then exit (false);
     @CVodeReInit := GetProcAddress (Handle, 'CVodeReInit'); if not Assigned (CVodeReInit) then exit (false);
     @CVodeSetUserData := GetProcAddress (Handle, 'CVodeSetUserData'); if not Assigned (CVodeSetUserData) then exit (false);
     @CvDense := GetProcAddress (Handle, 'CVDense'); if not Assigned (CvDense) then exit (false);
     @CVodeRun := GetProcAddress (Handle, 'CVode');  if not Assigned (CVodeRun) then exit (false);
     @CVodeSStolerances := GetProcAddress (Handle, 'CVodeSStolerances');     if not Assigned (CVodeSStolerances) then exit (false);
     @CVodeFreeMem := GetProcAddress (Handle, 'CVodeFree');                  if not Assigned (CVodeFreeMem) then exit (false);
     @CVodeSetMaxNumSteps := GetProcAddress (Handle, 'CVodeSetMaxNumSteps'); if not Assigned (CVodeSetMaxNumSteps) then exit (false);
     @CVodeSetErrHandlerF := GetProcAddress (Handle, 'CVodeSetErrHandlerFn');

//     @CVodeSetMaxOrder := GetProcAddress (Handle, 'SetMaxOrder');
//     @SetMaxErrTestFails := GetProcAddress (Handle, 'SetMaxErrTestFails');
//     @SetMaxConvFails := GetProcAddress (Handle, 'SetMaxConvFails');
     end
  else
     begin
     result := false;
     exit;
     end;

  result := true;
end;


procedure Release;
begin
  FreeLibrary(Handle);
end;

// -----------------------------------------------------------------------

constructor TUserData.Create (n : integer; model : TModelFcn);
begin
  inherited Create;
  self.n := n;
  Self.model := model;
  setLength(yv, n);
  setLength(ydotv, n);
end;


destructor TUserData.destroy;
begin
  setLength(yv, 0);
  setLength(ydotv, 0);
  inherited Destroy;
end;


// -----------------------------------------------------------------------


procedure setNVectorValue (v : pTGeneric_N_Vector; index : Integer; value : double);
var d : pdArray;
begin
  d := CVode_GetArrayPointer (v);
  d[index] := value;
end;

function getNVectorValue (v : pTGeneric_N_Vector; index : Integer) : double;
var d : pdArray;
begin
  d := CVode_GetArrayPointer (v);
  result := d[index];
end;

// -----------------------------------------------------------------------


// This is called from the Cvode DLL to get the dy/dt values
function callBackFcn (Time : Double; y, ydot : N_Vector; userData : pointer) : integer cdecl;
var i : integer;
    fData : TUserData;
    model : TModelfcn;
    yx : double;
begin
  fData := userData;

  for i := 0 to fData.n - 1 do
      fdata.yv[i] := getNVectorValue(y, i);

  model := TModelFcn (fData.Model);
  model (Time, fdata.yv, fdata.ydotv);
  for i := 0 to  fData.n - 1 do
      setNVectorValue(ydot, i, fdata.ydotv[i]);

  result := 0; // Important!
end;

// -----------------------------------------------------------------------


procedure cvodeErrorFunction (errorCode : integer; module, functionName, msg : PAnsiChar; errorData : pointer); cdecl;
var strModule, strFunctionName, strMsg : string;
begin
  strModule := string (module);
  strFunctionName := string (functionName);
  strMsg := string (msg);
  //writeln ('[' + strModule + ': ' + strFunctionName + ']' + strMsg);
  //raise Exception.Create ('[' + strModule + ': ' + strFunctionName + ']' + strMsg);
end;


constructor TCvode.Create (n : Integer);
var i : integer;
    ierr : integer;
    y : array of double;
    errorData : pointer;
begin
  nSize := n;
  Active := True;
  if @CvodeRun = nil then
     begin
     Active := False;
     Exit;
     end;

  try
    cv_y := NewCvode_Vector (n);
    for i := 0 to length (y) - 1 do
        setNVectorValue (cv_y, i, 0.0);

    userData := TUserData.Create (n, nil);
    cvodeMem := CVodeCreate (CV_BDF, CV_NEWTON);
    if cvodeMem = nil then
       raise Exception.Create ('Nil returned from CVodeCreate');

    ierr := CVodeSetErrHandlerF (cvodeMem, cvodeErrorFunction, errorData);
    checkForError (ierr);
    ierr := CVodeInit (cvodeMem, callBackFcn, 0.0, cv_y);
    checkForError (ierr);

    ierr := CVodeSetUserData (cvodeMem, userData);
    checkForError (ierr);

  except
    on E:exception do
       writeln ('Error in Execute: ' + e.message);
  end;

  ierr := CVodeSStolerances (cvodeMem, persistent_cvode_relativeTolerance, persistent_cvode_absoluteTolerance);
  checkForError (ierr);

  //CVodeSetMaxNumSteps (Cvode_mem, persistent_cvode_maximumSteps);

  ierr := CVDense (cvodeMem, nSize);
  case ierr of
      CVDLS_SUCCESS   : exit;
      CVDLS_MEM_NULL  : raise Exception.Create ('CVODE DENSE: Internal error, memory not allocated correctly for CVODE');
      CVDLS_ILL_INPUT : raise Exception.Create ('CVODE DENSE: Internal error, illegal input to CVDense');
      CVDLS_MEM_FAIL  : raise Exception.Create ('CVODE DENSE: Internal error, unable to allocate sufficient memory');
      CVDLS_JACFUNC_UNRECVR : raise Exception.Create ('The Jacobian function failed in an unrecoverable manner');
      CVDLS_JACFUNC_RECVR  : raise Exception.Create ('The Jacobian function had a recoverable error');
  else
       raise Exception.Create ('CVODE DENSE: Internal error, unknown error reported in CVDense: ' + inttostr (ierr));
  end;
end;


constructor TCvode.Create (n : Integer; model : TModelFcn; eps : Double);
var i : integer;
    ierr : integer;
begin
  nSize := n;
  Active := True;
  if @CvodeRun = nil then
     begin
     active := False;
     Exit;
     end;

  cv_y := NewCvode_Vector (n);

  userData := TUserData.Create (n, model);
  cvodeMem := CVodeCreate (CV_BDF, CV_NEWTON);
  ierr := CVodeSetUserData (cvodeMem, userData);
  checkForError (ierr);

  ierr := CVodeInit (cvodeMem, callBackFcn, 0.0, cv_y);
  checkForError (ierr);

  ierr := CVodeSStolerances (cvodeMem, persistent_cvode_relativeTolerance, persistent_cvode_absoluteTolerance);
  checkForError (ierr);

  //CVodeSetMaxNumSteps (Cvode_mem, persistent_cvode_maximumSteps);

  ierr := CVDense (cvodeMem, nSize);
  case ierr of
      CVDLS_SUCCESS   : exit;
      CVDLS_MEM_NULL  : raise Exception.Create ('CVODE DENSE: Internal error, memory not allocated correctly for CVODE');
      CVDLS_ILL_INPUT : raise Exception.Create ('CVODE DENSE: Internal error, illegal input to CVDense');
      CVDLS_MEM_FAIL  : raise Exception.Create ('CVODE DENSE: Internal error, unable to allocate sufficient memory');
      CVDLS_JACFUNC_UNRECVR : raise Exception.Create ('The Jacobian function failed in an unrecoverable manner');
      CVDLS_JACFUNC_RECVR  : raise Exception.Create ('The Jacobian function had a recoverable error');
  else
       raise Exception.Create ('CVODE DENSE: Internal error, unknown error reported in CVDense: ' + inttostr (ierr));
  end;
end;


destructor TCvode.Destroy;
begin
  freeMemory;
  inherited destroy;
end;


// Checks for main cvode module error codes
procedure TCvode.checkForError (ierr : integer);
begin
  if ierr < 0 then
     begin
     case ierr of
        CV_MEM_NULL      : raise Exception.Create ('CVODE: The cvode memory block was not initialized through a previous call to CVodeCreate');
        CV_MEM_FAIL      : raise Exception.Create ('CVODE: A memory allocation request has failed');
        CV_ILL_INPUT     : raise Exception.Create ('CVODE: Illegal input to a CVODE function');

        CV_NO_MALLOC     : raise Exception.Create ('CVODE: Memory not correctly allocated during initialization');
        CV_TOO_MUCH_WORK : raise Exception.Create ('CVODE: Too many internal steps required to reach end time point');
        CV_TOO_MUCH_ACC  : raise Exception.Create ('CVODE: Too much accuracy demanded by user, reduce tolerance');
        CV_ERR_FAILURE   : raise Exception.Create ('CVODE: Too many internal test errors occured during integration');
        CV_CONV_FAILURE  : raise Exception.Create ('CVODE: The linear solver''s initialization routine failed in an unrecoverable manner');
        CV_LSETUP_FAIL   : raise Exception.Create ('CVODE: The linear solver''s setup routine failed in an unrecoverable manner');
        CV_LSOLVE_FAIL   : raise Exception.Create ('CVODE: The linear solver''s solver routine failed in an unrecoverable manner');
        CV_RHSFUNC_FAIL  : raise Exception.Create ('The right-hand side function failed in an unrecoverable manner');
        CV_FIRST_RHSFUNC_ERR : raise Exception.Create ('The right-hand side function failed at the first call.');
        CV_REPTD_RHSFUNC_ERR : raise Exception.Create ('The right-hand side function had repetead recoverable errors.');
        CV_UNREC_RHSFUNC_ERR : raise Exception.Create (' The right-hand side function had a recoverable error, but no recovery is possible');
        CV_RTFUNC_FAIL : raise Exception.Create ('The rootfinding function failed in an unrecoverable manner');
        CV_BAD_K : raise Exception.Create ('The derivative order k is larger than the order used.');
        CV_BAD_T : raise Exception.Create ('The time t is outside the last step taken');
        CV_BAD_DKY : raise Exception.Create ('The output derivative vector is NULL.');
        CV_TOO_CLOSE : raise Exception.Create ('The output and initial times are too close to each other');
      else
       raise Exception.Create ('An unknown error has occured in CVODE: ' + inttostr (ierr));
     end;
     end;
end;


function TCvode.execute (y : TDArray; StartTime, HStep : Double) : Double;
var timeOut : Double; i, ierr : integer;
begin
  timeOut := startTime + Hstep;
  if length (y) = 0 then
     exit;

  for i := 0 to length (y) - 1 do
      setNVectorValue (cv_y, i, y[i]);

  ierr := CVodeReInit (cvodeMem, startTime, cv_y);
  checkForError (ierr);

  ierr := CvodeRun (cvodeMem, timeOut, cv_y, timeReached, CV_NORMAL);
  checkForError (ierr);

  result := timeReached;
  for i := 0 to length (y) - 1 do
      y[i] := getNVectorValue (cv_y, i);
end;


procedure TCvode.Initialise (y : TDArray; timeStart : double);
var i : integer;
    ierr : integer;
begin
  for i := 0 to length (y) - 1 do
      setNVectorValue (cv_y, i, y[i]);

  ierr := CVodeReInit (cvodeMem, timeStart, cv_y);
  checkForError (ierr);
end;


procedure TCVode.freeMemory;
begin
  try
    if cv_y <> nil then FreeCvode_Vector (cv_y);
    if CvodeMem <> nil then CVodeFreeMem (CvodeMem);
    userData.Free;
  finally
    cv_y := nil;
    CvodeMem := nil;
  end;
end;


procedure TCvode.setModel (model : TModelfcn);
begin
  userData.Model := model;
end;


procedure TCvode.Reset;
var errCode : integer;
begin
  errCode := CVodeReInit (cvodeMem, timeReached, cv_y);
  checkForError (errCode);
end;


procedure TCvode.SetRTol (i : integer; Value : Double);
begin
  persistent_cvode_relativeTolerance := value;
end;


function TCvode.GetRtol (i : integer) : Double;
begin
  Result := persistent_cvode_absoluteTolerance;
end;


procedure TCvode.SetATol (i : integer; Value : Double);
begin
  //setNVectorValue (abstol, i, value);
  persistent_cvode_absoluteTolerance := value;
end;


function TCvode.GetAtol (i : integer) : Double;
begin
  //result := getNVectorValue (abstol, i);
  result := persistent_cvode_absoluteTolerance;
end;


procedure TCvode.SetMaxSteps (Value : integer);
begin
  CVodeSetMaxNumSteps (CvodeMem, value);
  persistent_cvode_maximumSteps := value;
end;


function TCvode.GetMaxSteps : integer;
begin
  result := persistent_cvode_maximumSteps;
end;


procedure TCvode.setMaxAdamsOrder (value : integer);
var flag : integer;
begin
  flag := CVodeSetMaxOrder(CvodeMem, value);
  persistent_cvode_maxAdamsOrder := value;
end;


function TCvode.getMaxAdamsOrder : integer;
begin
  result := persistent_cvode_maxAdamsOrder;
end;


procedure TCvode.setMaxBDFOrder (value : integer);
var flag : integer;
begin
  flag := CVodeSetMaxOrder (CvodeMem, value);
  persistent_cvode_BDFOrder := value;
end;


function TCvode.getMaxBDFOrder : integer;
begin
  result := persistent_cvode_BDFOrder;
end;


end.


