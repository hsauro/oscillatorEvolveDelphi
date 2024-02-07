unit nVector_h;

// Part of Cvode.pas Declares a dynamic vector type used to pass vectors to and from
// C interfaces 

interface

type
  dArray = array[0..100000] of double;
  pdArray = ^dArray;

  TN_VectorContent_Serial = record
     length : LongInt;
     own_data : integer;
     data : pdArray;
   end;  pTN_VectorContent_Serial = ^TN_VectorContent_Serial;
  TGeneric_N_Vector = record
     content : pTN_VectorContent_Serial;
     dummy : pointer;
  end;
  pTGeneric_N_Vector = ^TGeneric_N_Vector;
  N_Vector = pTGeneric_N_Vector;

implementation

end.
