unit uObjectiveFunction;

interface

Uses Classes, SysUtils, uGlobal_Types, System.StrUtils, System.Types;

function setObjectiveFunction : TObjectiveFunctionData;

implementation

function setObjectiveFunction : TObjectiveFunctionData;
var f : TextFile;
    i : integer;
    astr : string;
    alist : TStringDynArray;
begin
  result :=  TObjectiveFunctionData.Create;
  AssignFile (f, 'objectivefunction.txt');
  Reset(f);
  readln(f, astr); // Dump the first comment line
  readln(f, astr);
  aList := splitString (astr, ' ');
  result.timeStart := strtofloat (aList[1]);

  readln(f, astr);
  aList := splitString (astr, ' ');
  result.stepSize := strtofloat (aList[1]);

  readln(f, astr);
  aList := splitString (astr, ' ');
  result.numberOfPoints := strtoint (aList[1]);
  setLength (result.outputData, result.numberOfPoints);
  for i := 0 to result.numberOfPoints - 1 do
      readln (f, result.outputData[i]);
  readln (f, astr);
  closefile (f);
end;

end.
