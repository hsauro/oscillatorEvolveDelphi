unit uGlobal_Types;

interface

Uses System.Generics.Collections, uIndividual;

type
   TPopulation =  TObjectList<TIndividual>;

   TObjectiveFunctionData = class (TObject)
      timeStart, stepSize : double;
      numberOfPoints : Integer;
      inputData, outputData : array of double;
   end;

implementation

end.
