unit uIndividual;

interface

Uses SysUtils, Classes, uListOfSpecies, uListOfReactions, Math, uGlobal;

type
   TIndividual = class (TObject)
      numberOfFloatingSpecies : integer;
      numberOfInputBoundarySpecies : Integer;
      numberOfOutputBoundarySpecies : Integer;
      listOfSpecies : TListOfSpecies;
      listOfReactions : TListOfReactions;
      fitness : double;
      badModel : boolean;
      speciesId : integer;

      constructor Create (experimentInfo : TExperimentInfo);
      destructor  Destroy; override;
      constructor clone (individual : TIndividual);
      constructor CreateTestNetwork;
      function convertToString : string;
      function getSpeciesName (index : integer) : string;
      procedure determineBoundaryNodes;
   end;

implementation

constructor TIndividual.clone (individual : TIndividual);
begin
  inherited Create;
  self.numberOfInputBoundarySpecies := individual.numberOfInputBoundarySpecies;
  self.numberOfOutputBoundarySpecies := individual.numberOfOutputBoundarySpecies;
  self.numberOfFloatingSpecies := individual.numberOfFloatingSpecies;
  self.listOfSpecies := copyListOfSpecies (individual.listOfSpecies);
  self.listOfReactions := copyListOfReactions (individual.listOfReactions);
  self.fitness := individual.fitness;
  self.badModel := individual.badModel;
  self.speciesId := individual.speciesId;
end;

function TIndividual.getSpeciesName (index : integer) : string;
begin
  if listOfSpecies[index].status = TSpeciesStatus.ssFloat then
     result := 'S' + inttostr (index)
  else
     result := '$S' + IntToStr (index);
end;


function TIndividual.convertToString : string;
var i : integer;
    S1, S2, P1, P2 : string;
begin
  result := '// Fitness of nodes:' + sLineBreak + '// ';
  for i := 0 to listOfSpecies.Count - 1 do
      result := result + format ('%8.3f, ', [listOfSpecies[i].nodefitscore]);
  result := result + sLineBreak + sLineBreak;

  for i := 0 to listOfReactions.Count - 1 do
      begin
      S1 := 'S' + inttostr (listOfReactions[i].connections[0]);
      S2 := 'S' + inttostr (listOfReactions[i].connections[1]);
      P1 := 'S' + inttostr (listOfReactions[i].connections[2]);
      P2 := 'S' + inttostr (listOfReactions[i].connections[3]);
      case listOfReactions[i].reactionType of
        rtUniUni :
           begin
           result := result + getSpeciesName (listOfReactions[i].connections[0]) + ' -> ' + getSpeciesName (listOfReactions[i].connections[1]) + ';  ' + S1 + '*' + floattostr (listOfReactions[i].rateConstant) + ';' + sLineBreak;
           end;
       rtBiUni :
           begin
           result := result + getSpeciesName (listOfReactions[i].connections[0]) + ' + ' + getSpeciesName (listOfReactions[i].connections[1]) + ' -> ' + getSpeciesName (listOfReactions[i].connections[2]) + ';  ' + S1 + '*' + S2 + '*'  + floattostr (listOfReactions[i].rateConstant) + ';' + sLineBreak;
           end;
       rtUniBi :
           begin
           result := result + getSpeciesName (listOfReactions[i].connections[0]) + ' -> ' + getSpeciesName (listOfReactions[i].connections[1]) + ' + ' + getSpeciesName (listOfReactions[i].connections[2]) + ';  ' + S1 + '*' + floattostr (listOfReactions[i].rateConstant) + ';' + sLineBreak;
           end;
       rtBiBi :
           begin
           result := result + getSpeciesName (listOfReactions[i].connections[0]) + ' + ' + getSpeciesName (listOfReactions[i].connections[1]) + ' -> ' + getSpeciesName (listOfReactions[i].connections[2]) + ' + ' + getSpeciesName (listOfReactions[i].connections[3]) + '; ' + S1 + '*' + S2 + '*' + floattostr (listOfReactions[i].rateConstant) + ';' + sLineBreak;
           end;
     end;
     end;
  result := result + sLineBreak;
  for i := 0 to listOfSpecies.Count - 1 do
      begin
      result := result + 'S' + inttostr (i) + ' = ' + floattostr (initialConcentrations[i]) + '; ';
      end;
end;

constructor TIndividual.CreateTestNetwork;
var i : Integer;
    rt : TReaction;
    sp : TSpecies;
begin
  inherited Create;

  badModel := false;
  listOfSpecies := TListOfSpecies.Create;
  listOfReactions := TListOfReactions.Create;
  for i := 0 to 7 do
      begin
      sp := TSpecies.Create;
      sp.status := TSpeciesStatus.ssFloat;
      sp.concentration := 5;
      listOfSpecies.Add(sp);
      end;
  numberOfInputBoundarySpecies := 2;
  numberOfOutputBoundarySpecies := 0;
  listOfSpecies[0].status := ssBoundaryInput;
  listOfSpecies[1].status := ssBoundaryInput;
  numberOfFloatingSpecies := 6;

  listOfSpecies[2].mappingToYVariable := 0;
  listOfSpecies[3].mappingToYVariable := 1;
  listOfSpecies[4].mappingToYVariable := 2;
  listOfSpecies[5].mappingToYVariable := 3;
  listOfSpecies[6].mappingToYVariable := 4;
  listOfSpecies[7].mappingToYVariable := 5;


  // s2 +  s4 ->  s7 +  s6; 1.8*s2*s4;
  rt := TReaction.Create;
  rt.rateConstant := 1.8;
  rt.reactionType := rtBiBi;
  rt.connections[0] := 2;  inc (listOfSpecies[2].leavings);
  rt.connections[1] := 4;  inc (listOfSpecies[4].leavings);
  rt.connections[2] := 7;  inc (listOfSpecies[7].enterings);
  rt.connections[3] := 6;  inc (listOfSpecies[6].enterings);
  listOfReactions.Add (rt);

  // s3 ->  s4; 15*s3;
  rt := TReaction.Create;
  rt.rateConstant := 15;
  rt.reactionType := rtUniUni;
  rt.connections[0] := 3;  inc (listOfSpecies[3].leavings);
  rt.connections[1] := 4;  inc (listOfSpecies[6].enterings);
  listOfReactions.Add (rt);

  // s5 +  s2 ->  s5; 4*s5*s2;
  rt := TReaction.Create;
  rt.rateConstant := 4;
  rt.reactionType := rtBiUni;
  rt.connections[0] := 5;  inc (listOfSpecies[5].leavings);
  rt.connections[1] := 2;  inc (listOfSpecies[2].leavings);
  rt.connections[2] := 5;  inc (listOfSpecies[5].enterings);
  listOfReactions.Add (rt);

  // $s1 -> s4 + s3; 2.5*s1;
  rt := TReaction.Create;
  rt.rateConstant := 2.5;
  rt.reactionType := rtUniBi;
  rt.connections[0] := 1;  inc (listOfSpecies[1].leavings);
  rt.connections[1] := 4;  inc (listOfSpecies[4].enterings);
  rt.connections[2] := 3;  inc (listOfSpecies[3].enterings);
  listOfReactions.Add (rt);

  // s5 + s6 ->  s7; 3.0*s5*s6;
  rt := TReaction.Create;
  rt.rateConstant := 3;
  rt.reactionType := rtBiUni;
  rt.connections[0] := 5;    inc (listOfSpecies[5].leavings);
  rt.connections[1] := 6;    inc (listOfSpecies[6].leavings);
  rt.connections[2] := 7;    inc (listOfSpecies[7].enterings);
  listOfReactions.Add (rt);

  // s7 +  s4 -> s2 + s3; 2*s7*s4;
  rt := TReaction.Create;
  rt.rateConstant := 2;
  rt.reactionType := rtBiBi;
  rt.connections[0] := 7;     inc (listOfSpecies[7].leavings);
  rt.connections[1] := 4;     inc (listOfSpecies[4].leavings);
  rt.connections[2] := 2;     inc (listOfSpecies[2].enterings);
  rt.connections[3] := 3;     inc (listOfSpecies[3].enterings);
  listOfReactions.Add (rt);

  // s2 ->  s2 +  s5; 12*s2;
  rt := TReaction.Create;
  rt.rateConstant := 12;
  rt.reactionType := rtUniBi;
  rt.connections[0] := 2;     inc (listOfSpecies[2].leavings);
  rt.connections[1] := 2;     inc (listOfSpecies[2].enterings);
  rt.connections[2] := 5;     inc (listOfSpecies[5].enterings);
  listOfReactions.Add (rt);

  // s2 +  s3 ->  s7; 2*s2*s3;
  rt := TReaction.Create;
  rt.rateConstant := 2;
  rt.reactionType := rtBiUni;
  rt.connections[0] := 2;      inc (listOfSpecies[2].leavings);
  rt.connections[1] := 3;      inc (listOfSpecies[3].leavings);
  rt.connections[2] := 7;      inc (listOfSpecies[7].enterings);
  listOfReactions.Add (rt);

  // $s0 ->  s3 +  s5; 2*s0;
  rt := TReaction.Create;
  rt.rateConstant := 2;
  rt.reactionType := rtUniBi;
  rt.connections[0] := 0;      inc (listOfSpecies[0].leavings);
  rt.connections[1] := 3;      inc (listOfSpecies[3].enterings);
  rt.connections[2] := 5;      inc (listOfSpecies[5].enterings);
  listOfReactions.Add (rt);

  // s6 ->  s4 +  s6; s6;
  rt := TReaction.Create;
  rt.rateConstant := 1;
  rt.reactionType := rtUniBi;
  rt.connections[0] := 6;      inc (listOfSpecies[6].leavings);
  rt.connections[1] := 4;      inc (listOfSpecies[4].enterings);
  rt.connections[2] := 6;      inc (listOfSpecies[6].enterings);
  listOfReactions.Add (rt);
end;


constructor TIndividual.Create (experimentInfo : TExperimentInfo);
var nSpecies, nReactions, i : Integer; rt : double;
    s1, s2, p1, p2 : integer;
    r : TReaction;
begin
  inherited Create;

  speciesId := createNewSpeciesId();

  badModel := false;
  listOfSpecies := TListOfSpecies.Create;
  listOfReactions := TListOfReactions.Create;
  nSpecies := Math.RandomRange (experimentInfo.minNumberOfSpecies, experimentInfo.maxNumberOfSpecies);
  for i := 0 to nSpecies - 1 do
      listOfSpecies.Add (TSpecies.Create);
  numberOfFloatingSpecies := nSpecies; // Reduced when we discover boundary nodes.
  nReactions := Math.RandomRange (experimentInfo.minNumberOfReactions, experimentInfo.maxNumberOfReactions);
  for i := 0 to nReactions - 1 do
      begin
      r := TReaction.Create;
      r.rateConstant := Random*experimentInfo.maxInitialRateConstant;
      listOfReactions.Add (r);
      rt := Random();
      if rt < experimentInfo.probabilityOfUniUni then
         begin
         // UniUni
         s1 := Math.RandomRange(0,listOfSpecies.Count);
         p1 := Math.RandomRange(0,listOfSpecies.Count);
         while s1 = p1 do
               p1 := Math.RandomRange(0,listOfSpecies.Count);

         r.reactionType := rtUniUni;
         r.connections[0] := s1; inc (listOfSpecies[s1].leavings);
         r.connections[1] := p1; inc (listOfSpecies[p1].enterings);
         end
      else if rt < experimentInfo.probabilityOfUniUni + experimentInfo.probabilityOfBiUni then
         begin
         // BiUni
         s1 := Math.RandomRange(0,listOfSpecies.Count);
         s2 := Math.RandomRange(0,listOfSpecies.Count);
         p1 := Math.RandomRange(0,listOfSpecies.Count);
         r.reactionType := rtBiUni;
         r.connections[0] := s1; inc (listOfSpecies[s1].leavings);
         r.connections[1] := s2; inc (listOfSpecies[s2].leavings);
         r.connections[2] := p1; inc (listOfSpecies[p1].enterings);
         end
      else if rt < experimentInfo.probabilityOfUniUni + experimentInfo.probabilityOfBiUni + experimentInfo.probabilityOfUniBi then
         begin
         // UniBi
         s1 := Math.RandomRange(0,listOfSpecies.Count);
         p1 := Math.RandomRange(0,listOfSpecies.Count);
         p2 := Math.RandomRange(0,listOfSpecies.Count);
         r.reactionType := rtUniBi;
         r.connections[0] := s1; inc (listOfSpecies[s1].leavings);
         r.connections[1] := p1; inc (listOfSpecies[p1].enterings);
         r.connections[2] := p2; inc (listOfSpecies[p2].enterings);
         end
      else
         begin
         // BiBi
         s1 := Math.RandomRange(0,listOfSpecies.Count);
         s2 := Math.RandomRange(0,listOfSpecies.Count);
         p1 := Math.RandomRange(0,listOfSpecies.Count);
         p2 := Math.RandomRange(0,listOfSpecies.Count);
         r.reactionType := rtBiBi;
         r.connections[0] := s1; inc (listOfSpecies[s1].leavings);
         r.connections[1] := s2; inc (listOfSpecies[s2].leavings);
         r.connections[2] := p1; inc (listOfSpecies[p1].enterings);
         r.connections[3] := p2; inc (listOfSpecies[p2].enterings);
         end;
     end;
  determineBoundaryNodes;
  if numberOfFloatingSpecies = 0 then
     badModel := true;
  if numberOfInputBoundarySpecies = 0 then
     badModel := true;
end;


destructor TIndividual.Destroy;
begin
  listOfSpecies.Free;
  listOfReactions.Free;
  inherited;
end;


procedure TIndividual.determineBoundaryNodes;
var i : Integer;
    yCount : integer;
begin
  // Determine floating or boundary nodes
  numberOfOutputBoundarySpecies := 0;
  numberOfInputBoundarySpecies := 0;
  // Assume all species are initially float.
  numberOfFloatingSpecies := listOfSpecies.Count;
  yCount := 0;
  for i := 0 to listOfSpecies.Count - 1 do
      begin
      listOfSpecies[i].status := TSpeciesStatus.ssFloat;
      if (listOfSpecies[i].leavings = 0) and (listOfSpecies[i].enterings = 0) then
          begin
          listOfSpecies[i].status := TSpeciesStatus.ssOrphan;
          Continue;
          end;

      if listOfSpecies[i].leavings = 0 then
         begin
         listOfSpecies[i].status := TSpeciesStatus.ssBoundaryOutput;
         listOfSpecies[i].concentration := 0;
         inc (numberOfOutputBoundarySpecies);
         dec (numberOfFloatingSpecies);
         Continue;
         end;
      if listOfSpecies[i].enterings = 0 then
         begin
         listOfSpecies[i].status := TSpeciesStatus.ssBoundaryInput;
         listOfSpecies[i].concentration := initialConcentrations[i];
         inc (numberOfInputBoundarySpecies);
         dec (numberOfFloatingSpecies);
         Continue;
         end;
      listOfSpecies[i].mappingToYVariable := yCount;
      Inc(yCount);
      end;

end;


end.
