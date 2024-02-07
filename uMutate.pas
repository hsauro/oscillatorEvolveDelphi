unit uMutate;

interface

Uses Classes, SysUtils, uGlobal, uIndividual, uListOfSpecies, uListOfReactions;

function  mutateIndividual (individual : TIndividual; var experimentInfo : TExperimentInfo; var numMutations : TNumMutations) : TIndividual;
procedure deleteNthReaction (individual : TIndividual; n : integer; var numDeletions : integer);
procedure removeNthReaction (individual : TIndividual; n : integer);



implementation

Uses Math;

procedure mutateRateConstant (individual : TIndividual; var numRateConstantMutations : integer);
var n : integer; newRate : double;
begin
  if individual.listOfReactions.Count = 0 then
     exit;

  newRate := (random()-0.5)*0.2;

  n := math.RandomRange(0, individual.listOfReactions.Count);
  individual.listOfReactions[n].rateConstant := individual.listOfReactions[n].rateConstant + newRate;
  inc (numRateConstantMutations);
end;


procedure addReaction (individual : TIndividual; var experimentInfo : TExperimentInfo; var numAdditions : integer);
var r: TReaction;
    rt : double;
    s1, s2, p1, p2 : integer;
begin
   r := TReaction.Create;
   r.rateConstant := Random;
   individual.listOfReactions.Add (r);
   rt := Random;
   if rt < experimentInfo.probabilityOfUniUni then
      begin
      // UniUni
      s1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      r.reactionType := rtUniUni;
      while s1 = p1 do
            p1 := Math.RandomRange(0, individual.listOfSpecies.Count);
      r.connections[0] := s1;
      r.connections[1] := p1;

      inc (individual.listOfSpecies[s1].leavings);
      inc (individual.listOfSpecies[p1].enterings);
      end
  else if rt < experimentInfo.probabilityOfUniUni + experimentInfo.probabilityOfBiUni then
      begin
      // BiUni
      s1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      s2 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      r.reactionType := rtBiUni;
      r.connections[0] := s1;
      r.connections[1] := s2;
      r.connections[2] := p1;
      inc (individual.listOfSpecies[s1].leavings);
      inc (individual.listOfSpecies[s2].leavings);
      inc (individual.listOfSpecies[p1].enterings);
      end
   else if rt < experimentInfo.probabilityOfUniUni + experimentInfo.probabilityOfBiUni + experimentInfo.probabilityOfUniBi then
      begin
      // UniBi
      s1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p2 := Math.RandomRange(0,individual.listOfSpecies.Count);
      r.reactionType := rtUniBi;
      r.connections[0] := s1;
      r.connections[1] := p1;
      r.connections[2] := p2;
      inc (individual.listOfSpecies[s1].leavings);
      inc (individual.listOfSpecies[p1].enterings);
      inc (individual.listOfSpecies[p2].enterings);
      end
   else
      begin
      // BiBi
      s1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      s2 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p1 := Math.RandomRange(0,individual.listOfSpecies.Count);
      p2 := Math.RandomRange(0,individual.listOfSpecies.Count);
      r.reactionType := rtBiBi;
      r.connections[0] := s1;
      r.connections[1] := s2;
      r.connections[2] := p1;
      r.connections[3] := p2;
      inc (individual.listOfSpecies[s1].leavings);
      inc (individual.listOfSpecies[s2].leavings);
      inc (individual.listOfSpecies[p1].enterings);
      inc (individual.listOfSpecies[p2].enterings);
      end;
    individual.determineBoundaryNodes;
    inc  (numAdditions);
    individual.speciesId := createNewSpeciesId;
end;


procedure removeNthReaction (individual : TIndividual; n : integer);
var  s1, s2, p1, p2 : integer;
begin
  if individual.listOfReactions.Count <= 2 then
     exit;

  try
  case individual.listOfReactions[n].reactionType of
     rtUniUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtBiUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtUniBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         p2 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
     rtBiBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         p2 := individual.listOfReactions[n].connections[3];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
  end;
  individual.listOfReactions.Delete (n);
  except
    writeln ('Fatal error in removeNthReaction');
  end;
end;


procedure deleteNthReaction (individual : TIndividual; n : integer; var numDeletions : integer);
var  s1, s2, p1, p2 : integer;
begin
  if individual.listOfReactions.Count <= 2 then
     exit;

  try
  case individual.listOfReactions[n].reactionType of
     rtUniUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtBiUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtUniBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         p2 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
     rtBiBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         p2 := individual.listOfReactions[n].connections[3];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
  end;
  individual.listOfReactions.Delete (n);
  except
    writeln ('Fatal error in deleteNthReaction');
  end;
  individual.determineBoundaryNodes;
  inc (numDeletions);
  individual.speciesId := createNewSpeciesId;
end;


procedure deleteReaction (individual : TIndividual; var numDeletions : integer);
var n : integer;
    s1, s2, p1, p2 : integer;
begin
  if individual.listOfReactions.Count <= 2 then
     exit;

  // Only delete reactions that don't result in a boundary node
  n := Math.RandomRange (0, individual.listOfReactions.Count);
  try
  case individual.listOfReactions[n].reactionType of
     rtUniUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtBiUni :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         end;
     rtUniBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         p1 := individual.listOfReactions[n].connections[1];
         p2 := individual.listOfReactions[n].connections[2];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
     rtBiBi :
         begin
         s1 := individual.listOfReactions[n].connections[0];
         s2 := individual.listOfReactions[n].connections[1];
         p1 := individual.listOfReactions[n].connections[2];
         p2 := individual.listOfReactions[n].connections[3];
         dec (individual.listOfSpecies[s1].leavings);
         dec (individual.listOfSpecies[s2].leavings);
         dec (individual.listOfSpecies[p1].enterings);
         dec (individual.listOfSpecies[p2].enterings);
         end;
  end;
  individual.listOfReactions.Delete (n);
  except
    writeln ('Stage A');
  end;
  individual.determineBoundaryNodes;
  inc (numDeletions);
  individual.speciesId := createNewSpeciesId;
end;


function mutateIndividual (individual : TIndividual; var experimentInfo : TExperimentInfo; var numMutations : TNumMutations) : TIndividual;
var choice : double; i : integer;
begin
  // Clone a copy
  result := TIndividual.clone(individual);
  if result.listOfReactions.Count = 0 then
     begin
     result.badModel := True;
     exit;
     end;

  for i := 0 to 1 do
      begin
      choice := Random();
      if choice < experimentInfo.rateConstantMutation then
         begin
         mutateRateConstant (result, numMutations.numRateConstantMutations);
         continue;
         end;

      if choice < experimentInfo.rateConstantMutation + experimentInfo.addReactionMutation then
         begin
         try
           addReaction (result, experimentInfo, numMutations.numAdditions);
         except
             on E: exception do
                writeln ('Exception in add reaction mutation: ' + e.Message);
         end;
         continue;
         end;

      if choice < experimentInfo.rateConstantMutation + experimentInfo.addReactionMutation + experimentInfo.deleteReactionMutation then
         begin
         try
           deleteReaction (result, numMutations.numDeletions);
         except
             on E: exception do
                writeln ('Exception in delete reaction mutation: ' + e.Message);
         end;
         continue;
         end;
      end;
end;


end.
