unit uPopulation;

interface

Uses uGlobal, uGlobal_Types, uIndividual, Generics.Defaults, SysUtils;

function createPopulation (var experimentInfo : TExperimentInfo) :  TPopulation;
procedure sortPopulation (population : TPopulation);

implementation

Uses uListOfSpecies, uListOfReactions;


function createPopulation (var experimentInfo : TExperimentInfo) : TPopulation;
var i : integer;
    t : TIndividual;
begin
  result := TPopulation.Create;
  for i := 0 to experimentInfo.maxPopulationSize - 1 do
      begin
      t := TIndividual.Create (experimentInfo);
      while t.badModel do
            begin
            t.Free;
            t := TIndividual.Create (experimentInfo);
            end;
      result.Add (t);
      end;
end;

procedure sortPopulation (population : TPopulation);
var i, j : integer;
    numberOfFloatingSpecies, numberOfOutputBoundarySpecies, numberOfInputBoundarySpecies : integer;
    fitness : double;  speciesId : integer;
    listOfSpecies : TListOfSpecies;
    listOfReactions : TListOfReactions;
begin
  for i := 0 to population.Count - 1 do
      for j := i+1 to population.Count - 1  do
          if population[i].fitness > population[j].fitness then
             begin
             try
              numberOfFloatingSpecies := population[i].  numberOfFloatingSpecies;
              numberOfOutputBoundarySpecies := population[i].numberOfOutputBoundarySpecies;
              numberOfInputBoundarySpecies := population[i].numberOfInputBoundarySpecies;
              listOfSpecies := population[i].listOfSpecies;
              listOfReactions := population[i].listOfReactions;
              fitness := population[i].fitness;
              speciesId := population[i].speciesId;

              population[i].numberOfFloatingSpecies := population[j].numberOfFloatingSpecies;
              population[i].numberOfOutputBoundarySpecies := population[j].numberOfOutputBoundarySpecies;
              population[i].numberOfFloatingSpecies := population[j].numberOfFloatingSpecies;
              population[i].listOfSpecies := population[j].listOfSpecies;
              population[i].listOfReactions := population[j].listOfReactions;
              population[i].fitness := population[j].fitness;
              population[i].speciesId := population[j].speciesId;


              population[j].numberOfFloatingSpecies := numberOfFloatingSpecies;
              population[j].numberOfOutputBoundarySpecies := numberOfOutputBoundarySpecies;
              population[j].numberOfFloatingSpecies := numberOfFloatingSpecies;
              population[j].listOfSpecies := listOfSpecies;
              population[j].listOfReactions := listOfReactions;
              population[j].fitness := fitness;
              population[j].speciesId := speciesId;

             except
               on E: exception do
                  writeln ('Error in Sort:', e.Message);
             end;
             end;
end;


initialization
end.
