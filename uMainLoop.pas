unit uMainLoop;

interface

Uses Classes, SysUtils, uGlobal, uGlobal_Types, ulibRK4, uIndividual, System.SyncObjs;

function startMainLoop (var population : TPopulation; objFunction : TObjectiveFunctionData; var experimentInfo : TExperimentInfo) : integer;
function simulate (individual : TIndividual; timeStart, stepSize : double; numberOfPoints : Integer) : T2DArray;
function computeFitnessOfIndividual (individual : TIndividual; objFunction : TObjectiveFunctionData) : double;
procedure computeFitnessOfPopulation (population : TPopulation; objFunction : TObjectiveFunctionData);


implementation

Uses uListOfReactions, uListOfSpecies, uCvode, uPopulation, Math, uMutate;

type
   TModelClass = class (TObject)
      class var individual : TIndividual;
      class procedure myfunc (time : double; y, dydt : TDArray);
   end;


class procedure TModelClass.myfunc (time : double; y, dydt : TDArray);
var i, j : integer;
    yCount : integer;
    s1, s2, p1, p2 : integer;
    concS1, concS2 : double;
    reactionRate : double;
begin
//  dydt[0] := -0.4*y[0];
//  dydt[1] := -0.4*y[1];
//  dydt[2] := -0.4*y[2];
//  exit;

  for i := 0 to individual.numberOfFloatingSpecies - 1 do dydt[i] := 0;

  for j := 0 to individual.listOfReactions.Count - 1 do
      case individual.listOfReactions[j].reactionType of
                 rtUniUni : begin
                            // Compute: k*s1, if s1 is boundary use listofspecies, else use y
                            s1 := individual.listOfReactions[j].connections[0];
                            p1 := individual.listOfReactions[j].connections[1];
                            // Get hold of the reactant concentation.
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               concS1 := y[individual.listOfSpecies[s1].mappingToYVariable]
                            else
                               concS1 := individual.listOfSpecies[individual.listOfReactions[j].connections[0]].concentration;

                            // Compute the reaction rate
                            reactionRate := individual.listOfReactions[j].rateConstant*concS1;

                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s1].mappingToYVariable] := dydt[individual.listOfSpecies[s1].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[p1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p1].mappingToYVariable] := dydt[individual.listOfSpecies[p1].mappingToYVariable] + reactionRate
                            end;
                  rtBiUni : begin
                            // Compute: k*s1*s2, if s1 is boundary use listofspecies, else use y
                            s1 := individual.listOfReactions[j].connections[0];
                            s2 := individual.listOfReactions[j].connections[1];
                            p1 := individual.listOfReactions[j].connections[2];

                            // Get hold of the two reactant concentations.
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               concS1 := y[individual.listOfSpecies[s1].mappingToYVariable]
                            else
                               concS1 := individual.listOfSpecies[individual.listOfReactions[j].connections[0]].concentration;

                            if individual.listOfSpecies[s2].status = TSpeciesStatus.ssFloat then
                               concS2 := y[individual.listOfSpecies[s2].mappingToYVariable]
                            else
                               concS2 := individual.listOfSpecies[individual.listOfReactions[j].connections[1]].concentration;

                            // Compute the reaction rate
                            reactionRate := individual.listOfReactions[j].rateConstant*concS1*concS2;

                            // Compute the dydys
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s1].mappingToYVariable] := dydt[individual.listOfSpecies[s1].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[s2].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s2].mappingToYVariable] := dydt[individual.listOfSpecies[s2].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[p1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p1].mappingToYVariable] := dydt[individual.listOfSpecies[p1].mappingToYVariable] + reactionRate;
                            end;
                  rtUniBi : begin
                            // Compute: k*s1*s2, if s1 is boundary use listofspecies, else use y
                            s1 := individual.listOfReactions[j].connections[0];
                            p1 := individual.listOfReactions[j].connections[1];
                            p2 := individual.listOfReactions[j].connections[2];

                            // Get hold of the reactant concentations
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               concS1 := y[individual.listOfSpecies[s1].mappingToYVariable]
                            else
                               concS1 := individual.listOfSpecies[individual.listOfReactions[j].connections[0]].concentration;

                            // Compute the reaction rate
                            reactionRate := individual.listOfReactions[j].rateConstant*concS1;

                            // Compute the dydys
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s1].mappingToYVariable] := dydt[individual.listOfSpecies[s1].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[p1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p1].mappingToYVariable] := dydt[individual.listOfSpecies[p1].mappingToYVariable] + reactionRate;

                            if individual.listOfSpecies[p2].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p2].mappingToYVariable] := dydt[individual.listOfSpecies[p2].mappingToYVariable] + reactionRate;
                            end;
                  rtBiBi : begin
                            // Compute: k*s1*s2, if s1 is boundary use listofspecies, else use y
                            s1 := individual.listOfReactions[j].connections[0];
                            s2 := individual.listOfReactions[j].connections[1];
                            p1 := individual.listOfReactions[j].connections[2];
                            p2 := individual.listOfReactions[j].connections[3];

                            // Get hold of the reactant concentations
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               concS1 := y[individual.listOfSpecies[s1].mappingToYVariable]
                            else
                               concS1 := individual.listOfSpecies[individual.listOfReactions[j].connections[0]].concentration;

                            // Get hold of the reactant concentations
                            if individual.listOfSpecies[s2].status = TSpeciesStatus.ssFloat then
                               concS2 := y[individual.listOfSpecies[s2].mappingToYVariable]
                            else
                               concS2 := individual.listOfSpecies[individual.listOfReactions[j].connections[1]].concentration;

                            // Compute the reaction rate
                            reactionRate := individual.listOfReactions[j].rateConstant*concS1*concS2;

                            // Compute the dydys
                            if individual.listOfSpecies[s1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s1].mappingToYVariable] := dydt[individual.listOfSpecies[s1].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[s2].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[s2].mappingToYVariable] := dydt[individual.listOfSpecies[s2].mappingToYVariable] - reactionRate;

                            if individual.listOfSpecies[p1].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p1].mappingToYVariable] := dydt[individual.listOfSpecies[p1].mappingToYVariable] + reactionRate;

                            if individual.listOfSpecies[p2].status = TSpeciesStatus.ssFloat then
                               dydt[individual.listOfSpecies[p2].mappingToYVariable] := dydt[individual.listOfSpecies[p2].mappingToYVariable] + reactionRate;
                            end;
      end;
end;

function simulate (individual : TIndividual; timeStart, stepSize : double; numberOfPoints : Integer) : T2DArray;
var i, j : integer; time : double;
    cvode : TCvode;
    y : TDArray;
begin
  setLength (result, numberOfPoints, individual.numberOfFloatingSpecies+1);  // +1 for time
  time:= timeStart;

  cvode := TCvode.Create (individual.numberOfFloatingSpecies);
  cvode.setModel (TModelClass.myfunc);
  TModelClass.individual := individual;

  setLength (y, individual.numberOfFloatingSpecies);
  for i := 0 to individual.numberOfFloatingSpecies - 1 do
      begin
      y[i] := initialConcentrations[i];
      result[0,i+1] := y[i];
      end;

  result[0,0] := time;

  for i := 1 to numberOfPoints - 1 do
      begin
      time := cvode.execute (y, time, stepSize);
      result[i,0] := time;
      for j := 0 to individual.numberOfFloatingSpecies - 1 do
          result[i,j+1] := y[j];
      end;
  cvode.Free;
end;


function computeFitnessOfIndividual (individual : TIndividual; objFunction : TObjectiveFunctionData) : double;
var cvode : TCvode;
    y : TDArray;
    i, j : integer;
    startTime, stepSize : double;
    numberOfPoints : integer;
    output : array of array of double;
    deviation : array of double;
    smallestDeviation : Double;
    estr : string;
    count : integer;
begin
  cvode := TCvode.Create (individual.numberOfFloatingSpecies);
  try
    cvode.setModel (TModelClass.myfunc);
    TModelClass.individual := individual;
    for i := 0 to individual.numberOfFloatingSpecies - 1 do
       cvode.rtol[i] := 1E-6;
    // Set up the y variable array
    setLength (y, individual.numberOfFloatingSpecies);
    for i := 0 to individual.numberOfFloatingSpecies - 1 do
        y[i] := initialConcentrations[i];

    // Run eval to get next time step and next values for y
    numberOfPoints := objFunction.numberOfPoints;
    stepSize := objFunction.stepSize;
    setLength (output, numberOfPoints, individual.numberOfFloatingSpecies);
    try
       startTime := 0.0;
       // Collect the first time point
       for j := 0 to individual.numberOfFloatingSpecies - 1 do
           output[0,j] := y[j];

       // Collect the rest of the points
       for i := 1 to numberOfPoints - 1 do
           begin
           startTime := cvode.execute (y, startTime, stepSize);
           for j := 0 to individual.numberOfFloatingSpecies - 1 do
               output[i,j] := y[j];
           end;

       // compute fitness with respect to each node
       setLength (deviation, individual.numberOfFloatingSpecies);
       smallestDeviation := 1E18;
       for j := 0 to individual.numberOfFloatingSpecies - 1 do
           begin
           deviation[j] := 0;
           for i := 0 to numberOfPoints - 1 do
               deviation[j] := deviation[j] + sqr (output[i, j] - objFunction.outputData[i]);
           count := 0;
           for i := 0 to individual.listOfSpecies.Count - 1 do
               begin
               if individual.listOfSpecies[i].status = TSpeciesStatus.ssFloat then
                  begin
                  individual.listOfSpecies[i].nodefitscore := deviation[count];
                  inc (count);
                  end;
               end;

           if smallestDeviation > deviation[j] then
              smallestDeviation := deviation[j];
           end;
    except
      on E: Exception do
         begin
         estr := e.Message;
         individual.fitness := 1E17;
         individual.badModel := True;
         Result := 1E12;
         exit;
         end;
    end;
    individual.fitness := smallestDeviation;
    result := smallestDeviation;

finally
   begin
   setlength (y, 0);
   cvode.Free;
   end;
end;
end;


// Also get rid of any bad models
procedure computeFitnessOfPopulation (population : TPopulation; objFunction : TObjectiveFunctionData);
var i : Integer;
    estr : string;
    d : double;
begin
  try
   for i := 0 to population.Count - 1 do
       begin
       if (population[i].numberOfFloatingSpecies <= 1) or (population[i].listOfReactions.Count = 0) then
           population[i].badModel := True
       else
           d := computeFitnessOfIndividual (population[i], objFunction);
        if d > 1E8 then
           population[i].badModel := True;
        if not population[i].badModel  then
           begin
           population[i].fitness := d;
           end;
        end;
      except
        on E: exception do
           writeln ('Fitness loop: ' + e.message);
      end;

      try
      for i := population.Count - 1 downto 0 do
          if population[i].badModel then
             population.Delete(i);
      except
        on E: exception do
          writeln ('Error in bad model delete, i = ', i);
      end;
end;


procedure saveSpeciesData (filename : string; generation : integer; population : TPopulation);
var i : integer;
    f : TextFile;
begin
  AssignFile (f, filename);
  if FileExists(fileName) then
     Append (f)
  else
     rewrite (f);

  write (f, generation, ',');
  for i := 0 to population.Count - 1 do
      write (f, population[i].speciesId, ',');

  writeln (f);
  CloseFile(f);

end;


function startMainLoop (var population : TPopulation; objFunction : TObjectiveFunctionData; var experimentInfo : TExperimentInfo) : Integer;
var i, n : integer; d : double;
    newPopulation : TPopulation;
    numberOfRouletteSelections : Integer;
    p1, p2, gen : Integer;
    individual : TIndividual;
    estr : string;
    f : TextFile;
    numMutations : TNumMutations;
    x : array of integer;
    y : integer;
begin
  numMutations.numRateConstantMutations := 0;
  numMutations.numAdditions := 0;
  numMutations.numDeletions := 0;

  if saveFitnessValuesToFile then
     begin
     AssignFile (f, experimentInfo.fileNamePrefix + 'fitness' + inttostr (experimentInfo.experimentId) + '.txt');
     rewrite (f);
     end;

  for gen := 0 to experimentInfo.maxGenerations - 1 do
      begin
      computeFitnessOfPopulation (population, objFunction);
      sortPopulation (population);

      if experimentInfo.savePopulationData then
         saveSpeciesData (experimentInfo.savePopulationFilename, gen, population);

      // Avoid running for too long if the fitness doesn't
      // decrease much as its not likely to decrease anymore
      if (population[0].fitness > 800) and (gen > 600) then
           begin
           result := gen;
           Break;
           end;

      // If the fitness is low and we are past 800 generations, then I recommend breaking out
      if (population[0].fitness < 30) and (gen > 800) then
           begin
           result := gen;
           Break;
           end;

      if gen mod intervalForOutput = 0 then
         begin
         if experimentInfo.showOutput then
            begin
            writeln ('Generation:', gen:3, ', fitness=', population[0].fitness:8:3,
            '  #P=', population.Count:3,
            ' #Fs = ', population[0].numberOfFloatingSpecies:2, ' #Rs=', population[0].listOfReactions.Count:2,
            ' #Ms=', numMutations.numRateConstantMutations:5, ' #As=', numMutations.numAdditions:5, ' #Ds=', numMutations.numDeletions:5);
            end
         else
            begin
            writeln ('Generation:', gen:4, ', fitness:', population[0].fitness:9:3);
            end;
         if saveFitnessValuesToFile then
            writeln (f, gen, population[0].fitness:9:3);
         end;

        // Default threshold is 10
        if population[0].fitness < fitnessThresholdToStop then
           begin
           result := gen;
           Break;
           end;

      // Create the next population

      newPopulation := TPopulation.Create;

      // Copy over the top n
      n := trunc (experimentInfo.percentageElite*population.Count);
      for i := 0 to n - 1 do
          newPopulation.Add (TIndividual.clone (population[i]));

     for i := experimentInfo.maxPopulationSize downto newPopulation.Count - 1 do
         begin
         p1 := math.RandomRange(0, population.Count);
         p2 := math.RandomRange(0, population.Count);
         if population[p1].fitness < population[p2].fitness then
            newPopulation.Add (TIndividual.clone (population[p1]))
         else
            newPopulation.Add (TIndividual.clone (population[p2]));
         end;

     if addRandomNetworks then
        begin
       // Create a few random ones
       for i := 0 to numRandomNetworks - 1 do
           newPopulation.Add(TIndividual.Create (experimentInfo));
       end;

     try
       for i := 1 to newPopulation.Count - 1 do
           begin
           // Don't mutate the best individual, hence we start at 1
           p1 := Math.RandomRange(1, newPopulation.Count);
           individual := mutateIndividual (newPopulation[p1], experimentInfo, numMutations);
           if not individual.badModel then
              begin
              newPopulation.Delete (p1);
              newPopulation.Add (individual);
              end;
           end;

       population.Free;
       population := newPopulation;
     except
       on E: Exception do
          writeln ('Mutation loop: ' + e.message)
     end;
  end;
  result := gen;
  if saveFitnessValuesToFile then
     closefile (f);
 end;

end.


      //numberOfRouletteSelections := population.Count div 2;
      //for i := 0 to numberOfRouletteSelections - 1 do
            //begin
            // Pick two random individuals (not including the first one
            // p1 := math.RandomRange(1, population.Count);
           // p2 := math.RandomRange(1, population.Count);
           // Copy over the best of the two
           // if population[p1].fitness < population[p2].fitness then
           //    newPopulation.Add (TIndividual.clone (population[p1]))
           // else
           //    newPopulation.Add (TIndividual.clone (population[p2]));
           // end;

