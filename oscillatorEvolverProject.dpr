program oscillatorEvolverProject;

{$APPTYPE CONSOLE}

{$R *.res}

uses
  Winapi.Windows,
  Classes,
  Vcl.clipbrd,
  System.SysUtils,
  System.Types,
  StrUtils,
  IOUtils,
  System.Diagnostics,
  System.Threading,
  nVector_h in 'nVector_h.pas',
  uGlobal in 'uGlobal.pas',
  uGlobal_Types in 'uGlobal_Types.pas',
  uIndividual in 'uIndividual.pas',
  uListOfReactions in 'uListOfReactions.pas',
  uListOfSpecies in 'uListOfSpecies.pas',
  uMainLoop in 'uMainLoop.pas',
  uModel in 'uModel.pas',
  uMutate in 'uMutate.pas',
  uObjectiveFunction in 'uObjectiveFunction.pas',
  uPopulation in 'uPopulation.pas',
  uArgumentParser in 'uArgumentParser.pas',
  uCvode in 'uCvode.pas';

var
  objFunction: TObjectiveFunctionData;
  i : integer;
  networkStr : string;
  rect : TSmallRect;
  Coord : TCoord;
  experimentId : integer;
  //cvode : TCvode;
  Parser: TArgumentParser;
  ParseResult: TParseResult;
  seed : LongInt;
  iterations : integer;
  reducedModel : TIndividual;
  I1 : integer;
  successFileHandle : TextFile;
  sucessfullCount : integer;
  sw : TStopwatch;
  experimentInfo : TExperimentInfo;
  tryLargePopulationFirst : boolean;

procedure WaitForInput;
var
  InputRec: TInputRecord;
  NumRead: Cardinal;
  OldKeyMode: DWORD;
  StdIn: THandle;
begin
  StdIn := GetStdHandle(STD_INPUT_HANDLE);
  GetConsoleMode(StdIn, OldKeyMode);
  SetConsoleMode(StdIn, 0);
  repeat
    ReadConsoleInput(StdIn, InputRec, 1, NumRead);
  until (InputRec.EventType and KEY_EVENT <> 0) and InputRec.Event.KeyEvent.bKeyDown;
  SetConsoleMode(StdIn, OldKeyMode);
end;

procedure loadSettings;
var f : TextFile;
    astr : string;
    alist : TStringDynArray;
begin
  if FileExists ('setup.txt') then
     begin
     AssignFile (f, 'setup.txt');
     reset (f);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_maxGenerations := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_maxPopulationSize := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_probabilityOfUniUni := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_probabilityOfBiUni := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_probabilityOfUniBi := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_probabilityOfBiBi := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_percentageElite := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     //initialConcentrations := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_maxNumberOfSpecies := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_minNumberOfSpecies := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_maxNumberOfReactions := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_minNumberOfReactions := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_maxInitialRateConstant := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_rateConstantMutation := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_addReactionMutation := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     default_deleteReactionMutation := strtofloat (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     intervalForOutput := strtoint (alist[1]);

     readln (f, astr);
     aList := splitString (astr, ' ');
     fitnessThresholdToStop := strtofloat (alist[1]);

     closeFile (f);
     end
  else
     writeln ('******* Settings file not found ********');
end;


procedure printHelp;
begin
  writeln ('Evolver ' + version);
  writeln (' -h This help');
  writeln (' -stress Stress cvode');
  writeln (' -s value Set seed');
  writeln (' -sg Save fitness values to a file');
  writeln (' -mp value set maximum population size');
  writeln (' -pre value String prefix to use for output files');
  writeln (' -e value Number of evolution experiments to run');
  writeln (' -res value Specify the name of the results directory');
  writeln (' -v Make output more verbose');
  writeln (' -pop filename Output data to the file filename that describes the species distriution over generations');
  writeln (' -w Wait for the user to hit a key to exit program');
end;


// Find reactions that don't contribute to the fitness and remove them
function reduceModel (individual: TIndividual) : TIndividual;
var fitness, newFitness : double; copyOfIndividual : TIndividual;
    i, j : integer; attempt : Boolean;
    oldRateConstant : double;
begin
  fitness := computeFitnessOfIndividual (individual, objFunction);
  copyOfIndividual := TIndividual.clone (individual);
  attempt := True;
  //while attempt do
        begin
        attempt := False;
        for j := 0 to copyOfIndividual.listOfReactions.Count - 1 do
            begin
            oldRateConstant := copyOfIndividual.listOfReactions[j].rateConstant;
            copyOfIndividual.listOfReactions[j].rateConstant := 0;
            newFitness := computeFitnessOfIndividual (copyOfIndividual, objFunction);
            if newFitness > 10*fitness then
               begin
               // keep reaction
               copyOfIndividual.listOfReactions[j].rateConstant := oldRateConstant;
               end
            else
               attempt := True;
            end;
        // Go through reactions deleting zero rate cosntant reactions
        for j := copyOfIndividual.listOfReactions.Count - 1 downto 0 do
            if copyOfIndividual.listOfReactions[j].rateConstant = 0 then
               removeNthReaction (copyOfIndividual, j);
         end;
  result := copyOfIndividual;
end;


procedure outputResults (population : TPopulation; sessionId : string; experimentInfo : TExperimentInfo);
begin
   try
     writeln;
     writeln ('Fitness of best individual = ', computeFitnessOfIndividual(population[0], objFunction):6:3);

     networkStr := 'import tellurium as te' + sLineBreak + sLineBreak;
     networkStr := networkStr + 'r = te.loada("""' + sLineBreak;
     networkStr := networkStr +  '// Version of software = ' + version + sLineBreak;
     networkStr := networkStr + '// Seed = ' + inttostr (experimentInfo.seedValue) + sLineBreak;
     networkStr := networkStr + '// maxPopulationSize = ' + inttostr (experimentInfo.maxPopulationSize) + sLineBreak;
     networkStr := networkStr + '// maxNumberOfSpecies = ' + inttostr (experimentInfo.maxNumberOfSpecies) + sLineBreak;
     networkStr := networkStr + '// minNumberOfSpecies = ' + inttostr (experimentInfo.minNumberOfSpecies) + sLineBreak;
     networkStr := networkStr + '// maxNumberOfReactions = ' + inttostr (experimentInfo.maxNumberOfReactions) + sLineBreak;
     networkStr := networkStr + '// minNumberOfReactions = ' + inttostr (experimentInfo.minNumberOfReactions) + sLineBreak;

     networkStr := networkStr + '// probabilityOfUniUni = ' + floattostr (experimentInfo.probabilityOfUniUni) + sLineBreak;
     networkStr := networkStr + '// probabilityOfBiUni = ' + floattostr (experimentInfo.probabilityOfBiUni) + sLineBreak;
     networkStr := networkStr + '// probabilityOfUniBi = ' + floattostr (experimentInfo.probabilityOfUniBi) + sLineBreak;
     networkStr := networkStr + '// probabilityOfBiBi = ' + floattostr (experimentInfo.probabilityOfBiBi) + sLineBreak;
     networkStr := networkStr + '// percentageClone = ' + floattostr (experimentInfo.percentageElite) + sLineBreak;
     networkStr := networkStr + '// initialConcentrations = ' + sLineBreak;
     networkStr := networkStr + '// ';
     for i := 0 to 19 do
         networkstr := networkStr + floattostr (initialConcentrations[i]) + ', ';
     networkStr := networkStr + sLineBreak;

     networkStr := networkStr + '// maxInitialRateConstant = ' + floattostr (experimentInfo.maxInitialRateConstant) + sLineBreak;
     networkStr := networkStr + '// rateConstantMutation = ' + floattostr (experimentInfo.rateConstantMutation) + sLineBreak;
     networkStr := networkStr + '// addReactionMutation = ' + floattostr (experimentInfo.addReactionMutation) + sLineBreak;
     networkStr := networkStr + '// deleteReactionMutation = ' + floattostr (experimentInfo.deleteReactionMutation) + sLineBreak;

     writeln ('Reducing size of network:');
     write ('Size before and after: ', IntToStr(population[0].listOfReactions.Count));
     reducedModel := reduceModel (population[0]);
     writeln (', ' + IntToStr(reducedModel.listOfReactions.Count) + sLineBreak);
     networkStr := networkStr + reducedModel.convertToString + sLineBreak;
     writeln (reducedModel.convertToString + sLineBreak);

     networkStr := networkStr + sLineBreak;
     networkStr := networkStr + '//------------------------------------------' + sLineBreak;
     networkStr := networkStr + '// Best fitness = ' + format ('%8.2f', [population[0].fitness]) + sLineBreak;
     networkStr := networkStr + '//------------------------------------------' + sLineBreak;

     networkStr := networkStr + sLineBreak + '// Number of reactions = ' + inttostr (reducedModel.listOfReactions.Count) + sLineBreak;

     networkStr := networkStr + '// Number of floating species = ' + inttostr (reducedModel.numberOfFloatingSpecies) + sLineBreak;
     networkStr := networkStr + '// Iterations carried out = ' + inttostr (iterations);
     networkStr := networkStr + sLinebreak + '""")' + sLineBreak;

     networkStr := networkStr + sLineBreak + 'm = r.simulate (0, 15, 100)' + sLineBreak;
     networkStr := networkStr +  'r.plot()' + sLineBreak;

     Clipboard.AsText := networkStr;

     // Don't bother outputing anything with a fitness greater than 100
     if population[0].fitness < 100 then
        TFile.WriteAllText ('./results/' + experimentInfo.fileNamePrefix + 'M' + sessionId + '_' + inttostr (experimentId) + '_' + inttostr (trunc(population[0].fitness)) + '.txt', networkStr);
   except
     writeln ('Error in writing out best individual');
   end;

   if experimentInfo.showOutput then
      begin
      writeln ('Best Fitness = ', format ('%8.2f', [population[0].fitness]));
      writeln ('Worst Fitness = ', format ('%8.2f', [population[population.Count-1].fitness]));
      writeln ('Number of iterations = ', iterations);
      end;
end;


// Run a single evolution experimeent
procedure runExperiment (experimentNumber : integer; seed : LongInt; experimentInfo : TExperimentInfo);
var
   population, tmpPop : TPopulation;
   tmp : integer;
begin
  if experimentInfo.showOutput then writeln ('Setup the random number seed');
  if seed = 0 then
     begin
     if experimentInfo.showOutput then writeln ('Create seed based on tickcount');
     I1 := (RandSeed shl 8);
     experimentInfo.seedValue := (Cardinal (I1) or GetCurrentProcessID) xor GetTickCount;
     RandSeed := integer (experimentInfo.seedValue);
     end
  else
     begin
     if experimentInfo.showOutput then writeln ('Use the seed from command line argument');
     experimentInfo.seedValue := seed;
     RandSeed := experimentInfo.seedValue;
     end;

  writeln;
  writeln ('Start of run: ' + inttostr (experimentId):3);
  writeln ('-----------------');
  if experimentInfo.showOutput then
     writeln ('Create population');

  tryLargePopulationFirst := False;
  if tryLargePopulationFirst then
     begin
     tmp := experimentInfo.maxPopulationSize;
     experimentInfo.maxPopulationSize := 10*tmp;
     tmpPop := createPopulation (experimentInfo);
     sortPopulation(tmpPop);
     population := TPopulation.Create;
     for var k := 0 to tmp - 1 do
         population.Add (TIndividual.clone(tmpPop[k]));
     experimentInfo.maxPopulationSize := tmp;
     end
  else
     population := createPopulation (experimentInfo);

  if injectTestingNetwork then
     begin
     population.Insert(0, TIndividual.CreateTestNetwork);
     population.Insert(0, TIndividual.CreateTestNetwork);
     population.Insert(0, TIndividual.CreateTestNetwork);
     population.Insert(0, TIndividual.CreateTestNetwork);
     population.Insert(0, TIndividual.CreateTestNetwork);
     end;

  if experimentInfo.showOutput then writeln ('Best initial network in the population is:');
  computeFitnessOfPopulation (population, objFunction);
  sortPopulation (population);
  if experimentInfo.showOutput then
     begin
     writeln (population[0].convertToString);
     writeln;
     writeln ('Fitness of best individual at the start = ', population[0].fitness:8:4);
     writeln;
     end;

  try
    if experimentInfo.showOutput then
       begin
       writeln ('#P(Population), #Fs(Species), #Rs(Reactions), #Ms(Mutationa), #As(Additions), #Ds(Deletions)');
       writeln;
       end;
    iterations := startMainLoop (population, objFunction, experimentInfo);
  except
    on E: exception do
       writeln ('Exception in startMainLoop: ' + e.Message);
  end;

  computeFitnessOfPopulation (population, objFunction);
  sortPopulation (population);

  OutputResults (population, generateSessionId, experimentInfo);

  // Record we found an oscillator, just a convenince, not essential
  if population[0].fitness < 12 then
     begin
     AssignFile (successFileHandle, 'successful.txt');
     rewrite(successFileHandle);

     inc (sucessfullCount);
     writeln (successFileHandle, sucessfullCount);
     Flush (successFileHandle);
     CloseFile (successFileHandle);
     end;

  population.Free;
end;



begin
  try
    Rect.Left := 2;
    Rect.Top := 1;
    Rect.Right := 130;
    Rect.Bottom := 60;
    Coord.X := Rect.Right + 1 - Rect.Left;
    Coord.y := Rect.Bottom + 1 - Rect.Top;
    seed := 0;
    sucessfullCount := 0;
    //SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), Coord);
    //SetConsoleWindowInfo(GetStdHandle(STD_OUTPUT_HANDLE), True, Rect);

    Parser := TArgumentParser.Create;
    Parser.AddArgument ('-h', saBool);
    Parser.AddArgument ('-stress', saBool);
    Parser.AddArgument ('-sg', saBool);

    Parser.AddArgument ('-s', saStore);
    Parser.AddArgument ('-pre', saStore);
    Parser.AddArgument ('-e', saStore);

    Parser.AddArgument ('-mp', saStore);
    Parser.AddArgument ('-res', saStore);

    Parser.AddArgument ('-v', saBool);
    Parser.AddArgument ('-pop', saStore);

    Parser.AddArgument ('-w', saBool);

    ParseResult := Parser.ParseArgs;
    try
      if ParseResult.HasArgument('h') then
         begin
         printHelp;
         exit;
         end;
      if ParseResult.HasArgument('stress') then
         begin
         //Stress Test
         ReportMemoryLeaksOnShutdown := True;
         writeln ('Running 1000,000 cvode instantiations with 20 state variables.....');
         for i := 0 to 1000000 do
             begin
             //cvode := TCvode.Create (20);
             //cvode.Free;
             end;
         writeln ('Done test, press any key to continue');
         WaitForInput;
         exit;
         end;
      if ParseResult.HasArgument('s') then
         begin
         seed := StrToInt64 (ParseResult.GetValue ('s'));
         Writeln('Setting seed to: ', seed);
         end;
      if ParseResult.HasArgument('mp') then
         begin
         default_maxPopulationSize := StrToInt (ParseResult.GetValue ('mp'));
         end;
      if ParseResult.HasArgument('pre') then
         begin
         default_fileNamePrefix := ParseResult.GetValue ('pre') + '_';
         end;
      if ParseResult.HasArgument('e') then
         begin
         numberOfExperiments := strtoint (ParseResult.GetValue ('e'));
         end;
      if ParseResult.HasArgument('v') then
         begin
         default_showOutput := True;
         end;
      if ParseResult.HasArgument('w') then
         begin
         default_waitForKeyPress := True;
         end;
      if ParseResult.HasArgument('pop') then
         begin
         default_savePopulationData := True;
         default_savePopulationFilename := ParseResult.GetValue ('pop');
         end;
      if ParseResult.HasArgument('sg') then
         saveFitnessValuesToFile := True;

      if ParseResult.HasArgument('res') then
         default_resultsFolderName := ParseResult.GetValue ('res');
    finally
      Parser.Free;
      ParseResult.Free;
    end;

    // Make sure there is a results folder, if not create one
    if not DirectoryExists('./' + default_resultsFolderName) then
       begin
       writeln ('Create results directory');
       TDirectory.CreateDirectory('./' + default_resultsFolderName);
       end;

    writeln;
    writeln ('Software version: ' + version);
    writeln ('---------------------');
    if not loadDLL() then
       begin
       writeln ('Failed to load cvode, type any key to exit');
       readln;
       exit;
       end
    else
       writeln ('Cvode loaded');

    loadSettings;
    writeln ('Settings loaded');

    objFunction := setObjectiveFunction;
    sw := TStopwatch.StartNew; // To time how long it takes

    if FileExists(default_savePopulationFilename) then
       DeleteFile(default_savePopulationFilename);

    //TParallel.For(0, numberOfExperiments - 1, procedure(experimentId : integer)
    for experimentId := 0 to numberOfExperiments - 1 do
        begin
        //var experimentInfo : TExperimentInfo;
        experimentInfo.fileNamePrefix := '';
        experimentInfo.showOutput := default_showOutput;
        experimentInfo.experimentId := experimentId;
        experimentInfo.maxGenerations := default_maxGenerations;
        experimentInfo.maxPopulationSize := default_maxPopulationSize;
        experimentInfo.maxNumberOfSpecies := default_minNumberOfSpecies;
        experimentInfo.maxNumberOfReactions := default_maxNumberOfReactions;
        experimentInfo.maxInitialRateConstant := default_maxInitialRateConstant;
        experimentInfo.minNumberOfSpecies := default_minNumberOfSpecies;
        experimentInfo.minNumberOfReactions := default_maxNumberOfReactions;
        experimentInfo.rateConstantMutation := default_rateConstantMutation;
        experimentInfo.addReactionMutation := default_addReactionMutation;
        experimentInfo.deleteReactionMutation := default_deleteReactionMutation;
        experimentInfo.percentageElite := default_percentageElite;
        experimentInfo.probabilityOfUniUni := default_probabilityOfUniUni;
        experimentInfo.probabilityOfUniBi := default_probabilityOfUniBi;
        experimentInfo.probabilityOfBiUni := default_probabilityOfBiUni;
        experimentInfo.probabilityOfBiBi := default_probabilityOfBiBi;
        experimentInfo.resultsFolderName := default_resultsFolderName;
        experimentInfo.savePopulationData := default_savePopulationData;
        experimentInfo.savePopulationFilename := default_savePopulationFilename;

        runExperiment (experimentId, seed, experimentInfo);
        //end);     // End of experiment parallel loop
        end;

    objFunction.Free;

    writeln ('Time taken in hours, minutes and seconds: ', FormatDateTime('hh:nn:ss',sw.ElapsedMilliseconds/MSecsPerDay));

    if default_waitForKeyPress then
       begin
       writeln ('Done. Hit any key to continue');
       readln;
       end;
  except
    on E: Exception do
      begin
      Writeln(E.ClassName, ': ', E.Message);
      //readln;
      end;
  end;
end.


//    output := simulate (population[0], 0, 0.05, 20);
//    for i := 0 to length (output) - 1 do
//        begin
//        for j := 0 to length (output[0]) - 1 do
//            write (output[i,j]:6:5, ', ');
//        writeln;
//        end;
//
//    networkStr := population[0].convertToString + sLineBreak;
//    writeln (population[0].convertToString);
//
//    for i := 0 to population[0].listOfSpecies.Count - 1 do
//        begin
//        networkStr := networkStr + 'S' + inttostr (i) + ' = 1;';
//        Writeln('S' + inttostr (i) + ' = 1;');
//        end;
//    Clipboard.AsText := networkStr;
//
//    Readln;
//    exit;
