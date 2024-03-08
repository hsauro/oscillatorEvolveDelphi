unit uGlobal;

interface

uses Generics.Collections;

type
    T2DArray = array of array of double;

    TNumMutations = record
        numRateConstantMutations : integer;
        numAdditions : integer;
        numDeletions : integer;
    end;

    TExperimentInfo = record

       experimentId : integer;
       fileNamePrefix : string;
       seedValue : Cardinal;

       maxGenerations : integer;

       maxPopulationSize : integer;
       maxNumberOfSpecies : integer;
       minNumberOfSpecies : integer;
       maxNumberOfReactions : integer;
       minNumberOfReactions : integer;
       maxInitialRateConstant : double;
       rateConstantMutation : double;
       addReactionMutation : double;
       deleteReactionMutation : double;

       percentageElite : double;

       probabilityOfUniUni : double;
       probabilityOfBiUni : double;
       probabilityOfUniBi : double;
       probabilityOfBiBi : double;

       showOutput : boolean;
       resultsFolderName : string;
       savePopulationData : boolean;
       savePopulationFilename : string;
    end;

var
   version : string = '1.51';

   default_maxPopulationSize : integer = 100;
   default_maxGenerations : integer = 1600;

   default_probabilityOfUniUni : double = 0.4;
   default_probabilityOfBiUni : double = 0.15;
   default_probabilityOfUniBi : double = 0.15;
   default_probabilityOfBiBi : double = 0.3;

   default_percentageElite : double = 0.4;

   initialConcentrations : array of double;
   default_maxNumberOfSpecies : integer = 5;
   default_minNumberOfSpecies : integer = 5;
   default_maxNumberOfReactions : integer = 4;
   default_minNumberOfReactions : integer = 4;
   default_maxInitialRateConstant : double = 50;

   default_rateConstantMutation : double = 0.5;
   default_addReactionMutation : double = 0.225;
   default_deleteReactionMutation : double = 0.275;

   default_showOutput : boolean = False;   // Shows additional runtime information if true
   default_savePopulationData : boolean = False;
   default_savePopulationFilename : string = '';   // Filename where program can dump all information about the population over time
   default_waitForKeyPress : boolean = False;

   saveFitnessValuesToFile : boolean = False;
   intervalForOutput  : integer = 10;
   fitnessThresholdToStop  : double = 50.0;

   numRandomNetworks : integer = 10;
   addRandomNetworks : boolean = False;
   injectTestingNetwork : boolean = False;   // Only for testing
   numberOfExperiments : integer = 1;
   default_resultsFolderName : string = 'results';  // All generated models go here
   default_fileNamePrefix : string = '';

   function createNewSpeciesId : integer;
   function generateSessionId : string;

implementation

var speciesCounter : integer;


function createNewSpeciesId : integer;
begin
  inc (speciesCounter);
  result := speciesCounter - 1;
end;


procedure initializeConcentrations;
var i : integer;
begin
  setlength (initialConcentrations, 20);
  initialConcentrations[0] :=  2.0;
  initialConcentrations[1] :=  5.0;
  initialConcentrations[2] :=  7.0;
  initialConcentrations[3] :=  10.0;
  initialConcentrations[4] :=  1.0;
  initialConcentrations[5] :=  0.5;
  initialConcentrations[6] :=  3.0;
  initialConcentrations[7] :=  7.0;
  initialConcentrations[8] :=  12.0;
  initialConcentrations[9] :=  4.0;
  initialConcentrations[10] :=  1.0;
  initialConcentrations[11] :=  3.0;
  initialConcentrations[12] :=  7.0;
  initialConcentrations[13] :=  9.0;
  initialConcentrations[14] :=  1.0;
  initialConcentrations[15] :=  7.0;
  initialConcentrations[16] :=  2.0;
  initialConcentrations[17] :=  8.0;
  initialConcentrations[18] :=  5.0;
  initialConcentrations[19] :=  10.0;
end;


// Each run of the program willtag the resulting model files with a four character sessionId
// Number of cobminations = (10 + 52)^6 = 56,800,235,584
function generateSessionId : string;
const
  letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789';
begin
  result := '';
  result := result + letters[random(length (letters)) + 1];
  result := result + letters[random(length (letters)) + 1];
  result := result + letters[random(length (letters)) + 1];
  result := result + letters[random(length (letters)) + 1];
  result := result + letters[random(length (letters)) + 1];
  result := result + letters[random(length (letters)) + 1];
  end;


initialization
  speciesCounter := 0;
  initializeConcentrations;
end.
