unit uListOfSpecies;

interface

Uses Generics.Collections, System.TypInfo;

type
  TSpeciesStatus = (ssBoundaryInput, ssBoundaryOutput, ssFloat, ssBlank, ssInput, ssOrphan);


  TSpecies = class

    concentration : double;
    status : TSpeciesStatus;
  	nodefitscore : double;
    // Used to determine if species is boundary or not once network is built
    leavings, enterings : integer;
    mappingToYVariable : Integer; // filled in when network is created, used by cvode

    constructor Create;
    //structor  Destroy; override;
    function statusToString : string;
  end;

  TListOfSpecies = TObjectList<TSpecies>;

  function copyListOfSpecies (listOfSpecies : TListOfSpecies) : TListOfSpecies;

implementation

constructor TSpecies.Create;
begin
  inherited Create;
  leavings := 0;
  enterings := 0;
  nodefitscore := 0.0;
end;

//destructor TSpecies.Destroy;
//begin
//  inherited;
//end;

function copyListOfSpecies (listOfSpecies : TListOfSpecies) : TListOfSpecies;
var i : integer;
    species : TSpecies;
begin
  result := TListOfSpecies.Create;
  for i := 0 to listOfSpecies.Count - 1 do
      begin
      species := TSpecies.Create;
      species.leavings := listOfSpecies[i].leavings;
      species.enterings := listOfSpecies[i].enterings;
      species.status := listOfSpecies[i].status;
      species.nodefitscore := listOfSpecies[i].nodefitscore;
      species.concentration := listOfSpecies[i].concentration;
      species.mappingToYVariable := listOfSpecies[i].mappingToYVariable;
      result.Add (species);
      end;
end;


function TSpecies.statusToString : string;
begin
  result := GetEnumName(TypeInfo(TSpeciesStatus), Integer(status));
end;

end.
