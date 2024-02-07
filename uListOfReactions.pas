unit uListOfReactions;

interface

Uses Classes, SysUtils, Generics.Collections;

type
   TReactionType = (rtUniUni, rtBiUni, rtUniBi, rtBiBi);

   TReaction = class (TObject)
      rateConstant : double;
      reactionType : TReactionType;
      connections : array[0..3] of integer;
      constructor Create;
      destructor Destroy; override;
   end;

   TListOfReactions = TObjectList<TReaction>;

   function copyListOfReactions (listOfReactions : TListOfReactions) : TListOfReactions;

implementation

constructor TReaction.Create;
begin
  inherited Create;
  rateConstant := 0.1;
end;


destructor TReaction.Destroy;
begin
  //writeln ('Free');
  inherited;
end;


function copyListOfReactions (listOfReactions : TListOfReactions) : TListOfReactions;
var i : integer;
    reaction : TReaction;
begin
  result := TListOfReactions.Create;
  for i := 0 to listOfReactions.Count - 1 do
      begin
      reaction := TReaction.Create;
      reaction.rateConstant := listOfReactions[i].rateConstant;
      reaction.reactionType := listOfReactions[i].reactionType;
      reaction.connections[0] := listOfReactions[i].connections[0];
      reaction.connections[1] := listOfReactions[i].connections[1];
      reaction.connections[2] := listOfReactions[i].connections[2];
      reaction.connections[3] := listOfReactions[i].connections[3];
      result.add (reaction);
      end;
end;

end.
