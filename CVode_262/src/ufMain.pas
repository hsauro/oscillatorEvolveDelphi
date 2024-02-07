unit ufMain;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, uCvode, Vcl.Grids,
  AdvObj, BaseGrid, AdvGrid;

type
  TfrmMain = class(TForm)
    btnRun: TButton;
    ListBox: TListBox;
    grid: TAdvStringGrid;
    procedure btnRunClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
    cvode : TCvode;
    procedure model (time : double; y, ydotv : TDArray);
  end;

var
  frmMain: TfrmMain;

implementation

{$R *.dfm}

procedure TfrmMain.model (time : double; y, ydotv : TDArray);
begin
  ydotv[0] := 77.27*(y[1]+y[0]*(1-2E-2*y[0]-y[1]));
  ydotv[1] := (y[2]-(1+y[0])*y[1])/77.27;
  ydotv[2] := 0.161*(y[0]-y[2])
end;

procedure TfrmMain.btnRunClick(Sender: TObject);
var y : TDArray;
    t : double;
    i, j : integer;
    n : Integer;
begin
  if not loadDLL then
     begin
     showmessage ('Failed to load cvode');
     exit;
     end;

  n := 3;
  setLength (y, n);
  y[0] := 1;
  y[1] := 0.5;
  y[2] := 0.6;

  cvode := TCvode.Create (3);
  cvode.setModel (model);
  t := 0;
  for i := 0 to 10 do
      begin
      grid.Cells[1,i+1] := FloatToStr(t);
      for j := 0 to n - 1 do
          grid.Cells[j+2, i+1] := FloatToStr(y[j]);
      listbox.Items.Add (FloatToStr(t) + ', ' + FloatToStr(y[0]) + ', ' + FloatToStr(y[1]) + ', ' + FloatToStr(y[2]));
      t := cvode.execute (y, t, 1);
      end;

  cvode.free;
end;

end.
