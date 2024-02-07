program cvode262_Project;

uses
  Vcl.Forms,
  nVector_h in 'src\nVector_h.pas',
  uCvode in 'src\uCvode.pas',
  ufMain in 'src\ufMain.pas' {frmMain};

{$R *.res}

begin
  Application.Initialize;
  Application.MainFormOnTaskbar := True;
  Application.CreateForm(TfrmMain, frmMain);
  Application.Run;
end.
