object Form1: TForm1
  Left = 291
  Top = 67
  Width = 1018
  Height = 940
  Caption = 'Form1'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Panel1: TPanel
    Left = 0
    Top = 0
    Width = 100
    Height = 901
    Align = alLeft
    TabOrder = 0
    object Button1: TButton
      Left = 8
      Top = 8
      Width = 75
      Height = 25
      Caption = 'Button1'
      TabOrder = 0
      OnClick = Button1Click
    end
    object Button2: TButton
      Left = 16
      Top = 72
      Width = 75
      Height = 25
      Caption = 'Button2'
      TabOrder = 1
    end
    object TrackBar1: TTrackBar
      Left = 24
      Top = 104
      Width = 45
      Height = 785
      Max = 255
      Orientation = trVertical
      Frequency = 1
      Position = 0
      SelEnd = 0
      SelStart = 0
      TabOrder = 2
      TickMarks = tmBottomRight
      TickStyle = tsAuto
      OnChange = TrackBar1Change
    end
  end
  object Panel2: TPanel
    Left = 100
    Top = 0
    Width = 902
    Height = 901
    Align = alClient
    TabOrder = 1
    object Image1: TImage
      Left = 1
      Top = 1
      Width = 900
      Height = 899
      Align = alClient
    end
  end
end
