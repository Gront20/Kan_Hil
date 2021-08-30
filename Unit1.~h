//---------------------------------------------------------------------------

#ifndef Unit1H
#define Unit1H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <ComCtrls.hpp>
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
        TPanel *Panel1;
        TButton *Button1;
        TPanel *Panel2;
        TImage *Image1;
        TButton *Button2;
        TTrackBar *TrackBar1;
        void __fastcall Button1Click(TObject *Sender);
        void __fastcall TrackBar1Change(TObject *Sender);
private:	// User declarations
public:		// User declarations
        int n;
        long double **vac;
        long double **d2vac;

        long double **sia;
        long double **d2sia;
        long double **s_phi_new;

        long double **fpso;
        long double **mu;
        long double **d2mu;
        long double **vac_new;
        long double **vac_grf;
        long double *m;

        void __fastcall DrawConc(int x, int y, int color);
        void __fastcall DrawAllConc(void); 

        __fastcall TForm1(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
