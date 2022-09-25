//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
 
#include "Unit1.h"
#include <math>
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 *Form1;

long double p=0.001; // коэффициент задающий начальную интенсивность концентраций дефектов
int h = 1;
long double dt = 1 * pow(10, -3); // дельта времени
double a = 1*pow(10,-2)/dt, b = 1.5*pow(10,-2)/dt; // коэффициенты в расчете химического потенциала
double P=1.5*pow(10,-4),K=1.5*pow(10,-3); // коэффициенты скорости генерации дефектов и рекомбинации
int s_b = 3*b; 
//double a = 0.5, b = 1;

double tm = 1000; // общее время процесса
long double t = tm / dt;
//---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner)
{

   // Инициализация массивов


   (long double *)m=(long double *)new long double[sizeof(long double)*(long(t)+1)];
   for(int i=0;i<long(t);i++) m[i]=0;

   for(int i=1;i<long(t);i++) m[i]=m[i-1]+dt;

   n=50; // размер матрицы концентраций

   (unsigned int *)vac=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)vac[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) vac[i][j]=0;


   (unsigned int *)d2vac=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)d2vac[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) d2vac[i][j]=0;


   (unsigned int *)sia=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)sia[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) sia[i][j]=0;

   (unsigned int *)d2sia=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)d2sia[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) d2sia[i][j]=0;


   (unsigned int *)fpso=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)fpso[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) fpso[i][j]=0;

   (unsigned int *)mu=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)mu[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) mu[i][j]=0;


   (unsigned int *)d2mu=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)d2mu[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) d2mu[i][j]=0;

   (unsigned int *)vac_new=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)vac_new[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) vac_new[i][j]=0;

   (unsigned int *)vac_grf=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) (unsigned int *)vac_grf[i]=(unsigned int *)new int[sizeof(int)*(n+1)];
   for(int i=0;i<n;i++) for(int j=0;j<n;j++) vac_grf[i][j]=0;

}
//---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
        randomize();
        
        //формирование матриц концентраций дефектов, интенсивностью р , случайным образом  

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
                        long double ee=random(1000000)/1000000.0;
                        vac[i][j] = p*ee;
                        ee=random(1000000)/1000000.0;
                        sia[i][j] = vac[i][j];
		}
	}

        DrawAllConc(); // отображение на экране

	for(long k=0;k<long(t);k++)
	{
                // Расчет втророй производной концентраций вакансий методом центральной разности
		for(int j=0;j<n;j++)
		{
			for(int i=0;i<n;i++)
			{
                                int dip=i+1,dim=i-1,djp=j+1,djm=j-1;
                                if(dip>n-1) dip=0;
                                if(dim<0) dim=n-1;
                                if(djp>n-1) djp=0;
                                if(djm<0) djm=n-1;
				d2vac[i][j]=(vac[dip][j]-4*vac[i][j]+vac[dim][j]+vac[i][djp]+vac[i][djm])/(h*h);
			}
		}
                // Расчет второй производной концентраций SIA методом центральной разности

        	for(int j=0;j<n;j++)
		{
			for(int i=0;i<n;i++)
			{
                                int dip=i+1,dim=i-1,djp=j+1,djm=j-1;
                                if(dip>n-1) dip=0;
                                if(dim<0) dim=n-1;
                                if(djp>n-1) djp=0;
                                if(djm<0) djm=n-1;
				d2sia[i][j]=(vac[dip][j]-4*vac[i][j]+vac[dim][j]+vac[i][djp]+vac[i][djm])/(h*h);
			}
		}

                // Расчет функции потенциала свободной энергии

		for(int i=0;i<n;i++)
		{
			for(int j=0;j<n;j++)
			{
				fpso[i][j]=2*b*(2*pow(vac[i][j],3)-3*pow(vac[i][j],2)+vac[i][j]);
                                //fpso[i][j] = 2 * b*(2 * pow(vac[i][j], 3) - 3 * pow(vac[i][j], 2) + vac[i][j]);
                                //fpso[i][j] = b*(pow(vac[i][j], 3) - vac[i][j]);
			}
		}

                // Расчет химического потенциала
 
		for(int i=0;i<n;i++) for(int j=0;j<n;j++) mu[i][j]=fpso[i][j]-a*(d2vac[i][j]);

                // Расчет второй производной химического потенциала методом центральной разности

		for(int j=0;j<n;j++)
		{
			for(int i=0;i<n;i++)
			{
                                int dip=i+1,dim=i-1,djp=j+1,djm=j-1;
                                if(dip>n-1) dip=0;
                                if(dim<0) dim=n-1;
                                if(djp>n-1) djp=0;
                                if(djm<0) djm=n-1;
				d2mu[i][j]=(mu[dip][j]-4*mu[i][j]+mu[dim][j]+mu[i][djp]+mu[i][djm])/(h*h);
			}
		}


                long double sum_vac=0,sum_sia=0;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
                                sum_vac+=vac[i][j];
                                sum_sia+=sia[i][j];
			}
		}
                sum_vac=sum_vac/(n*n);
                sum_sia=sum_sia/(n*n);


                // Расчет нового значения концентрации с учетом скорости генерации дефектов и рекомбинации


		for (int j = 0; j < n; j++)
		{
			for (int i = 0; i < n; i++)
			{
                                vac_new[i][j] = P*(1-vac[i][j])+vac[i][j]-K*vac[i][j]*sia[i][j] + dt*d2mu[i][j];
                                sia[i][j] = P*(1-sia[i][j])+sia[i][j]-K*vac[i][j]*sia[i][j] + dt*d2sia[i][j];
                                vac[i][j] = vac_new[i][j];
			}
		}

                if(!(k%1000)) DrawAllConc();

        }

}

//---------------------------------------------------------------------------
void __fastcall TForm1::DrawConc(int x, int y, int color)
{
    int dx,dy;

    if(x==0) dx=4;
    else dx=x*9;

    if(y==0) dy=Image1->Height-4;
    else dy=Image1->Height-y*9;

    Image1->Canvas->Pen->Color=RGB(color, color, 0);
    Image1->Canvas->Brush->Color = RGB(color, color, 0);
    Image1->Canvas->Rectangle(dx-4,dy-4,dx+4,dy+4);
}

//---------------------------------------------------------------------------
void __fastcall TForm1::DrawAllConc(void)
{
        long double min=99999999;
        long double max=0;
        long double minus=0;


	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if(vac[i][j]<minus) minus=vac[i][j];
		}
	}

        if(fabs(minus)>0)
        {
  	   for (int i = 0; i < n; i++)
  	   {
		for (int j = 0; j < n; j++)
		{
			vac_grf[i][j]=vac[i][j]+fabs(minus);
		}
	   }
        } else
        {
  	   for (int i = 0; i < n; i++)
  	   {
		for (int j = 0; j < n; j++)
		{
			vac_grf[i][j]=vac[i][j];
		}
	   }
        }

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if(vac_grf[i][j]<min) min=vac_grf[i][j];
                        if(vac_grf[i][j]>max) max=vac_grf[i][j];
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
                        long double temp=(max-min);
			vac_grf[i][j]=255*(vac_grf[i][j]-min)/(max-min);
                        DrawConc(i,j,int(vac_grf[i][j]));
		}
	}
        Application->ProcessMessages();

}

//---------------------------------------------------------------------------
void __fastcall TForm1::TrackBar1Change(TObject *Sender)
{
       int x=50;
       int y=50;

       DrawConc(x,y,TrackBar1->Position);
}
//---------------------------------------------------------------------------
