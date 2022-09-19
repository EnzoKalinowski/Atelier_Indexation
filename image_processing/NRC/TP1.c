#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include <math.h>


#define TRESHOLD 15

int main(void){

	long nrh,nrl,nch,ncl;
	long nrh2,nrl2,nch2,ncl2;
	byte **I;
	byte **R;
	byte **Ix;
	byte **Iy;
	byte **G;
	byte **C;

	rgb8 **Icolor;

	byte **Rr;
	byte **Rg;
	byte **Rb;

	byte **IxR;
	byte **IyR;
	byte **Gr;
	byte **Cr;


 

	int H[3][3]={{1,1,1}, {1,1,1}, {1,1,1}};
	int SobelH[3][3]={{-1,0,1},{-2,0,2},{-1,0,1} };
	int SobelV[3][3]={{-1,-2,-1},{0,0,0},{1,2,1} };

	int total,totalR,totalG,totalB;

	int i,j;
	int u,v;

	I=LoadPGM_bmatrix("cubes.pgm",&nrl,&nrh,&ncl,&nch);
	Icolor=LoadPPM_rgb8matrix("trees.ppm",&nrl2,&nrh2,&ncl2,&nch2);

	R=bmatrix(nrl,nrh,ncl,nch);
	Ix=bmatrix(nrl,nrh,ncl,nch);
	Iy=bmatrix(nrl,nrh,ncl,nch);
	G=bmatrix(nrl,nrh,ncl,nch);
	C=bmatrix(nrl,nrh,ncl,nch);

	Rr=bmatrix(nrl,nrh,ncl,nch);
	Rg=bmatrix(nrl,nrh,ncl,nch);
	Rb=bmatrix(nrl,nrh,ncl,nch);

	

	for(i=nrl+1;i<nrh;i++)
	{	
		for(j=ncl+1;j<nch;j++)
		{
			//convolution reponse impulsionnelle
			total=I[i-1][j-1]*H[0][0] + I[i][j-1]*H[1][0] + I[i+1][j-1]*H[2][0] + I[i-1][j]*H[0][1] + I[i][j]*H[1][1] + I[i+1][j]*H[2][1] + I[i-1][j+1]*H[0][2] + I[i+1][j+1]*H[1][2] + I[i+1][j+1]*H[2][2];
			R[i][j]=total/9;

			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelH[0][0] + I[i][j-1]*SobelH[1][0] + I[i+1][j-1]*SobelH[2][0] + I[i-1][j]*SobelH[0][1] + I[i][j]*SobelH[1][1] + I[i+1][j]*SobelH[2][1] + I[i-1][j+1]*SobelH[0][2] + I[i+1][j+1]*SobelH[1][2] + I[i+1][j+1]*SobelH[2][2];
			Ix[i][j]=abs(total)/4;

			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelV[0][0] + I[i][j-1]*SobelV[1][0] + I[i+1][j-1]*SobelV[2][0] + I[i-1][j]*SobelV[0][1] + I[i][j]*SobelV[1][1] + I[i+1][j]*SobelV[2][1] + I[i-1][j+1]*SobelV[0][2] + I[i+1][j+1]*SobelV[1][2] + I[i+1][j+1]*SobelV[2][2];
			Iy[i][j]=abs(total)/4;

			//-------------------PPM(couleur)------------------------
			//Red
			//convolution Sobel horizontal
			totalR=Icolor[i-1][j-1].r*SobelH[0][0] + Icolor[i][j-1].r*SobelH[1][0] + Icolor[i+1][j-1].r*SobelH[2][0] + Icolor[i-1][j].r*SobelH[0][1] + Icolor[i][j].r*SobelH[1][1] + Icolor[i+1][j].r*SobelH[2][1] + Icolor[i-1][j+1].r*SobelH[0][2] + Icolor[i+1][j+1].r*SobelH[1][2] + Icolor[i+1][j+1].r*SobelH[2][2];
			IxR[i][j]=abs(totalR)/4;

			//convolution Sobel horizontal
			totalR=Icolor[i-1][j-1].r*SobelV[0][0] + Icolor[i][j-1].r*SobelV[1][0] + Icolor[i+1][j-1].r*SobelV[2][0] + Icolor[i-1][j].r*SobelV[0][1] + Icolor[i][j].r*SobelV[1][1] + Icolor[i+1][j].r*SobelV[2][1] + Icolor[i-1][j+1].r*SobelV[0][2] + Icolor[i+1][j+1].r*SobelV[1][2] + Icolor[i+1][j+1].r*SobelV[2][2];
			IyR[i][j]=abs(totalR)/4;

		}
	}

	for(i=nrl+1;i<nrh;i++)
	{	
		for(j=ncl+1;j<nch;j++)
		{
			//norme gradient
			total=sqrt(Ix[i][j]*Ix[i][j] +Iy[i][j]*Iy[i][j]);
			G[i][j]=total;

			//detection contour
			if(G[i][j]>TRESHOLD)
			{
				C[i][j]=255;
			}else
			{
				C[i][j]=0;
			}
				
		}
	}

	SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"./img_test/cube_conv_reponse_impulsionnelle.pgm");
	SavePGM_bmatrix(Ix,nrl,nrh,ncl,nch,"./img_test/cube_conv_SobelH.pgm");
	SavePGM_bmatrix(Iy,nrl,nrh,ncl,nch,"./img_test/cube_conv_SobelV.pgm");
	SavePGM_bmatrix(G,nrl,nrh,ncl,nch,"./img_test/cube_normeGradient.pgm");
	SavePGM_bmatrix(C,nrl,nrh,ncl,nch,"./img_test/cube_contour.pgm");




	
	free_bmatrix(I,nrl,nrh,ncl,nch);
	free_bmatrix(R,nrl,nrh,ncl,nch);
	free_bmatrix(Ix,nrl,nrh,ncl,nch);
	free_bmatrix(Iy,nrl,nrh,ncl,nch);
	free_bmatrix(G,nrl,nrh,ncl,nch);
	free_bmatrix(C,nrl,nrh,ncl,nch);





	return 1;

}
