#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include <math.h>

#define TRESHOLD 15

void sobel(byte **I, double **Ix, double **Iy, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	long total;
	int SobelH[3][3]={{-1,0,1},{-2,0,2},{-1,0,1} };
	int SobelV[3][3]={{-1,-2,-1},{0,0,0},{1,2,1} };

	for(i=nrl+1;i<nrh;i++)
	{	
		for(j=ncl+1;j<nch;j++)
		{
			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelH[0][0] + I[i][j-1]*SobelH[1][0] + I[i+1][j-1]*SobelH[2][0] + I[i-1][j]*SobelH[0][1] + I[i][j]*SobelH[1][1] + I[i+1][j]*SobelH[2][1] + I[i-1][j+1]*SobelH[0][2] + I[i][j+1]*SobelH[1][2] + I[i+1][j+1]*SobelH[2][2];
			Ix[i][j]=total/4;

			//convolution Sobel horizontal
			total=I[i-1][j-1]*SobelV[0][0] + I[i][j-1]*SobelV[1][0] + I[i+1][j-1]*SobelV[2][0] + I[i-1][j]*SobelV[0][1] + I[i][j]*SobelV[1][1] + I[i+1][j]*SobelV[2][1] + I[i-1][j+1]*SobelV[0][2] + I[i][j+1]*SobelV[1][2] + I[i+1][j+1]*SobelV[2][2];
			Iy[i][j]=total/4;
		}
	}

}

void convert_dmatrix_bmatrix(double **D, byte **B, long nrl, long nrh, long ncl, long nch)
{

	for(int i=nrl+1;i<nrh;i++)
	{	
		for(int j=ncl+1;j<nch;j++)
		{
			B[i][j]=abs(D[i][j]);
		}
	}
}

void norm_gradient(double **SobelX, double **SobelY, double **Sobel, long nrl, long nrh, long ncl, long nch)
{
	for(int i=nrl;i<nrh;i++)
	{	
		for(int j=ncl;j<nch;j++)
		{
			Sobel[i][j]=sqrt(pow(SobelX[i][j],2) + pow(SobelY[i][j],2));
		}
	}	
}

void binarize(byte **I, byte **B, int treshold, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			if (I[i][j] > treshold)
			{
				B[i][j] = 255;
			}
			else
			{
				B[i][j] = 0;
			}
		}
	}
}

void print_histogram(int *H)
{
	int max=0;
	for(int i=0;i<=255;i++)
	{
		if(H[i]>max)
		{
			max=H[i];
		}
	}
	
	for(int i=0;i<=255;i++)
	{
		printf("%d\t|",i);
		for(int j=0;j<=(H[i]*100)/max;j++)
		{
			printf("■");
		}
		printf(" %d \n",H[i]);
	}

}

void convert_rgb8_to_byte(rgb8 **I, byte **B, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	int moy;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			moy= (I[i][j].r +I[i][j].g + I[i][j].b)/3;

				B[i][j] = moy;
			
		}
	}
}


int main()
{

	long nrh,nrl,nch,ncl;
	rgb8 **I;
	byte **grayscale;
	double ** SobelX;
	double **SobelY;
	double **Sobel;
	byte **R;
	byte **Ix;
	byte **Iy;
	byte **binarized;
	int *histogram;
	I=LoadPPM_rgb8matrix("./img_test/bus1.ppm",&nrl,&nrh,&ncl,&nch);

	grayscale=bmatrix(nrl,nrh,ncl,nch);

	convert_rgb8_to_byte(I,grayscale,nrl,nrh,ncl,nch);

	histogram=ivector0(0,255);
	
	SobelX=dmatrix(nrl,nrh,ncl,nch);
	SobelY=dmatrix(nrl,nrh,ncl,nch);
	Sobel=dmatrix(nrl,nrh,ncl,nch);

	Ix=bmatrix(nrl,nrh,ncl,nch);
	Iy=bmatrix(nrl,nrh,ncl,nch);

	binarized=bmatrix(nrl,nrh,ncl,nch);
	R=bmatrix(nrl,nrh,ncl,nch);

	sobel(grayscale,SobelX,SobelY,nrl,nrh,ncl,nch);
	norm_gradient(SobelX, SobelY, Sobel, nrl, nrh, ncl, nch);

	convert_dmatrix_bmatrix(SobelX,Ix,nrl,nrh,ncl,nch);
	convert_dmatrix_bmatrix(SobelY,Iy,nrl,nrh,ncl,nch);
	convert_dmatrix_bmatrix(Sobel,R,nrl,nrh,ncl,nch);

	binarize(R,binarized,TRESHOLD,nrl,nrh,ncl,nch);
	histogram_bmatrix(grayscale,nrl,nrh,ncl,nch,histogram);
	
	
	print_histogram(histogram);
	SavePPM_rgb8matrix(I,nrl,nrh,ncl,nch,"./img_test/imageTest.ppm");
	SavePGM_bmatrix(grayscale,nrl,nrh,ncl,nch,"./img_test/grayscaleTest.pgm");
	SavePGM_bmatrix(Ix,nrl,nrh,ncl,nch,"./img_test/sobelXTest.pgm");
	SavePGM_bmatrix(Iy,nrl,nrh,ncl,nch,"./img_test/sobelYTest.pgm");
	SavePGM_bmatrix(R,nrl,nrh,ncl,nch,"./img_test/SobelTest.pgm");
	SavePGM_bmatrix(binarized,nrl,nrh,ncl,nch,"./img_test/SobelBinary.pgm");

	free_rgb8matrix(I,nrl,nrh,ncl,nch);
	free_bmatrix(grayscale,nrl,nrh,ncl,nch);
	free_bmatrix(Ix,nrl,nrh,ncl,nch);
	free_bmatrix(Iy,nrl,nrh,ncl,nch);
	free_bmatrix(R,nrl,nrh,ncl,nch);
	free_bmatrix(binarized,nrl,nrh,ncl,nch);


	free_dmatrix(SobelX,nrl,nrh,ncl,nch);
	free_dmatrix(SobelY,nrl,nrh,ncl,nch);
	free_dmatrix(Sobel,nrl,nrh,ncl,nch);

	return 0;

}
