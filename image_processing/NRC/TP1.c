#include<stdio.h>
#include<stdlib.h>

#include "def.h"
#include "nrio.h"
#include "nrarith.h"
#include "nralloc.h"
#include <math.h>

#define TRESHOLD 20
#define COLOR_DIFFERENCE 10
#define NB_IMAGES 500
//#define TEST


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

	for(int i=nrl;i<nrh;i++)
	{	
		for(int j=ncl;j<nch;j++)
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

int nb_pixel_contour(byte **B, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	int sum=0;
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			
			sum+= B[i][j];
		}
	}

	sum/=255;
	return sum;
}

void count_rgb(rgb8 **I,int count[3],int colorDiff, long nrl, long nrh, long ncl, long nch)
{
	int i,j;
	rgb8 Iij;
	count[0]=0;
	count[1]=0;
	count[2]=0;
	
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			Iij=I[i][j];
			if((Iij.r>Iij.b+colorDiff) && (Iij.r>Iij.g+colorDiff))
			{
				count[0]++;
			}
			else 
			{
				if((Iij.g>Iij.r+colorDiff) && (Iij.g>Iij.b+colorDiff))
				{
					count[1]++;

				}
				else
				{
					if((Iij.b>Iij.r+colorDiff) && (Iij.b>Iij.g+colorDiff))
					{
						count[2]++;
					}	
				}
			}
				
		}
	}
}

void avg_color(rgb8 **I,rgb8 *avgColor, long nrl, long nrh, long ncl, long nch)
{	
	int i, j;
	int r=0,g=0,b=0;
	
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			r+=I[i][j].r;
			g+=I[i][j].g;
			b+=I[i][j].b;
		}
	}

	r/=(nrh*nch);
	g/=(nrh*nch);
	b/=(nrh*nch);
	
	avgColor->r=r;
	avgColor->g=g;
	avgColor->b=b;
}

void avg_norm_gradient(byte **I,byte *avgNormGradient, long nrl, long nrh, long ncl, long nch)
{	
	int i, j;
	int avg=0;
	
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			avg+=I[i][j];
		}
	}

	avg/=(nrh*nch);

	*avgNormGradient=avg;

}

boolean is_colorful(rgb8 **I,long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	
	for (i = nrl; i < nrh; i++)
	{
		for (j = ncl; j < nch; j++)
		{
			if(I[i][j].r!=I[i][j].b || I[i][j].r!=I[i][j].g || I[i][j].b!=I[i][j].g)
			{
				return TRUE;
			}
		}
	}
	return FALSE;
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
	int nbPixelContour;
	int imgSize;
	int count[3]={0,0,0};
	rgb8 averageColor;
	boolean isColorful;
	byte avgNormGradient;	
	char imgFolderPath[100]="./archive500ppm/";
	char imgPath[100];
	char txtFolderPath[100]="./img_caracteristics/";
	char txtPath[100];
#ifndef TEST

	for (int i=1;i<=NB_IMAGES;i++)
	{
		sprintf(imgPath,"%s%d.%s",imgFolderPath,i,"ppm");

		sprintf(txtPath,"%s%d.%s",txtFolderPath,i,"txt");

		I=LoadPPM_rgb8matrix(imgPath,&nrl,&nrh,&ncl,&nch);

		imgSize=nrh*nch;
		grayscale=bmatrix(nrl,nrh,ncl,nch);

		isColorful=is_colorful(I,nrl,nrh,ncl,nch);
		count_rgb(I,&count,COLOR_DIFFERENCE,nrl,nrh,ncl,nch);

		avg_color(I,&averageColor,nrl,nrh,ncl,nch);
		convert_rgb8_to_byte(I,grayscale,nrl,nrh,ncl,nch);

		histogram=ivector0(0,255);
		
		SobelX=dmatrix(nrl,nrh,ncl,nch);
		SobelY=dmatrix(nrl,nrh,ncl,nch);
		Sobel=dmatrix(nrl,nrh,ncl,nch);

		Ix=bmatrix(nrl,nrh,ncl,nch);
		Iy=bmatrix(nrl,nrh,ncl,nch);
		R=bmatrix(nrl,nrh,ncl,nch);
		binarized=bmatrix(nrl,nrh,ncl,nch);

		sobel(grayscale,SobelX,SobelY,nrl,nrh,ncl,nch);
		norm_gradient(SobelX, SobelY, Sobel, nrl, nrh, ncl, nch);

		convert_dmatrix_bmatrix(SobelX,Ix,nrl,nrh,ncl,nch);
		convert_dmatrix_bmatrix(SobelY,Iy,nrl,nrh,ncl,nch);
		convert_dmatrix_bmatrix(Sobel,R,nrl,nrh,ncl,nch);

		avg_norm_gradient(R,&avgNormGradient,nrl,nrh,ncl,nch);

		binarize(R,binarized,TRESHOLD,nrl,nrh,ncl,nch);
		nbPixelContour=nb_pixel_contour(binarized,nrl,nrh,ncl,nch);

		histogram_bmatrix(grayscale,nrl,nrh,ncl,nch,histogram);

		//création du fichier
			FILE *fp = fopen(txtPath, "w");
			if (fp == NULL)
			{
				printf("Error opening the file %s", txtPath);
			}

			//écrire plus haute valeur d'étiquette
			fprintf(fp, "%d\n", imgSize);
			fprintf(fp, "%d\n", (nbPixelContour));
			fprintf(fp, "%d\n", avgNormGradient);
			fprintf(fp, "%d\n", isColorful);
			fprintf(fp, "%d,%d,%d\n", averageColor.r,averageColor.g,averageColor.b);
			fprintf(fp, "%d,%d,%d\n", count[0],count[1],count[2]);
			
			for(int j=0;j<256;j++)
			{
				fprintf(fp, "%d,",histogram[j]);
			}
			fprintf(fp, "\n");



		free_rgb8matrix(I,nrl,nrh,ncl,nch);
		free_bmatrix(grayscale,nrl,nrh,ncl,nch);
		free_bmatrix(Ix,nrl,nrh,ncl,nch);
		free_bmatrix(Iy,nrl,nrh,ncl,nch);
		free_bmatrix(R,nrl,nrh,ncl,nch);
		free_bmatrix(binarized,nrl,nrh,ncl,nch);

		free_dmatrix(SobelX,nrl,nrh,ncl,nch);
		free_dmatrix(SobelY,nrl,nrh,ncl,nch);
		free_dmatrix(Sobel,nrl,nrh,ncl,nch);
			
		}


#endif
#ifdef TEST
	I=LoadPPM_rgb8matrix("./img_test/bus1.ppm",&nrl,&nrh,&ncl,&nch);

	imgSize=nrh*nch;
	printf("Image size : %d\n",imgSize);
	grayscale=bmatrix(nrl,nrh,ncl,nch);

	isColorful=is_colorful(I,nrl,nrh,ncl,nch);
	printf("Is colorful: %s\n", isColorful ? "true" : "false");
	count_rgb(I,&count,COLOR_DIFFERENCE,nrl,nrh,ncl,nch);
	printf("Count of main colors : r: %d g: %d b: %d\n",count[0],count[1],count[2]);

	avg_color(I,&averageColor,nrl,nrh,ncl,nch);
	printf("Average color:  r: %d g: %d b: %d\n",averageColor.r,averageColor.g,averageColor.b);
	convert_rgb8_to_byte(I,grayscale,nrl,nrh,ncl,nch);

	histogram=ivector0(0,255);
	
	SobelX=dmatrix(nrl,nrh,ncl,nch);
	SobelY=dmatrix(nrl,nrh,ncl,nch);
	Sobel=dmatrix(nrl,nrh,ncl,nch);

	Ix=bmatrix(nrl,nrh,ncl,nch);
	Iy=bmatrix(nrl,nrh,ncl,nch);
	R=bmatrix(nrl,nrh,ncl,nch);
	binarized=bmatrix(nrl,nrh,ncl,nch);

	

	sobel(grayscale,SobelX,SobelY,nrl,nrh,ncl,nch);
	norm_gradient(SobelX, SobelY, Sobel, nrl, nrh, ncl, nch);

	convert_dmatrix_bmatrix(SobelX,Ix,nrl,nrh,ncl,nch);
	convert_dmatrix_bmatrix(SobelY,Iy,nrl,nrh,ncl,nch);
	convert_dmatrix_bmatrix(Sobel,R,nrl,nrh,ncl,nch);

	avg_norm_gradient(R,&avgNormGradient,nrl,nrh,ncl,nch);
	printf("Average norm gradient value: %d\n",avgNormGradient);

	binarize(R,binarized,TRESHOLD,nrl,nrh,ncl,nch);
	nbPixelContour=nb_pixel_contour(binarized,nrl,nrh,ncl,nch);
	printf("Contour pixel number: %d\n",nbPixelContour);

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

#endif

	return 0;

}
