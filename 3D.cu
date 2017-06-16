#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <fstream>
#include <cmath>
#include <string>
#include <iostream>

using namespace std;

typedef struct {
   double x,y,z;
} XYZ;

typedef struct {
   XYZ p[8];
   double val[8];
} GRIDCELL;

typedef struct {
   XYZ p[3];         /* Vertices */
   XYZ c;            /* Centroid */
   XYZ n[3];         /* Normal   */
} TRIANGLE;

#define ABS(x) (x < 0 ? -(x) : (x))

// Prototypes
//__global__
//int PolygoniseCube(GRIDCELL,double,TRIANGLE *);
//XYZ VertexInterp(double,XYZ,XYZ,double,double);
/*
#define NX 200
#define NY 160
#define NZ 160
*/
int NX = 100;
int NY = 160;
int NZ = 160;

void fillMatrix(XYZ* a, int n)
{
   int i;
   for (i = 0; i < n; ++i)
   {
        a[i].x = 3;
        a[i].y = 2;
        a[i].z = 5;//rand()%5;
   }
}
__global__
void matrixAdition(XYZ * b, XYZ *a,int n)
{
   	int ij = threadIdx.x + blockDim.x * blockIdx.x;

	if(ij<n)
	{
		b[ij].x = a[ij].x+2;
		b[ij].y = a[ij].y+3;
		b[ij].z = a[ij].z+0;
		//printf("da %d \n" , b[ij].x);
	}
}
void printMatrix(string s, XYZ *a , int tam){
	cout<<s;
	for(int i=0;i<tam;i++)
	{
		cout<<a[i].x<<" "<<a[i].y<<" "<<a[i].z<<" ";
		cout<<endl;
	}
}
void assingMem(int *** data)
{
	int i,j;
	data = (int ***)malloc(NX*sizeof(short int **));
	for (i=0;i<NX;i++)
		data[i] = (int **)malloc(NY*sizeof(short int *));
		for (i=0;i<NX;i++)
			for (j=0;j<NY;j++)
				data[i][j] = (int *)malloc(NZ*sizeof(short int));

}
void readFile(FILE *fptr, const char * namefile , int themin , int themax, int *** data)
{
	int i,j,k,c;
	fprintf(stderr,"Load data ...\n");
	if ((fptr = fopen(namefile,"rb")) == NULL) {
		fprintf(stderr,"Error al leer archivo\n");
		exit(-1);
	}
	for (k=0;k<NZ;k++) {
		for (j=0;j<NY;j++) {
			for (i=0;i<NX;i++) {
				if ((c = fgetc(fptr)) == EOF) {
					fprintf(stderr,"Error en tamaño\n");
					exit(-1);
				}
				data[i][j][k] = c;
				cout<<"leyendo :"<<c<<endl;
				if (c > themax)
					themax = c;
				if (c < themin)
					themin = c;
			}
		}
	}
	fclose(fptr);
	fprintf(stderr,"Rango del volumen: %d -> %d\n",themin,themax);
}

void constructCubes(GRIDCELL * vectGrids, int *** data, int gtam)
{
		int i,j,k;
		//fprintf(stderr,"Construyendo Cubos ...\n");
		int cont=0;
		for (i=0;i<NX-1;i++) {
			//cout<<i<<endl;
			//if (i % (NX/10) == 0)
				//fprintf(stderr,"   Slice %d de %d\n",i,NX);
			for (j=0;j<NY-1;j++) {
				for (k=0;k<NZ-1;k++) {
					GRIDCELL grid;
					grid.p[0].x = i;
					grid.p[0].y = j;
		         	grid.p[0].z = k;
						grid.val[0] = data[i][j][k];
		            grid.p[1].x = i+1;
		            grid.p[1].y = j;
		            grid.p[1].z = k;
						grid.val[1] = data[i+1][j][k];
		            grid.p[2].x = i+1;
		            grid.p[2].y = j+1;
		            grid.p[2].z = k;
						grid.val[2] = data[i+1][j+1][k];
		            grid.p[3].x = i;
		            grid.p[3].y = j+1;
		            grid.p[3].z = k;
						grid.val[3] = data[i][j+1][k];
		            grid.p[4].x = i;
		            grid.p[4].y = j;
		            grid.p[4].z = k+1;
						grid.val[4] = data[i][j][k+1];
		            grid.p[5].x = i+1;
		            grid.p[5].y = j;
		            grid.p[5].z = k+1;
						grid.val[5] = data[i+1][j][k+1];
		            grid.p[6].x = i+1;
		            grid.p[6].y = j+1;
		            grid.p[6].z = k+1;
						grid.val[6] = data[i+1][j+1][k+1];
		            grid.p[7].x = i;
		            grid.p[7].y = j+1;
		            grid.p[7].z = k+1;
						grid.val[7] = data[i][j+1][k+1];
					vectGrids[i+j*NY+k*NY*NZ]=grid;
					cont++;
					//cout<<cont<<endl;
				}
			}
		}

}
__device__
XYZ VertexInterp(double isolevel,XYZ p1,XYZ p2,double valp1,double valp2)
{
   double mu;
   XYZ p;
   if (ABS(isolevel-valp1) < 0.00001)
      return(p1);
   if (ABS(isolevel-valp2) < 0.00001)
      return(p2);
   if (ABS(valp1-valp2) < 0.00001)
      return(p1);
   mu = (isolevel - valp1) / (valp2 - valp1);
   p.x = p1.x + mu * (p2.x - p1.x);
   p.y = p1.y + mu * (p2.y - p1.y);
   p.z = p1.z + mu * (p2.z - p1.z);
   return p;
}

__device__
void copyXYZ(XYZ &a, XYZ &b)
{
	a.x=b.x ; a.y=b.y ; a.z = b.z;
}

__device__
XYZ defect()
{
	XYZ a; 
	a.x=300 ; a.y=300 ; a.z = 300;
	return a;
}


__global__
void coyGRID(GRIDCELL * a, GRIDCELL * b, int x, int y, int z)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = threadIdx.y + blockDim.y * blockIdx.y;
    int k = threadIdx.z + blockDim.z * blockIdx.z;
	
	/*if(i<x && j<y && k<z)
	{
		a[ij].p = b[ij].p;
		a[ij].val = b[ij].val;
	}*/
}

__global__
void copyGRID1(GRIDCELL * a, GRIDCELL * b, int x, int y, int z)
{
	int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = threadIdx.y + blockDim.y * blockIdx.y;
    int k = threadIdx.z + blockDim.z * blockIdx.z;
	
	if(i<x && j<y && k<z)
	{
		for(int w=0;w<8;w++)
		{
			a[i+j*y+k*y*z].p[w] = b[i+j*y+k*y*z].p[w];
			a[i+j*y+k*y*z].val[w] = b[i+j*y+k*y*z].val[w];
		}
	}
}

/*
__global__
void PolygoniseCube(XYZ * vertlist ,GRIDCELL * g ,double iso, int x ,int y , int z)
*/

__global__
void PolygoniseCube(XYZ * vertlist ,GRIDCELL * g ,double iso, int x ,int y , int z)
{
	//printf("g %d \n",iso);
	int i = threadIdx.x + blockDim.x * blockIdx.x;
    int j = threadIdx.y + blockDim.y * blockIdx.y;
    int k = threadIdx.z + blockDim.z * blockIdx.z;
	if(i<x && j<y && k<z)
	{

		//printf("thread %d \n", g[i].p[7].x);
		int cubeindex;
		//int tamVert=12;
		//XYZ vertlist[12];
		int edgeTable[256]={
		0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
		0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
		0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
		0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
		0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
		0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
		0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
		0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
		0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
		0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
		0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
		0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
		0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
		0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
		0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
		0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
		0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
		0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
		0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
		0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
		0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
		0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
		0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
		0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
		0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
		0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
		0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
		0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
		0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
		0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
		0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
		0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };


		//int i,ntri = 0;

		cubeindex = 0;

		if (g[i+j*y+k*y*z].val[0] < iso) cubeindex |= 1;
		if (g[i+j*y+k*y*z].val[1] < iso) cubeindex |= 2;
		if (g[i+j*y+k*y*z].val[2] < iso) cubeindex |= 4;
		if (g[i+j*y+k*y*z].val[3] < iso) cubeindex |= 8;
		if (g[i+j*y+k*y*z].val[4] < iso) cubeindex |= 16;
		if (g[i+j*y+k*y*z].val[5] < iso) cubeindex |= 32;
		if (g[i+j*y+k*y*z].val[6] < iso) cubeindex |= 64;
		if (g[i+j*y+k*y*z].val[7] < iso) cubeindex |= 128;
		

		//XYZ a;
		//a.x=20 ; a.y=50; a,z=0;
		//vertlist[i+j*y+k*y*z+0].x=g[i+j*y+k*y*z].val[6];
		//vertlist[i+j*y+k*y*z+0].y=10;
		//vertlist[i+j*y+k*y*z+0].z=10;*/
	   /* Cube is entirely in/out of the surface */
	   if (edgeTable[cubeindex] == 0)
	      return;
	   /* Find the vertices where the surface intersects the cube */
	   if (edgeTable[cubeindex] & 1) {
	      vertlist[i+j*y+k*y*z+0] = VertexInterp(iso,g[i+j*y+k*y*z].p[0],g[i+j*y+k*y*z].p[1],g[i+j*y+k*y*z].val[0],g[i+j*y+k*y*z].val[1]);
	   }
	   if (edgeTable[cubeindex] & 2) {
	      vertlist[i+j*y+k*y*z+1] = VertexInterp(iso,g[i+j*y+k*y*z].p[1],g[i+j*y+k*y*z].p[2],g[i+j*y+k*y*z].val[1],g[i+j*y+k*y*z].val[2]);
	   }
	   if (edgeTable[cubeindex] & 4) {
	      vertlist[i+j*y+k*y*z+2] = VertexInterp(iso,g[i+j*y+k*y*z].p[2],g[i+j*y+k*y*z].p[3],g[i+j*y+k*y*z].val[2],g[i+j*y+k*y*z].val[3]);
	   }
	   if (edgeTable[cubeindex] & 8) {
	      vertlist[i+j*y+k*y*z+3] = VertexInterp(iso,g[i+j*y+k*y*z].p[3],g[i+j*y+k*y*z].p[0],g[i+j*y+k*y*z].val[3],g[i+j*y+k*y*z].val[0]);
	   }
	   if (edgeTable[cubeindex] & 16) {
	      vertlist[i+j*y+k*y*z+4] = VertexInterp(iso,g[i+j*y+k*y*z].p[4],g[i+j*y+k*y*z].p[5],g[i+j*y+k*y*z].val[4],g[i+j*y+k*y*z].val[5]);
	   }
	   if (edgeTable[cubeindex] & 32) {
	      vertlist[i+j*y+k*y*z+5] = VertexInterp(iso,g[i+j*y+k*y*z].p[5],g[i+j*y+k*y*z].p[6],g[i+j*y+k*y*z].val[5],g[i+j*y+k*y*z].val[6]);
	   }
	   if (edgeTable[cubeindex] & 64) {
	      vertlist[i+j*y+k*y*z+6] = VertexInterp(iso,g[i+j*y+k*y*z].p[6],g[i+j*y+k*y*z].p[7],g[i+j*y+k*y*z].val[6],g[i+j*y+k*y*z].val[7]);
	   }
	   if (edgeTable[cubeindex] & 128) {
	      vertlist[i+j*y+k*y*z+7] = VertexInterp(iso,g[i+j*y+k*y*z].p[7],g[i+j*y+k*y*z].p[4],g[i+j*y+k*y*z].val[7],g[i+j*y+k*y*z].val[4]);
	   }
	   if (edgeTable[cubeindex] & 256) {
	      vertlist[i+j*y+k*y*z+8] = VertexInterp(iso,g[i+j*y+k*y*z].p[0],g[i+j*y+k*y*z].p[4],g[i+j*y+k*y*z].val[0],g[i+j*y+k*y*z].val[4]);
	   }
	   if (edgeTable[cubeindex] & 512) {
	      vertlist[i+j*y+k*y*z+9] = VertexInterp(iso,g[i+j*y+k*y*z].p[1],g[i+j*y+k*y*z].p[5],g[i+j*y+k*y*z].val[1],g[i+j*y+k*y*z].val[5]);
	   }
	   if (edgeTable[cubeindex] & 1024) {
	      vertlist[i+j*y+k*y*z+10] = VertexInterp(iso,g[i+j*y+k*y*z].p[2],g[i+j*y+k*y*z].p[6],g[i+j*y+k*y*z].val[2],g[i+j*y+k*y*z].val[6]);
	   }
	   if (edgeTable[cubeindex] & 2048) {
	      vertlist[i+j*y+k*y*z+11] = VertexInterp(iso,g[i+j*y+k*y*z].p[3],g[i+j*y+k*y*z].p[7],g[i+j*y+k*y*z].val[3],g[i+j*y+k*y*z].val[7]);
	   }
	  // printf("hasta aqui llega \n");
	   
	}
}


void printGrid(string a, GRIDCELL * g, int tam)
{
	cout<<a;
	for(int i =0; i<tam ;i++)
		for(int j=0;j<8;j++)
			//printf("%f  %f  %f \n", g[i].p[j].x ,g[i].p[j].y,g[i].p[j].z);
		      printf("%f \n", g[i].val[j]);		
}

int main(int argc, char *argv[])
{
	int i,j,k,c;
	int ***data;
	FILE *fptr;
	int N= (NX*NY*NZ);
	cout<<N<<endl; //return 1;
	int THREADS_PER_BLOCK =8;
	int themin=255;
	int themax=0;
	int isolevel=80;
	//const char* FILENAME = "mri.raw";
	//assingMem(data);
	//readFile(fptr,FILENAME,themin, themax,data);

	// Malloc the volumetric data, hardwired size!
	data = (int***)malloc(NX*sizeof(int **));
	for (i=0;i<NX;i++)
		data[i] = (int**)malloc(NY*sizeof(int *));
	for (i=0;i<NX;i++)
		for (j=0;j<NY;j++)
			data[i][j] = (int*)malloc(NZ*sizeof(int));

	//cout<<data[199][60][0]<<endl;
	// Open and read the raw data
	fprintf(stderr,"Reading data ...\n");
	if ((fptr = fopen(argv[argc-1],"rb")) == NULL) {
		fprintf(stderr,"File open failed\n");
		exit(-1);
	}
	cout<<"llega"<<endl;
	for (k=0;k<NZ;k++) {
		for (j=0;j<NY;j++) {
			for (i=0;i<NX;i++) {
				if ((c = fgetc(fptr)) == EOF) {
					fprintf(stderr,"Unexpected end of file\n");
					exit(-1);
				}
				data[i][j][k] = c;
				//cout<<i<<" "<<j <<" "<<k <<" data : "<<data[i][j][k]<<endl;
				if (c > themax)
					themax = c;
				if (c < themin)
					themin = c;
			}
		}
	}
	fclose(fptr);
	fprintf(stderr,"Volumetric data range: %d -> %d\n",themin,themax);

	int sizeGRID = N*sizeof(GRIDCELL);
	cout<<"pasa"<<endl;
	int sizeXYZ  = N*12*sizeof(XYZ);

	cout<<"sizeGRID "<<sizeGRID<<endl;
	cout<<"sizeXYZ "<<sizeXYZ<<endl;
	
	//cudaMalloc((void **)&d_b, size);

	GRIDCELL * vectGrids;
	GRIDCELL * d_vectGrid;
	XYZ * d_points;
	XYZ * points;
	points = (XYZ *)malloc(sizeXYZ);
	vectGrids = (GRIDCELL *)malloc(sizeGRID);
	constructCubes(vectGrids,data,N);
	/*
		typedef struct {
		double x,y,z;
		} XYZ;

		typedef struct {
		XYZ p[8];
		double val[8];
		} GRIDCELL;	
	*/	
	XYZ * d_p; double * d_val;

	size_t available, total;
	cudaMemGetInfo(&available, &total);
	cout<<"available:  " << available<<" total:  "<<total <<endl;
	cudaMalloc((void **)&d_vectGrid, sizeGRID);

	/*
	for(int i=0;i<N;i++)
	{
		cudaMalloc((void**)&d_p,8*sizeof(XYZ));
		//cudaMemGetInfo(&available, &total);
		//cout<<"available:  " << available<<" total:  "<<total <<endl;
		cudaMalloc((void**)&d_val,8*sizeof(double));
		//cudaMemGetInfo(&available, &total);
		//cout<<"available:  " << available<<" total:  "<<total <<endl;
		cudaMemcpy(d_p,vectGrids[i].p,8*sizeof(XYZ),cudaMemcpyHostToDevice);
		//for(int w=0;w<8;w++)
		//{
			cout<<vectGrids[i].p[w].y<<endl;
		//}
		//cudaMemGetInfo(&available, &total);
		//cout<<"available:  " << available<<" total:  "<<total <<endl;
		cudaMemcpy(d_val,vectGrids[i].val,8*sizeof(double),cudaMemcpyHostToDevice);
		//cudaMemGetInfo(&available, &total);
		//cout<<"available:  " << available<<" total:  "<<total <<endl
		cudaMemGetInfo(&available, &total);
		//cout<<"available:  " << available<<" total:  "<<total <<endl;
		cudaMemcpy(d_vectGrid[i].val, d_val, 8*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy(d_vectGrid[i].p, d_p, 8*sizeof(XYZ),cudaMemcpyHostToDevice);
 	 }*/
	cudaMemcpy(d_vectGrid,vectGrids, sizeGRID, cudaMemcpyHostToDevice);	
	cout<<"termino de asignar memoria"<<endl;
 	 XYZ * d_a, * d_sal;
 	GRIDCELL * d_res;
 	d_sal=(XYZ *)malloc(sizeXYZ);
 	cudaMalloc((void **)&d_res, sizeGRID);
 	cudaMalloc((void **)&d_a, sizeXYZ);
	cudaMalloc((void **)&d_points, sizeXYZ);
	//cout<<"grid  "<<vectGrids<<endl;
	//cout<<"point "<<points<<endl;

	//fillMatrix(points, N);
	printMatrix("imprimiendo pruevba",points, 10);

	cudaMemcpy(d_points, points, sizeXYZ, cudaMemcpyHostToDevice);
	
	cout<<"grid "<<d_vectGrid<<endl;
	cout<<"pointsssss "<<d_points<<endl;
	//printf("dir %d \n",*d_points);
	cout<<"separa memoria sin problemas"<<endl;
	//printGrid("imprimiendo Grid inicial en Host \n ",vectGrids,N);
	cudaEvent_t start, stop;
	float elapsedTime;
	cudaEventCreate(&start);


	int x = NX; int y = NY ; int z = NZ;
	int blockX= (NX + THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK;
	int blockY= (NY + THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK;
	int blockZ= (NZ + THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK;
	cout<<"blocks : "<<blockX<<" threds:  "<<THREADS_PER_BLOCK<<endl;
	cout<<"blocks : "<<blockY<<" threds:  "<<THREADS_PER_BLOCK<<endl;
	cout<<"blocks : "<<blockZ<<" threds:  "<<THREADS_PER_BLOCK<<endl;
	//int blocks= (10 + THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK;
	/*cout<<"blocks : \n"<<blocks<<"\n threds: \n "<<THREADS_PER_BLOCK<<endl; */

	dim3 dimGrid(blockX, blockY, blockZ);
	dim3 dimBlock(THREADS_PER_BLOCK,THREADS_PER_BLOCK, THREADS_PER_BLOCK);
	cudaEventRecord(start,0);
	isolevel=10;

		//copyGRID1<<<dimGrid,dimBlock>>>(d_res,d_vectGrid,x,y,z);
		PolygoniseCube<<<dimGrid,dimBlock>>>(d_points,d_vectGrid,isolevel,x,y,z);
		//PolygoniseCube<<<blocks,THREADS_PER_BLOCK>>>(d_points,d_vectGrids,isolevel);
		//matrixAdition<<<blocks,THREADS_PER_BLOCK>>>(d_a, d_points,10);
		//matrixAditionCol<<<blocks2,THREADS_PER_BLOCK>>>( d_c, d_a, d_b,N);
	cudaEventCreate(&stop);
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start,stop);
	printf("Elapsed time : %f ms\n" ,elapsedTime);
	cudaMemcpy(points,d_points, sizeXYZ, cudaMemcpyDeviceToHost);
	//GRIDCELL * res;
	//res = (GRIDCELL *)malloc(sizeGRID);
	//cudaMemcpy(res,d_vectGrid, sizeGRID, cudaMemcpyDeviceToHost);
	//printGrid("imprimiendo Grid final despues de la copia \n ",res,N);

	printMatrix("Printing Matrix  A \n",points,N);
	/*/printMatrix("Printing Matrix B \n",b,N);
	//printMatrix("Printing Matrix C \n",c,N);
	*/
	free(points); free(vectGrids);
	cudaFree(d_points); cudaFree(d_vectGrid);
	return 0;
}