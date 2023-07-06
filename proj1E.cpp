/*
Krishna Patel
CS441 Intro to Computer Graphics
Project 1E
*/

#include <iostream>
#include <fstream>
using namespace std;
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define UNUSED __attribute__((unused))
#define NORMALS
#define M 1000
#define N 1000

double C441(double f)
{
    return ceil(f-0.00001);
}

double F441(double f)
{
    return floor(f+0.00001);
}


typedef struct 
{
    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
} LightingParameters;


LightingParameters 
GetLighting(Camera c)
{
    LightingParameters lp;
    lp.Ka = 0.3;
    lp.Kd = 0.7;
    lp.Ks = 2.8;
    lp.alpha = 50.5;
    lp.lightDir[0] = c.position[0]-c.focus[0];
    lp.lightDir[1] = c.position[1]-c.focus[1];
    lp.lightDir[2] = c.position[2]-c.focus[2];
    double mag = sqrt(lp.lightDir[0]*lp.lightDir[0]
                    + lp.lightDir[1]*lp.lightDir[1]
                    + lp.lightDir[2]*lp.lightDir[2]);
    if (mag > 0)
    {
        lp.lightDir[0] /= mag;
        lp.lightDir[1] /= mag;
        lp.lightDir[2] /= mag;
    }

    return lp;
}

typedef struct
{
    double          A[4][4];     // A[i][j] means row i, column j
} Matrix;


void
PrintMatrix(Matrix m)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        printf("(%.7f %.7f %.7f %.7f)\n", m.A[i][0], m.A[i][1], m.A[i][2], m.A[i][3]);
    }
}

Matrix
ComposeMatrices(Matrix M1, Matrix M2)
{
    Matrix m_out;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            m_out.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                m_out.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }
    return m_out;
}

void 
TransformPoint(Matrix m, const double *ptIn, double *ptOut)
{  
    ptOut[0] = ptIn[0]*m.A[0][0]
             + ptIn[1]*m.A[1][0]
             + ptIn[2]*m.A[2][0]
             + ptIn[3]*m.A[3][0];
    ptOut[1] = ptIn[0]*m.A[0][1]
             + ptIn[1]*m.A[1][1]
             + ptIn[2]*m.A[2][1]
             + ptIn[3]*m.A[3][1];
    ptOut[2] = ptIn[0]*m.A[0][2]
             + ptIn[1]*m.A[1][2]
             + ptIn[2]*m.A[2][2]
             + ptIn[3]*m.A[3][2];
    ptOut[3] = ptIn[0]*m.A[0][3]
             + ptIn[1]*m.A[1][3]
             + ptIn[2]*m.A[2][3]
             + ptIn[3]*m.A[3][3];
}



typedef struct
{
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
} Camera;


double SineParameterize(int curFrame, int nFrames, int ramp)
{  
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {        
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }        
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
} 

Camera       
GetCamera(int frame, int nframes)
{            
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0; 
    c.focus[1] = 0; 
    c.focus[2] = 0;
    c.up[0] = 0;    
    c.up[1] = 1;    
    c.up[2] = 0;    
    return c;       
}

double dot_product(double *a, double *b, int size) {
    double result = 0.0;
    for (int i = 0; i < size; i++) {
        result += a[i] * b[i];
    }
    return result;
}

void cross_product(double a[3], double b[3], double result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

void subtract_vectors(double a[3], double b[3], double result[3]){
	result[0] = a[0] - b[0];
	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
}

void normalize(double *v){
	double length = sqrt((v[0] * v[0]) + (v[1] * v[1]) + (v[2] * v[2]));
	v[0] = v[0]/length; v[1] = v[1]/length; v[2] = v[2]/length;
	}


Matrix GetViewTransform(Camera c)
{
	//printf("View transform:\n");

    Matrix m;
	m.A[0][0] = 1/tan(c.angle/2); m.A[0][1] = 0; m.A[0][2] = 0; m.A[0][3] = 0;
	m.A[1][0] = 0; m.A[1][1] = 1/tan(c.angle/2); m.A[1][2] = 0; m.A[1][3] = 0;
	m.A[2][0] = 0; m.A[2][1] = 0; m.A[2][2] = (c.far + c.near)/(c.far - c.near); m.A[2][3] = -1;
	m.A[3][0] = 0; m.A[3][1] = 0; m.A[3][2] = (2*(c.far * c.near))/(c.far - c.near); m.A[3][3] = 0;
	//PrintMatrix(m);
    return m;
}


Matrix
GetCameraTransform(Camera c)
{   
    Matrix m;
	//printf("Camera transform:\n");
	double zeroVector[3] = {0.0,0.0,0.0};
	double w[3]; 
	w[0] = c.position[0]-c.focus[0];
	w[1] = c.position[1]-c.focus[1];
	w[2] = c.position[2]-c.focus[2];
	normalize(w);
	double v[3];
	double u[3]; cross_product(c.up,w,u);
	normalize(u);
	cross_product(w,u,v);
	//printf("%f, %f, %f\n", v[0], v[1], v[2]);
	normalize(v);
	double t[3]; 
	t[0] = zeroVector[0]-c.position[0]; t[1] = zeroVector[1]-c.position[1]; t[2] = zeroVector[2]-c.position[2]; 
	double ut = dot_product(u,t,3); double vt = dot_product(v,t,3); double wt = dot_product(w,t,3);
	m.A[0][0] = u[0]; m.A[0][1] = v[0]; m.A[0][2]= w[0]; m.A[0][3] = 0;
	m.A[1][0] = u[1]; m.A[1][1] = v[1]; m.A[1][2]= w[1]; m.A[1][3] = 0;
	m.A[2][0] = u[2]; m.A[2][1] = v[2]; m.A[2][2]= w[2]; m.A[2][3] = 0;
	m.A[3][0] = ut; m.A[3][1] = vt; m.A[3][2]= wt; m.A[3][3] = 1;
	//PrintMatrix(m);
    return m;
}

Matrix
GetDeviceTransform()
{   
	//printf("Device transform:\n");

    Matrix m;
	m.A[0][0] = N/2; m.A[0][1] = 0; m.A[0][2] = 0; m.A[0][3] = 0;
	m.A[1][0] = 0; m.A[1][1] = N/2; m.A[1][2] = 0; m.A[1][3] = 0;
	m.A[2][0] = 0; m.A[2][1] = 0; m.A[2][2] = 1; m.A[2][3] = 0;
	m.A[3][0] = N/2; m.A[3][1] = M/2; m.A[3][2] = 0; m.A[3][3] = 1;
	//PrintMatrix(m);
    return m;
}



typedef struct
{
   double         X[3];
   double         Y[3];
   double         Z[3];
   double         color[3][3]; // color[2][0] is for V2, red channel
#ifdef NORMALS
   double         normals[3][3]; // normals[2][0] is for V2, x-component
#endif
} Triangle;

typedef struct
{
   int numTriangles;
   Triangle *triangles;
} TriangleList;

char *
Read3Numbers(char *tmp, double *v1, double *v2, double *v3)
{
    *v1 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v2 = atof(tmp);
    while (*tmp != ' ')
       tmp++;
    tmp++; /* space */
    *v3 = atof(tmp);
    while (*tmp != ' ' && *tmp != '\n')
       tmp++;
    return tmp;
}

TriangleList *
Get3DTriangles()
{
   FILE *f = fopen("ws_tris.txt", "r");
   if (f == NULL)
   {
       fprintf(stderr, "You must place the ws_tris.txt file in the current directory.\n");
       exit(EXIT_FAILURE);
   }
   fseek(f, 0, SEEK_END);
   int numBytes = ftell(f);
   fseek(f, 0, SEEK_SET);
   if (numBytes != 3892295)
   {
       fprintf(stderr, "Your ws_tris.txt file is corrupted.  It should be 3892295 bytes, but you have %d.\n", numBytes);
       exit(EXIT_FAILURE);
   }

   char *buffer = (char *) malloc(numBytes);
   if (buffer == NULL)
   {
       fprintf(stderr, "Unable to allocate enough memory to load file.\n");
       exit(EXIT_FAILURE);
   }
   
   fread(buffer, sizeof(char), numBytes, f);

   char *tmp = buffer;
   int numTriangles = atoi(tmp);
   while (*tmp != '\n')
       tmp++;
   tmp++;
 
   if (numTriangles != 14702)
   {
       fprintf(stderr, "Issue with reading file -- can't establish number of triangles.\n");
       exit(EXIT_FAILURE);
   }

   TriangleList *tl = (TriangleList *) malloc(sizeof(TriangleList));
   tl->numTriangles = numTriangles;
   tl->triangles = (Triangle *) malloc(sizeof(Triangle)*tl->numTriangles);

   for (int i = 0 ; i < tl->numTriangles ; i++)
   {
       for (int j = 0 ; j < 3 ; j++)
       {
           double x, y, z;
           double r, g, b;
           double normals[3];
/*
 * Weird: sscanf has a terrible implementation for large strings.
 * When I did the code below, it did not finish after 45 minutes.
 * Reading up on the topic, it sounds like it is a known issue that
 * sscanf fails here.  Stunningly, fscanf would have been faster.
 *     sscanf(tmp, "(%lf, %lf), (%lf, %lf), (%lf, %lf) = (%d, %d, %d)\n%n",
 *              &x1, &y1, &x2, &y2, &x3, &y3, &r, &g, &b, &numRead);
 *
 *  So, instead, do it all with atof/atoi and advancing through the buffer manually...
 */
           tmp = Read3Numbers(tmp, &x, &y, &z);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, &r, &g, &b);
           tmp += 3; /* space+slash+space */
           tmp = Read3Numbers(tmp, normals+0, normals+1, normals+2);
           tmp++;    /* newline */

           tl->triangles[i].X[j] = x;
           tl->triangles[i].Y[j] = y;
           tl->triangles[i].Z[j] = z;
           tl->triangles[i].color[j][0] = r;
           tl->triangles[i].color[j][1] = g;
           tl->triangles[i].color[j][2] = b;
		   
#ifdef NORMALS
           tl->triangles[i].normals[j][0] = normals[0];
           tl->triangles[i].normals[j][1] = normals[1];
           tl->triangles[i].normals[j][2] = normals[2];
#endif	

       }
	   
   }
   fclose(f);
   free(buffer);
   return tl;
}

    // code to read all of the triangles: 
    //    TriangleList *tl = GetTriangles(0);

    // code to read just the first 100 triangles: 
    //    TriangleList *tl = GetTriangles(0);

typedef struct pixel{
	unsigned char red;
	unsigned char green;
	unsigned char blue;
	double depth;
}Pixel;

typedef struct image{
	Pixel **pixels;
}Image;

void addHeader(int width, int height, ofstream& file){
	file << "P6\n";
	file << width << " " << height << endl;
	file << 255 << endl;
}

Pixel addPixel(unsigned char r, unsigned char g, unsigned char b){
	Pixel p;
	p.red = r; p.green = g; p.blue = b;
	return p;
}
double calculateDepth(double fa, double fb, double b, double a, double y){
	double t= (y-a)/(b-a);
	double result = fa + (t*(fb-fa));
	if(b == a){
		return b;
	}
	return result;
}

double findEndpoint(double topY, double topX, double bottomY, double bottomX, double curY){

	double y = topY - bottomY;
	double x = topX - bottomX;
	double slope = y/x;
	double b = topY - (slope * topX);
	double endPoint = ((curY-b)/slope);
	if(topX == bottomX)
		return topX;
	return endPoint;
}

double findenddepth(double topY, double topX, double bottomY, double bottomX, double curY){

	double y = topY - bottomY;
	double x = topX - bottomX;
	double slope = y/x;
	double b = topY - (slope * topX);
	double endPoint = ((curY-b)/slope);
	if(topX == bottomX)
		return topX;
	if(topY == bottomY)
		return bottomY;	
	return endPoint;
}



void assignPixels(Image *img, int colMin, int colMax, 
int rowMin, UNUSED int rowMax, double *leftEndColor, double *rightEndColor, 
double leftDepth, double rightDepth, double rightEnd, double leftEnd){
	//for(int x = rowMin; x<=rowMax; x++){
		for(int y = colMin; y<=colMax; y++){
			
			double currDepth = calculateDepth(leftDepth, rightDepth, rightEnd, leftEnd,y);
			double R = findenddepth(leftEnd, leftEndColor[0], rightEnd, rightEndColor[0], y);
			double G = findenddepth(leftEnd, leftEndColor[1],  rightEnd,rightEndColor[1], y);
			double B = findenddepth(leftEnd,  leftEndColor[2], rightEnd, rightEndColor[2], y);
			/*printf("====\n");
			printf("leftR: %f, leftG: %f, leftB: %f\n", leftEndColor[0], leftEndColor[1], leftEndColor[2]);
			printf("rightR: %f, rightG: %f, rightB: %f\n", rightEndColor[0], rightEndColor[1], rightEndColor[2]);
			printf("leftDepth = %f, rightDepth = %f, currDepth = %f\n", leftDepth, rightDepth, currDepth);
			printf("red: %f green: %f, blue %f\n", R, G, B);
			printf("row: %d, col: %d\n", rowMin, y);
			printf("====\n");*/
		if((rowMin < 1000 && y < 1000) && (rowMin >= 0 && y >= 0)){
			Pixel p; p.red = (unsigned char)C441(255 * R); p.green = (unsigned char)C441(255 * G); p.blue = (unsigned char)C441(255 * B); p.depth = currDepth;
			if(currDepth >= img->pixels[rowMin][y].depth){
				img->pixels[rowMin][y] = p;
			}
		}
	}
	}
//}

int findMax (double *v, int size){ //returns index
	double temp = 0.0;
	int index = 0;
	for(int i = 0; i < size; i++){
		if(v[i] > temp){
			temp = v[i];
			index = i;
		}
	}
	return index;
}




int findMin(double *v, int size){
	double temp = v[0];
	int index = 0;
	for(int i = 0; i < size; i++){
		if(v[i] < temp){
			index = i;
			temp = v[i];
		}
	}
	return index;
}

int findmidVertex(int v1, int v2){
	if(v1 != 2 and v2 != 2){
		return 2;
	}
	if(v1 != 1 and v2 != 1){
		return 1;
	}else{
		return 0;
	}
}

void swap(double *a, int aIndex, double *b, int bIndex){
	double temp = a[aIndex];
	a[aIndex] = b[bIndex];
	b[bIndex] = temp;
}

void transformTriangle(Triangle * t, Matrix m){
	double v[4] = {t->X[0], t->Y[0], t->Z[0], 1};
	double out[4];

	TransformPoint(m, v, out);
	//printf("out v0: %f, %f, %f\n", out[0], out[1], out[2]);
	t->X[0] = out[0]/out[3]; t->Y[0] = out[1]/out[3]; t->Z[0] = out[2]/out[3];

	v[0] = t->X[1]; v[1] = t->Y[1]; v[2] = t->Z[1]; v[3] = 1;
	TransformPoint(m, v, out);
	t->X[1] = out[0]/out[3]; t->Y[1] = out[1]/out[3]; t->Z[1] = out[2]/out[3];

	v[0] = t->X[2]; v[1] = t->Y[2]; v[2] = t->Z[2]; v[3] = 1;
	TransformPoint(m, v, out);
	t->X[2] = out[0]/out[3]; t->Y[2] = out[1]/out[3]; t->Z[2] = out[2]/out[3];
}


Image *drawImage(Camera c){
	Image *im = (Image *)malloc(sizeof(image));
	im->pixels = (Pixel **)malloc(sizeof(Pixel*)*1000);
	for(int x = 0; x < 1000; x++){
		im->pixels[x] = (Pixel *)malloc(sizeof(Pixel)*1000);
	}
	
 
	for(int y = 0; y < 1000; y++){
		for(int x = 0; x<1000; x++){
			Pixel p; p.blue = 0; p.red = 0;
			p.green = 0; p.depth = -1.0;
			im->pixels[y][x] = p;
		}
	}

	Matrix m = ComposeMatrices(GetCameraTransform(c), GetViewTransform(c));
	m = ComposeMatrices(m, GetDeviceTransform());
	

	TriangleList *tl = Get3DTriangles();
	for(int t = 0; t < tl->numTriangles; t++){
		/*printf("Triangle #: %d\n", t);
		printf("==================\n");
		printf("before v0: %f, %f, %f\n", tl->triangles[0].X[0], tl->triangles[0].Y[0], tl->triangles[0].Z[0]);*/

		transformTriangle(&(tl->triangles[t]), m);
		//printf("transformed v0: %f, %f, %f\n", tl->triangles[0].X[0], tl->triangles[0].Y[0], tl->triangles[0].Z[0]);
		double rightEndColor[3]; double leftEndColor[3];
     	int maxY = findMax(tl->triangles[t].Y, 3); int minY = findMin(tl->triangles[t].Y, 3);
		double rowMax = ((tl->triangles[t].Y[maxY])); double rowMin = (tl->triangles[t].Y[minY]);
		double rightEnd; double leftEnd; double leftDepth; double rightDepth;
		int midVertex = findmidVertex(maxY, minY);
		for (int i = C441(rowMin); i <= F441(rowMax); i++){
			//printf("i: %d\n", i);
			//printf("midpoint Y: %f\n", tl->triangles[t].Y[midVertex]);
			if(i>= 0 && i < 1000){
				 if(i < tl->triangles[t].Y[midVertex]){
					//printf("going up\n");
					rightEnd = findEndpoint(tl->triangles[t].Y[minY], tl->triangles[t].X[minY], tl->triangles[t].Y[maxY], tl->triangles[t].X[maxY], i);
					
					rightDepth = (findenddepth(tl->triangles[t].Y[midVertex], tl->triangles[t].Z[midVertex], tl->triangles[t].Y[minY], tl->triangles[t].Z[minY], i));

					rightEndColor[0] = findenddepth(tl->triangles[t].Y[minY], tl->triangles[t].color[minY][0], tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][0], i);
					rightEndColor[1] = findenddepth(tl->triangles[t].Y[minY], tl->triangles[t].color[minY][1], tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][1], i);
					rightEndColor[2] = findenddepth(tl->triangles[t].Y[minY], tl->triangles[t].color[minY][2], tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][2], i);

					leftEnd = findEndpoint(tl->triangles[t].Y[midVertex], tl->triangles[t].X[midVertex], tl->triangles[t].Y[minY], tl->triangles[t].X[minY], i);

					leftDepth = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].Z[maxY], tl->triangles[t].Y[minY], tl->triangles[t].Z[minY], i);

					leftEndColor[0] = (findenddepth(tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][0], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][0], i));
					leftEndColor[1] = (findenddepth(tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][1], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][1], i));
					leftEndColor[2] = (findenddepth(tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][2], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][2], i));
					
				}else{
					rightEnd = (findEndpoint(tl->triangles[t].Y[maxY], tl->triangles[t].X[maxY], tl->triangles[t].Y[minY], tl->triangles[t].X[minY], i));
					
					rightDepth = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].Z[maxY], tl->triangles[t].Y[midVertex], tl->triangles[t].Z[midVertex], i);

					rightEndColor[0] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][0], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][0], i);
					rightEndColor[1] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][1], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][1], i);
					rightEndColor[2] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][2], tl->triangles[t].Y[minY], tl->triangles[t].color[minY][2], i);
					
					leftEnd = (findEndpoint(tl->triangles[t].Y[maxY], tl->triangles[t].X[maxY], tl->triangles[t].Y[midVertex], tl->triangles[t].X[midVertex], i));

					leftDepth = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].Z[maxY], tl->triangles[t].Y[minY], tl->triangles[t].Z[minY], i);
					leftEndColor[0] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][0], tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][0], i);
					leftEndColor[1] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][1], tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][1], i);
					leftEndColor[2] = findenddepth(tl->triangles[t].Y[maxY], tl->triangles[t].color[maxY][2], tl->triangles[t].Y[midVertex], tl->triangles[t].color[midVertex][2], i);
				}
            //^^uncomment for going down triangles
				if(leftEnd > rightEnd){
					double temp = leftEnd;
					leftEnd = rightEnd;
					rightEnd = temp;
					temp = rightDepth; rightDepth = leftDepth; leftDepth = rightDepth;
					swap(leftEndColor, 0, rightEndColor, 0);
					swap(leftEndColor, 1, rightEndColor, 1);
					swap(leftEndColor, 2, rightEndColor, 2);
				}
				//printf("leftR: %f, leftG: %f, leftB: %f\n", leftEndColor[0], leftEndColor[1], leftEndColor[2]);
			//printf("rightR: %f, rightG: %f, rightB: %f\n", rightEndColor[0], rightEndColor[1], rightEndColor[2]);

				assignPixels(im, C441(leftEnd), F441(rightEnd), i, i, 
				leftEndColor, rightEndColor,
				leftDepth, rightDepth, rightEnd, leftEnd);
			}else{
				continue;
			}
		}
		//maxY = topVertex minY = bottomVertex
		
		


	}
	free(tl->triangles);
	free(tl);
	return im;
}
void writeImage(ofstream &file, Image *i){
	for(int x = 999; x>=0; x--){
		for(int y = 0; y<1000; y++){
				file.write((char *)&(i->pixels[x][y]), 3);
		}
	}
	return;
}

void freeImage(Image *i){
	for(int x = 0; x < 1000; x++){
			free((i->pixels[x]));
		}
	free(i->pixels); free(i);
}

int main(){
	for(int x = 0; x < 1000; x++){
		if(x%250 != 0){
			continue;
		}
		char fileName[BUFSIZ];
		sprintf(fileName, "proj1E_frame%04d.pnm", x);
		ofstream file(fileName);
		addHeader(1000,1000,file);
		Camera c = GetCamera(x, 1000);
		Image *i = drawImage(c);
		writeImage(file, i);
		file.close();
		freeImage(i);
	}
	return 0;
}
