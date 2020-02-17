/* Computes the gravity generated by a set of rectangles over a set of exterior
 * points. This function is called internally by grav2d_GravityComputation.m */
#include<mex.h>
#include<math.h>
/*Macro for extracting element from a matrix stored in column major order*/
#define PCMO(totalRows,row,col) ((col)*(totalRows)+(row))
/*Function code*/
void mexFunction(int nlhs,mxArray* plhs[],int nrhs,const mxArray* prhs[])
{
    /*Inputs and their dimensions*/
    double* rect=NULL;
    double* rho=NULL;
    double* points=NULL;
    int nr=0,cr=0,nrho=0,crho=0,np=0,cp=0;
    /*Output*/
    double* grav=NULL;
    /*Auxiliary variables*/
    int i=0,j=0;
    double xpt=0.0,zpt=0.0,A=0.0,B=0.0,C=0.0,D=0.0,A2=0.0,B2=0.0,C2=0.0,D2=0.0;
    double xm=0.0,xM=0.0,zm=0.0,zM=0.0,p1=0.0,p2=0.0,p3=0.0,p4=0.0,rhor=0.0;
    /*Gravity constant*/
    /* https://physics.nist.gov/cgi-bin/cuu/Value?bg|search_for=universal_in! */
    double G=6.67408e-11;
    /*////////////////////////////////////////////////////////////////////////*/
    /*////////////////////////////////////////////////////////////////////////*/
    /*Rectangles geometry*/
    rect = mxGetPr(prhs[0]);
    nr = mxGetM(prhs[0]);
    cr = mxGetN(prhs[0]);
    /*Density vector*/
    rho = mxGetPr(prhs[1]);
    nrho = mxGetM(prhs[1]);
    crho = mxGetN(prhs[1]);
    /*Points coordinates*/
    points = mxGetPr(prhs[2]);
    np = mxGetM(prhs[2]);
    cp = mxGetN(prhs[2]);
    /*Output vector*/
    plhs[0] = mxCreateNumericMatrix(np,1,mxDOUBLE_CLASS,mxREAL);
    /*Link to the pointer to have access in the code*/
    grav = mxGetPr(plhs[0]);
    for(i=0;i<np;i++)
    {
        grav[PCMO(np,i,0)] = 0.0;
    }
    /*////////////////////////////////////////////////////////////////////////*/
    /*////////////////////////////////////////////////////////////////////////*/
    /*OpenMP parallel*/
#if defined(_OPENMP)
#pragma omp parallel for default(none) \
 shared(np,points,nr,rho,rect,grav,G) \
 private(i,xpt,zpt,j,rhor,xm,xM,zm,zM,A,B,C,D,A2,B2,C2,D2,p1,p2,p3,p4)
#endif
    /*Loop through all points*/
    for(i=0;i<np;i++)
    {
        /*Working point coordinates*/
        xpt = points[PCMO(np,i,0)];
        zpt = points[PCMO(np,i,1)];
        /*Loop through all rectangles*/
        for(j=0;j<nr;j++)
        {
            /*Density*/
            rhor = rho[PCMO(nr,j,0)];
            /*Rectangle extreme coordinates*/
            xm = rect[PCMO(nr,j,0)];
            xM = rect[PCMO(nr,j,1)];
            zm = rect[PCMO(nr,j,2)];
            zM = rect[PCMO(nr,j,3)];
            /*Constants for the attraction formula*/
            A = xpt-xm;
            B = xpt-xM;
            C = zpt-zM;
            D = zpt-zm;
            /*Squared constants*/
            A2 = A*A;
            B2 = B*B;
            C2 = C*C;
            D2 = D*D;
            /*Formula parts initalization*/
            p1 = 0.0;
            p2 = 0.0;
            p3 = 0.0;
            p4 = 0.0;
            /*The computation is done only if the rectangle dimensions are !=0*/
            if((xm!=xM)&&(zm!=zM))
            {
                /*Computing part 1*/
                if(((A!=0.0)||(C!=0.0))&&((A!=0.0)||(D!=0.0)))
                {
                    p1 = A*log((A2+D2)/(A2+C2));
                }
                /*Computing part 2*/
                if(((B!=0.0)||(C!=0.0))&&((B!=0.0)||(D!=0.0)))
                {
                    p2 = B*log((B2+D2)/(B2+C2));
                }
                /*Computing part 3*/
                if(D!=0.0)
                {
                    p3 = 2.0*D*(atan(A/D)-atan(B/D));
                }
                /*Computing part 4*/
                if(C!=0.0)
                {
                    p4 = 2.0*C*(atan(A/C)-atan(B/C));
                }
            }
            /*Add the gravity due to the prism*/
            grav[PCMO(np,i,0)] += G*rhor*(p1-p2+p3-p4);
        }
    } /* --> End of #pragma omp parallel for */
}
/*
Copyright (c) 2018, J.L.G. Pallero, jgpallero@gmail.com,
                    J.L. Fernández Martínez, jlfm@uniovi.es
                    Z. Fernández Muñiz, zulima@uniovi.es
                    Sylvain Bonvalot, sylvain.bonvalot@ird.fr

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.
- Neither the name of the copyright holders nor the names of its contributors
  may be used to endorse or promote products derived from this software without
  specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
