// Omid55

#include "mex.h"


/// My Function Code
void MyCFunction(double *alpha, double **data, size_t m, size_t n, double **net, size_t N, double *initOp, double *alphaBar, double * items)
{
    size_t u;
    size_t i;
    double item;
    double* netu;
    for(u=0;u<N;u++)
    {
        item = exp(alpha[u]);
        for(i=0;i<N;i++)
        {
            if(net[u][i] != 0)    // this means that i is the adjacent of u
            {
                mexPrintf("sallam ");
            }
        }
        
        items[u] = item;
    }
}


/// Gateway Function
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* variable declarations here */
    /// Declarating my variables for code
    double *alpha;
    double **data;
    double **net;
    double *initOp;
    double *alphaBar;
    size_t m;
    size_t n;
    size_t N;
    double *items;
    
    /// Checking number of input and output variables to the function is correct
    /* check for proper number of arguments */
    if(nrhs != 5)  // number of input variables
    {
        mexErrMsgIdAndTxt("MyToolbox:MyArrayProduct:nrhs","Five inputs for the function is required not %d .", nrhs);         //  << CHECK HERE >>
    }
    if(nlhs != 1)  // number of output variables
    {
        mexErrMsgIdAndTxt("MyToolbox:MyArrayProduct:nlhs","One output for the function is required not %d .", nlhs);          //  << CHECK HERE >>
    }
    
// // //     /// Checking type of each variable is correct
// // //     if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0])!=1 )   // first variable should be scaler
// // //     {
// // //         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input multiplier must be a scalar.");
// // //     }
// // //     if(mxGetM(prhs[1])!=1)                                                                                                           // second variable should be a row vector
// // //     {
// // //         mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input must be a row vector.");
// // //     }
    
    
    /* code here */
    /// Loading the input values    
    alpha = mxGetPr(prhs[0]);
    data = mxGetPr(prhs[1]);
    m = mxGetM(prhs[1]);
    n = mxGetN(prhs[1]);
    net = mxGetPr(prhs[2]);
    N = mxGetM(prhs[2]);
    initOp = mxGetPr(prhs[3]);
    alphaBar = mxGetPr(prhs[4]);
    
    /// Creating output values
    plhs[0] = mxCreateDoubleMatrix(1,N,mxREAL);
    items = mxGetPr(plhs[0]);
    
    /// Call my function and calculate the output
    MyCFunction(alpha,data,m,n,net,N,initOp,alphaBar,items);
}


 