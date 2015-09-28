#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <string.h>
#include <math.h>

//FUNCTON DECLARTIONS
/******************************************************************************/
int main ( int argc, char **argv);
double run(int version);
double Average(double ary[]);
int check();
void assemble (int version);
double ff ( double x );
void geometry (int version);
void init ();
void output ();
void phi ( int il, double x, double *phii, double *phiix, double xleft, 
  double xrite );
double pp ( double x );
void prsys ();
double qq ( double x );
void solve ();
void timestamp ( void );

//GLOBALS
/******************************************************************************/
    long int NSUB;
    # define NL 20

	FILE *fp_out;
	FILE *fp_sol;

    double *adiag;
    double *aleft;
    double *arite;
    double *f;
    double *h;
    int ibc;
    int *indx;
    int *node;
    int nquad;
    int nu;
    double ul;
    double ur;
    double xl;
    double *xn;
    double *xquad;
    double xr;

/******************************************************************************/
int main(int argc, char **argv)
{
	int i;
	int threads;
	int trials;

	//get NSUB, threads and trails from argument
	if(argc != 4){
		printf("Invalid number of arguments.\n");
		printf("Usage: ./fem1d [SUB_SIZE] [NUM_THREADS] [TRIALS]\n");
		exit(EXIT_FAILURE);
	} 
	else if((NSUB = atoi(argv[1])) == 0) {
		printf("Invalid subdivison size.\n");
		printf("Usage: ./fem1d [SUB_SIZE] [NUM_THREADS] [TRIALS]\n");
		exit(EXIT_FAILURE);
	} else if ((threads = atoi(argv[2])) == 0){
		printf("Invalid number of threads.\n");
		printf("Usage: ./fem1d [SUB_SIZE] [NUM_THREADS] [TRIALS]\n");
		exit(EXIT_FAILURE);
	} else if ((trials = atoi(argv[3])) == 0){
		printf("Invalid number of trails.\n");
		printf("Usage: ./fem1d [SUB_SIZE] [NUM_THREADS] [TRIALS]\n");
		exit(EXIT_FAILURE);
	}

	//set number of threads
	omp_set_num_threads(threads);

	//wipe solutions files
	fp_sol = fopen("old_sol.txt","w");
	fclose(fp_sol);
	fp_sol = fopen("new_sol.txt","w");
	fclose(fp_sol);

	//time spent on execution for each trail
	double time_spent[trials];

	//RUN OLD VERSION TRAILS
	printf("**************************** RUNNING OLD VERSION ******************************\n");
	for(i=0;i<trials;i++){
		printf("Running old version, trail %d...\n",i+1);
		time_spent[i] = run(0);
		printf("Old version, trail %d succesfully complete.\n",i+1);
		printf("Time taken: %fsec\n",time_spent[i]);
	}
	//store old time for later display
	double old_time = Average(time_spent);
	printf("**************************** RUNNING NEW VERSION ******************************\n");

	//RUN NEW VERSION TRAILS
	for(i=0;i<trials;i++){
		printf("Running new version, trail %d...\n",i+1);
		time_spent[i] = run(1);
		printf("New version, trail %d succesfully complete.\n", i+1); 
		printf("Time taken: %fsec\n",time_spent[i]);
	}

	printf("********************************** RESULTS ************************************\n");

	double new_time = Average(time_spent);
	printf("New version Average Time: %fsec\n",new_time);
	printf("Old version Average Time: %fsec\n",old_time);
	printf("Speed up: %fsec\n", (old_time - new_time));

	//CHECK FOR CORRCTNESS
	if(check()==0){
		printf("CORRECT: Outputs are identical.\n");
	} else {
		printf("INCORRECT: Outputs are not identical.\n");
	}

	return 0;
}
/******************************************************************************/
/**
*			Computes and returns average of an array of doubles
*
*/
double Average(double ary[]){
	double ave = 0;
	int i;
	int length = (int)(sizeof(ary)/sizeof(double));
	for(i=0;i < length;i++){
		ave += ary[i];
	}
	return ave / (double)length;
}
/******************************************************************************/

	//TODO take print statments out of parallel for loops

/******************************************************************************/
/*
  Purpose:

    RUN is the main program for FEM1D.

  Discussion:

    FEM1D solves a one dimensional ODE using the finite element method.

    The differential equation solved is

      - d/dX (P dU/dX) + Q U  =  F

    The finite-element method uses piecewise linear basis functions.

    Here U is an unknown scalar function of X defined on the
    interval [XL,XR], and P, Q and F are given functions of X.

    The values of U or U' at XL and XR are also specified.

    The interval [XL,XR] is "meshed" with NSUB+1 points,

    XN(0) = XL, XN(1)=XL+H, XN(2)=XL+2*H, ..., XN(NSUB)=XR.

    This creates NSUB subintervals, with interval number 1
    having endpoints XN(0) and XN(1), and so on up to interval
    NSUB, which has endpoints XN(NSUB-1) and XN(NSUB).

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    double ADIAG(NU), the "diagonal" coefficients.  That is, ADIAG(I) is
    the coefficient of the I-th unknown in the I-th equation.

    double ALEFT(NU), the "left hand" coefficients.  That is, ALEFT(I) 
    is the coefficient of the (I-1)-th unknown in the I-th equation.
    There is no value in ALEFT(1), since the first equation
    does not refer to a "0-th" unknown.

    double ARITE(NU).
    ARITE(I) is the "right hand" coefficient of the I-th
    equation in the linear system.  ARITE(I) is the coefficient
    of the (I+1)-th unknown in the I-th equation.  There is
    no value in ARITE(NU) because the NU-th equation does not
    refer to an "NU+1"-th unknown.

    double F(NU).
    ASSEMBLE stores into F the right hand side of the linear
    equations.
    SOLVE replaces those values of F by the solution of the
    linear equations.

    double H(NSUB)
    H(I) is the length of subinterval I.  This code uses
    equal spacing for all the subintervals.

    int IBC.
    IBC declares what the boundary conditions are.
    1, at the left endpoint, U has the value UL,
       at the right endpoint, U' has the value UR.
    2, at the left endpoint, U' has the value UL,
       at the right endpoint, U has the value UR.
    3, at the left endpoint, U has the value UL,
       and at the right endpoint, U has the value UR.
    4, at the left endpoint, U' has the value UL,
       at the right endpoint U' has the value UR.

    int INDX[NSUB+1].
    For a node I, INDX(I) is the index of the unknown
    associated with node I.
    If INDX(I) is equal to -1, then no unknown is associated
    with the node, because a boundary condition fixing the
    value of U has been applied at the node instead.
    Unknowns are numbered beginning with 1.
    If IBC is 2 or 4, then there is an unknown value of U
    at node 0, which will be unknown number 1.  Otherwise,
    unknown number 1 will be associated with node 1.
    If IBC is 1 or 4, then there is an unknown value of U
    at node NSUB, which will be unknown NSUB or NSUB+1,
    depending on whether there was an unknown at node 0.

    int NL.
    The number of basis functions used in a single
    subinterval.  (NL-1) is the degree of the polynomials
    used.  For this code, NL is fixed at 2, meaning that
    piecewise linear functions are used as the basis.

    int NODE[NL*NSUB].
    For each subinterval I:
    NODE[0+I*2] is the number of the left node, and
    NODE[1+I*2] is the number of the right node.

    int NQUAD.
    The number of quadrature points used in a subinterval.
    This code uses NQUAD = 1.

    int NSUB.
    The number of subintervals into which the interval [XL,XR] is broken.

    int NU.
    NU is the number of unknowns in the linear system.
    Depending on the value of IBC, there will be NSUB-1,
    NSUB, or NSUB+1 unknown values, which are the coefficients
    of basis functions.

    double UL.
    If IBC is 1 or 3, UL is the value that U is required
    to have at X = XL.
    If IBC is 2 or 4, UL is the value that U' is required
    to have at X = XL.

    double UR.
    If IBC is 2 or 3, UR is the value that U is required
    to have at X = XR.
    If IBC is 1 or 4, UR is the value that U' is required
    to have at X = XR.

    double XL.
    XL is the left endpoint of the interval over which the
    differential equation is being solved.

    double XN(0:NSUB).
    XN(I) is the location of the I-th node.  XN(0) is XL,
    and XN(NSUB) is XR.

    double XQUAD(NSUB)
    XQUAD(I) is the location of the single quadrature point
    in interval I.

    double XR.
    XR is the right endpoint of the interval over which the
    differential equation is being solved.
*/
  double run(int version){

   if(version == 0){
      if((fp_out = fopen("old_out.txt", "w+")) == NULL || 
      	(fp_sol = fopen("old_sol.txt", "a")) == NULL){
      		printf("Old Version files not found.\n");
      		exit(EXIT_FAILURE);
      }
    }
    else if (version == 1){
      if((fp_out = fopen("new_out.txt", "w+")) == NULL || 
      	(fp_sol = fopen("new_sol.txt", "a")) == NULL){
      		printf("New Version files not found.\n");
      		exit(EXIT_FAILURE);
      }
    }

  	//START TIMER//
  	double begin, end, time_spent;
  	begin = omp_get_wtime();

  	//Allocate array memory
  	adiag = (double *)malloc(sizeof(double)*(double)(NSUB+1));
    aleft = (double *)malloc(sizeof(double)*(double)(NSUB+1));
    arite = (double *)malloc(sizeof(double)*(double)(NSUB+1));
    f = (double *)malloc(sizeof(double)*(double)(NSUB+1));
    h = (double *)malloc(sizeof(double)*(double)(NSUB));
    indx = (int *)malloc(sizeof(int)*(int)(NSUB+1));
    node = (int *)malloc(sizeof(int)*((int)NL*(int)NSUB));
    xn = (double *)malloc(sizeof(double)*(double)(NSUB+1));
    xquad = (double *)malloc(sizeof(double)*(double)(NSUB));  	

    timestamp ();

    fprintf (fp_out, "\n" );
    fprintf (fp_out, "FEM1D\n" );
    fprintf (fp_out, "  C version\n" );
    fprintf (fp_out, "\n" );
    fprintf (fp_out, "  Solve the two-point boundary value problem\n" );
    fprintf (fp_out, "\n" );
    fprintf (fp_out, "  - d/dX (P dU/dX) + Q U  =  F\n" );
    fprintf (fp_out, "\n" );
    fprintf (fp_out, "  on the interval [XL,XR], specifying\n" );
    fprintf (fp_out,"  the value of U or U' at each end.\n" );
    fprintf (fp_out, "\n" );
    fprintf (fp_out,"  The interval [XL,XR] is broken into NSUB = %ld subintervals\n", NSUB );
    fprintf (fp_out, "  Number of basis functions per element is NL = %d\n", NL );

  /*
    Initialize the data.
  */
    init ();
  /*
    Compute the geometric quantities.
  */
    geometry (version);
  /*
    Assemble the linear system.
  */

    assemble (version);
  /*
    Print out the linear system.
  */
    prsys ();
  /*
    Solve the linear system.
  */
    solve ();

  /*
    Print out the solution.
  */
    output ();
  /*
    Terminate.
  */
    fprintf (fp_out, "\n" );
    fprintf (fp_out,"FEM1D:\n" );
    fprintf (fp_out, "  Normal end of execution.\n" );

    fprintf ( fp_out,"\n" );
    timestamp ( );

    //CLOSE STREAMS
    fclose(fp_out);
    fclose(fp_sol);

    //FREE MEMORY
    free(adiag); 
    free(aleft);
    free(arite); 
    free(f); 
    free(h); 
    free(indx); 
    free(node); 
    free(xn); 
    free(xquad);

    //END TIMER//
    end = omp_get_wtime();
    time_spent = end - begin;


    return time_spent;
}

/******************************************************************************/
int check(){

	FILE *fp1 = fopen("old_sol.txt","r");
	FILE *fp2 = fopen("new_sol.txt","r");

	int ch1, ch2;

   if (fp1 == NULL) {
      printf("Cannot open for reading ld_out.txt");
      exit(1);
   } else if (fp2 == NULL) {
      printf("Cannot open for reading new_out.txt");
      exit(1);
   } else {
      ch1 = getc(fp1);
      ch2 = getc(fp2);
      //printf("%d, %d\n",ch1,ch2);
 
      while ((ch1 != EOF) && (ch2 != EOF) && (ch1 == ch2)) {
         ch1 = getc(fp1);
         ch2 = getc(fp2);
      }
 
      if (ch1 == ch2){
        fclose(fp1);
        fclose(fp2);
        return 0;
      } else if (ch1 != ch2) {
        fclose(fp1);
        fclose(fp2);
        return 1;
      }
   }
   return 1;
}



/******************************************************************************/

/*
  Purpose:

    ASSEMBLE assembles the matrix and right-hand-side of the linear system.

  Discussion:

    The linear system has the form:

      K * C = F

    that is to be solved for the coefficients C.

    Numerical integration is used to compute the entries of K and F.

    Note that a 1 point quadrature rule, which is sometimes used to
    assemble the matrix and right hand side, is just barely accurate
    enough for simple problems.  If you want better results, you
    should use a quadrature rule that is more accurate.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Output, double ADIAG(NU), the "diagonal" coefficients.  That is, 
    ADIAG(I) is the coefficient of the I-th unknown in the I-th equation.

    Output, double ALEFT(NU), the "left hand" coefficients.  That is, 
    ALEFT(I) is the coefficient of the (I-1)-th unknown in the I-th equation.
    There is no value in ALEFT(1), since the first equation
    does not refer to a "0-th" unknown.

    Output, double ARITE(NU).
    ARITE(I) is the "right hand" coefficient of the I-th
    equation in the linear system.  ARITE(I) is the coefficient
    of the (I+1)-th unknown in the I-th equation.  There is
    no value in ARITE(NU) because the NU-th equation does not
    refer to an "NU+1"-th unknown.

    Output, double F(NU).
    ASSEMBLE stores into F the right hand side of the linear
    equations.
    SOLVE replaces those values of F by the solution of the
    linear equations.

    Input, double H(NSUB)
    H(I) is the length of subinterval I.  This code uses
    equal spacing for all the subintervals.

    Input, int INDX[NSUB+1].
    For a node I, INDX(I) is the index of the unknown
    associated with node I.
    If INDX(I) is equal to -1, then no unknown is associated
    with the node, because a boundary condition fixing the
    value of U has been applied at the node instead.
    Unknowns are numbered beginning with 1.
    If IBC is 2 or 4, then there is an unknown value of U
    at node 0, which will be unknown number 1.  Otherwise,
    unknown number 1 will be associated with node 1.
    If IBC is 1 or 4, then there is an unknown value of U
    at node NSUB, which will be unknown NSUB or NSUB+1,
    depending on whether there was an unknown at node 0.

    Input, int NL.
    The number of basis functions used in a single
    subinterval.  (NL-1) is the degree of the polynomials
    used.  For this code, NL is fixed at 2, meaning that
    piecewise linear functions are used as the basis.

    Input, int NODE[NL*NSUB].
    For each subinterval I:
    NODE[0+I*2] is the number of the left node, and
    NODE[1+I*2] is the number of the right node.

    Input, int NU.
    NU is the number of unknowns in the linear system.
    Depending on the value of IBC, there will be NSUB-1,
    NSUB, or NSUB+1 unknown values, which are the coefficients
    of basis functions.

    Input, int NQUAD.
    The number of quadrature points used in a subinterval.
    This code uses NQUAD = 1.

    Input, int NSUB.
    The number of subintervals into which the interval [XL,XR] is broken.

    Input, double UL.
    If IBC is 1 or 3, UL is the value that U is required
    to have at X = XL.
    If IBC is 2 or 4, UL is the value that U' is required
    to have at X = XL.

    Input, double UR.
    If IBC is 2 or 3, UR is the value that U is required
    to have at X = XR.
    If IBC is 1 or 4, UR is the value that U' is required
    to have at X = XR.

    Input, double XL.
    XL is the left endpoint of the interval over which the
    differential equation is being solved.

    Input, double XR.
    XR is the right endpoint of the interval over which the
    differential equation is being solved.
*/
 void assemble (int version){

    if(version == 0){

    double aij;
    double he;
    int i;
    int ie;
    int ig;
    int il;
    int iq;
    int iu;
    int jg;
    int jl;
    int ju;
    double phii;
    double phiix;
    double phij;
    double phijx;
    double x;
    double xleft;
    double xquade;
    double xrite;
	  /*
	    Zero out the arrays that hold the coefficients of the matrix
	    and the right hand side.
	  */
	    for ( i = 0; i < nu; i++ )
	    {
	      f[i] = 0.0;
	    }
	    for ( i = 0; i < nu; i++ )
	    {
	      adiag[i] = 0.0;
	    }
	    for ( i = 0; i < nu; i++ )
	    {
	      aleft[i] = 0.0;
	    }
	    for ( i = 0; i < nu; i++ )
	    {
	      arite[i] = 0.0;
	    }

	  /*
	    For interval number IE,
	  */
	    for ( ie = 0; ie < NSUB; ie++ )
	    {
	      he = h[ie];
	      xleft = xn[node[0+ie*2]];
	      xrite = xn[node[1+ie*2]];
	  /*
	    consider each quadrature point IQ,
	  */
	      for ( iq = 0; iq < nquad; iq++ )
	      {
	        xquade = xquad[ie];
	  /*
	    and evaluate the integrals associated with the basis functions
	    for the left, and for the right nodes.
	  */
	        for ( il = 1; il <= NL; il++ )
	        {
	          ig = node[il-1+ie*2];
	          iu = indx[ig] - 1;

	          if ( 0 <= iu )
	          {
	            phi ( il, xquade, &phii, &phiix, xleft, xrite );
	            f[iu] = f[iu] + he * ff ( xquade ) * phii;
	  /*
	    Take care of boundary nodes at which U' was specified.
	  */
	            if ( ig == 0 )
	            {
	              x = 0.0;
	              f[iu] = f[iu] - pp ( x ) * ul;
	            }
	            else if ( ig == NSUB )
	            {
	              x = 1.0;
	              f[iu] = f[iu] + pp ( x ) * ur;
	            }



	  /*
	    Evaluate the integrals that take a product of the basis
	    function times itself, or times the other basis function
	    that is nonzero in this interval.
	  */
	            for ( jl = 1; jl <= NL; jl++ )
	            {

	              jg = node[jl-1+ie*2];


	              ju = indx[jg] - 1;

	              phi ( jl, xquade, &phij, &phijx, xleft, xrite );

	              aij = he * ( pp ( xquade ) * phiix * phijx 
	                         + qq ( xquade ) * phii  * phij   );
	  /*
	    If there is no variable associated with the node, then it's
	    a specified boundary value, so we multiply the coefficient
	    times the specified boundary value and subtract it from the
	    right hand side.
	  */

	              if ( ju < 0 )
	              {
	                if ( jg == 0 )
	                {
	                  f[iu] = f[iu] - aij * ul;
	                }
	                else if ( jg == NSUB )
	                {               
	                  f[iu] = f[iu] - aij * ur;
	                }
	              }

	  /*
	    Otherwise, we add the coefficient we've just computed to the
	    diagonal, or left or right entries of row IU of the matrix.
	  */
	              else
	              {
	                if ( iu == ju )
	                {
	                  adiag[iu] = adiag[iu] + aij;
	                }
	                else if ( ju < iu )
	                {
	                  aleft[iu] = aleft[iu] + aij;
	                }
	                else
	                {
	                  arite[iu] = arite[iu] + aij;
	                }
	              }
	            }
	           // printf("%d\n",jl );
	          }
	        }
	      }
	    }

	    return;
//****************************** PARALLELIZED CODE **********************************
	} else {

    double aij;
    double he;
    int i;
    int ie;
    int ig;
    int il;
    int iq;
    int iu;
    int jg;
    int jl;
    int ju;
    double phii;
    double phiix;
    double phij;
    double phijx;
    double x;
    double xleft;
    double xquade;
    double xrite;
      /*
        Zero out the arrays that hold the coefficients of the matrix
        and the right hand side.
      */
        #pragma omp parallel for
        for ( i = 0; i < nu; i++ )
        {
          f[i] = 0.0;
          adiag[i] = 0.0;
          aleft[i] = 0.0;
          arite[i] = 0.0;
        }

      /*
        For interval number IE,
      */
        for ( ie = 0; ie < NSUB; ie++ )
        {
          he = h[ie];
          xleft = xn[node[0+ie*2]];
          xrite = xn[node[1+ie*2]];
      /*
        consider each quadrature point IQ,
      */
          for ( iq = 0; iq < nquad; iq++ )
          {
            xquade = xquad[ie];
      /*
        and evaluate the integrals associated with the basis functions
        for the left, and for the right nodes.
      */
            for ( il = 1; il <= NL; il++ )
            {
              ig = node[il-1+ie*2];
              iu = indx[ig] - 1;

              if ( 0 <= iu )
              {
                phi ( il, xquade, &phii, &phiix, xleft, xrite );
                f[iu] = f[iu] + he * ff ( xquade ) * phii;
      /*
        Take care of boundary nodes at which U' was specified.
      */
                if ( ig == 0 )
                {
                  x = 0.0;
                  f[iu] = f[iu] - pp ( x ) * ul;
                }
                else if ( ig == NSUB )
                {
                  x = 1.0;
                  f[iu] = f[iu] + pp ( x ) * ur;
                }



      /*
        Evaluate the integrals that take a product of the basis
        function times itself, or times the other basis function
        that is nonzero in this interval.
      */
                for ( jl = 1; jl <= NL; jl++ )
                {

                  jg = node[jl-1+ie*2];


                  ju = indx[jg] - 1;

                  phi ( jl, xquade, &phij, &phijx, xleft, xrite );

                  aij = he * ( pp ( xquade ) * phiix * phijx 
                             + qq ( xquade ) * phii  * phij   );
      /*
        If there is no variable associated with the node, then it's
        a specified boundary value, so we multiply the coefficient
        times the specified boundary value and subtract it from the
        right hand side.
      */

                  if ( ju < 0 )
                  {
                    if ( jg == 0 )
                    {
                      f[iu] = f[iu] - aij * ul;
                    }
                    else if ( jg == NSUB )
                    {               
                      f[iu] = f[iu] - aij * ur;
                    }
                  }

      /*
        Otherwise, we add the coefficient we've just computed to the
        diagonal, or left or right entries of row IU of the matrix.
      */
                  else
                  {
                    if ( iu == ju )
                    {
                      adiag[iu] = adiag[iu] + aij;
                    }
                    else if ( ju < iu )
                    {
                      aleft[iu] = aleft[iu] + aij;
                    }
                    else
                    {
                      arite[iu] = arite[iu] + aij;
                    }
                  }
                }
               // printf("%d\n",jl );
              }
            }
          }
        }

        return;
	}
}
  
/******************************************************************************/

/*
  Purpose:

    FF evaluates the right hand side function.

  Discussion:

    This routine evaluates the function F(X) in the differential equation.

      -d/dx (p du/dx) + q u  =  f

    at the point X.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double FF, the value of the function.
*/
  double ff ( double x ){
    double value;

    value = 0.0;

    return value;
  }
/******************************************************************************/

/*
  Purpose: 

    GEOMETRY sets up the geometry for the interval [XL,XR].

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Output, double H(NSUB)
    H(I) is the length of subinterval I.  This code uses
    equal spacing for all the subintervals.

    Input, int IBC.
    IBC declares what the boundary conditions are.
    1, at the left endpoint, U has the value UL,
       at the right endpoint, U' has the value UR.
    2, at the left endpoint, U' has the value UL,
       at the right endpoint, U has the value UR.
    3, at the left endpoint, U has the value UL,
       and at the right endpoint, U has the value UR.
    4, at the left endpoint, U' has the value UL,
       at the right endpoint U' has the value UR.

    Output, int INDX[NSUB+1].
    For a node I, INDX(I) is the index of the unknown
    associated with node I.
    If INDX(I) is equal to -1, then no unknown is associated
    with the node, because a boundary condition fixing the
    value of U has been applied at the node instead.
    Unknowns are numbered beginning with 1.
    If IBC is 2 or 4, then there is an unknown value of U
    at node 0, which will be unknown number 1.  Otherwise,
    unknown number 1 will be associated with node 1.
    If IBC is 1 or 4, then there is an unknown value of U
    at node NSUB, which will be unknown NSUB or NSUB+1,
    depending on whether there was an unknown at node 0.

    Input, int NL.
    The number of basis functions used in a single
    subinterval.  (NL-1) is the degree of the polynomials
    used.  For this code, NL is fixed at 2, meaning that
    piecewise linear functions are used as the basis.

    Output, int NODE[NL*NSUB].
    For each subinterval I:
    NODE[0+I*2] is the number of the left node, and
    NODE[1+I*2] is the number of the right node.

    Input, int NSUB.
    The number of subintervals into which the interval [XL,XR] is broken.

    Output, int *NU.
    NU is the number of unknowns in the linear system.
    Depending on the value of IBC, there will be NSUB-1,
    NSUB, or NSUB+1 unknown values, which are the coefficients
    of basis functions.

    Input, double XL.
    XL is the left endpoint of the interval over which the
    differential equation is being solved.

    Output, double XN(0:NSUB).
    XN(I) is the location of the I-th node.  XN(0) is XL,
    and XN(NSUB) is XR.

    Output, double XQUAD(NSUB)
    XQUAD(I) is the location of the single quadrature point
    in interval I.

    Input, double XR.
    XR is the right endpoint of the interval over which the
    differential equation is being solved.
*/
  void geometry (int version){
  if(version == 0){
	  long int i;
	  /*
	    Set the value of XN, the locations of the nodes.
	  */
	    fprintf ( fp_out,"\n" );
	    fprintf (fp_out, "  Node      Location\n" );
	    fprintf (fp_out, "\n" );

	    for ( i = 0; i <= NSUB; i++ )
	    {
	      xn[i]  =  ( ( double ) ( NSUB - i ) * xl 
	                + ( double )          i   * xr ) 
	                / ( double ) ( NSUB );
	      fprintf (fp_out, "  %8ld  %14f \n", i, xn[i] );
	    }
	  /*
	    Set the lengths of each subinterval.
	  */
	    fprintf (fp_out, "\n" );
	    fprintf (fp_out, "Subint    Length\n" );
	    fprintf ( fp_out,"\n" );

	    for ( i = 0; i < NSUB; i++ )
	    {
	      h[i] = xn[i+1] - xn[i];
	      fprintf (fp_out, "  %8ld  %14f\n", i+1, h[i] );
	    }
	  /*
	    Set the quadrature points, each of which is the midpoint
	    of its subinterval.
	  */
	    fprintf (fp_out, "\n" );
	    fprintf (fp_out, "Subint    Quadrature point\n" );
	    fprintf ( fp_out,"\n" );

	    for ( i = 0; i < NSUB; i++ )
	    {
	      xquad[i] = 0.5 * ( xn[i] + xn[i+1] );
	      fprintf ( fp_out,"  %8ld  %14f\n", i+1, xquad[i] );
	    }
	  /*
	    Set the value of NODE, which records, for each interval,
	    the node numbers at the left and right.
	  */
	    fprintf ( fp_out,"\n" );
	    fprintf ( fp_out,"Subint  Left Node  Right Node\n" );
	    fprintf (fp_out, "\n" );

	    for ( i = 0; i < NSUB; i++ )
	    {
	      node[0+i*2] = i;
	      node[1+i*2] = i + 1;
	      fprintf (fp_out, "  %8ld  %8d  %8d\n", i+1, node[0+i*2], node[1+i*2] );
	    }
	  /*
	    Starting with node 0, see if an unknown is associated with
	    the node.  If so, give it an index.
	  */
	    nu = 0;
	  /*
	    Handle first node.
	  */
	    i = 0;
	    if ( ibc == 1 || ibc == 3 )
	    {
	      indx[i] = -1;
	    }
	    else
	    {
	      nu = nu + 1;
	      indx[i] = nu;
	    }
	  /*
	    Handle nodes 1 through nsub-1
	  */

	    for ( i = 1; i < NSUB; i++ )
	    {
	      nu = nu + 1;
	      indx[i] = nu;
	    }
	  /*
	    Handle the last node.
	  /*/
	    i = NSUB;

	    if ( ibc == 2 || ibc == 3 )
	    {
	      indx[i] = -1;
	    }
	    else
	    {
	      nu = nu + 1;
	      indx[i] = nu;
	    }

	    fprintf ( fp_out,"\n" );
	    fprintf ( fp_out,"  Number of unknowns NU = %8d\n", nu );
	    fprintf (fp_out, "\n" );
	    fprintf (fp_out, "  Node  Unknown\n" );
	    fprintf (fp_out, "\n" );
	    for ( i = 0; i <= NSUB; i++ )
	    {
	      fprintf (fp_out, "  %8ld  %8d\n", i, indx[i] );
	    }

	    return;
//***************************  PARALLELIZED CODE **************************
	} else {
		 long int i;
      /*
        Set the value of XN, the locations of the nodes.
      */

        #pragma omp parallel for
        for ( i = 0; i <= NSUB; i++ )
        {
          xn[i]  =  ( ( double ) ( NSUB - i ) * xl 
                    + ( double )          i   * xr ) 
                    / ( double ) ( NSUB );
        }
      /*
        Set the lengths of each subinterval.
      */

        #pragma omp parallel for
        for ( i = 0; i < NSUB; i++ )
        {
          h[i] = xn[i+1] - xn[i];
          xquad[i] = 0.5 * ( xn[i] + xn[i+1] );
          node[0+i*2] = i;
          node[1+i*2] = i + 1;
          
        }

        fprintf ( fp_out,"\n" );
        fprintf (fp_out, "  Node      Location\n" );
        fprintf (fp_out, "\n" );
        for (i=0;i<=NSUB;i++)
        {
            fprintf (fp_out, "  %8ld  %14f \n", i, xn[i] );
        }
        fprintf (fp_out, "\n" );
        fprintf (fp_out, "Subint    Length\n" );
        fprintf ( fp_out,"\n" );
        for (i=0;i<=NSUB;i++)
        {
            fprintf (fp_out, "  %8ld  %14f\n", i+1, h[i] );
        }
        fprintf (fp_out, "\n" );
        fprintf (fp_out, "Subint    Quadrature point\n" );
        fprintf ( fp_out,"\n" );
        for (i=0;i<=NSUB;i++)
        {
            fprintf ( fp_out,"  %8ld  %14f\n", i+1, xquad[i] );
        }
        fprintf ( fp_out,"\n" );
        fprintf ( fp_out,"Subint  Left Node  Right Node\n" );
        fprintf (fp_out, "\n" );
        for (i=0;i<=NSUB;i++)
        {
            fprintf (fp_out, "  %8ld  %8d  %8d\n", i+1, node[0+i*2], node[1+i*2] );
        }


      /*
        Starting with node 0, see if an unknown is associated with
        the node.  If so, give it an index.
      */
        nu = 0;
      /*
        Handle first node.
      */
        i = 0;
        if ( ibc == 1 || ibc == 3 )
        {
          indx[i] = -1;
        }
        else
        {
          nu = nu + 1;
          indx[i] = nu;
        }
      /*
        Handle nodes 1 through nsub-1
      */

        /* cannot parallelize due to nu */
        for ( i = 1; i < NSUB; i++ )
        {
          nu = nu + 1;
          indx[i] = nu;
        }
      /*
        Handle the last node.
      /*/
        i = NSUB;

        if ( ibc == 2 || ibc == 3 )
        {
          indx[i] = -1;
        }
        else
        {
          nu = nu + 1;
          indx[i] = nu;
        }

        fprintf ( fp_out,"\n" );
        fprintf ( fp_out,"  Number of unknowns NU = %8d\n", nu );
        fprintf (fp_out, "\n" );
        fprintf (fp_out, "  Node  Unknown\n" );
        fprintf (fp_out, "\n" );
        for ( i = 0; i <= NSUB; i++ )
        {
          fprintf (fp_out, "  %8ld  %8d\n", i, indx[i] );
        }

        return;
	}
  }

/******************************************************************************/
/*
  Purpose: 

    INIT assigns values to variables which define the problem.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Output, int *IBC.
    IBC declares what the boundary conditions are.
    1, at the left endpoint, U has the value UL,
       at the right endpoint, U' has the value UR.
    2, at the left endpoint, U' has the value UL,
       at the right endpoint, U has the value UR.
    3, at the left endpoint, U has the value UL,
       and at the right endpoint, U has the value UR.
    4, at the left endpoint, U' has the value UL,
       at the right endpoint U' has the value UR.

    Output, int *NQUAD.
    The number of quadrature points used in a subinterval.
    This code uses NQUAD = 1.

    Output, double *UL.
    If IBC is 1 or 3, UL is the value that U is required
    to have at X = XL.
    If IBC is 2 or 4, UL is the value that U' is required
    to have at X = XL.

    Output, double *UR.
    If IBC is 2 or 3, UR is the value that U is required
    to have at X = XR.
    If IBC is 1 or 4, UR is the value that U' is required
    to have at X = XR.

    Output, double *XL.
    XL is the left endpoint of the interval over which the
    differential equation is being solved.

    Output, double *XR.
    XR is the right endpoint of the interval over which the
    differential equation is being solved.
*/
  void init (){
  /*
    IBC declares what the boundary conditions are.
  */
    ibc = 1;
  /*
    NQUAD is the number of quadrature points per subinterval.
    The program as currently written cannot handle any value for
    NQUAD except 1.
  */
    nquad = 1;
  /*
    Set the values of U or U' at the endpoints.
  */
    ul = 0.0;
    ur = 1.0;
  /*
    Define the location of the endpoints of the interval.
  */
    xl = 0.0;
    xr = 1.0;
  /*
    Print out the values that have been set.
  */
    fprintf (fp_out, "\n" );
    fprintf ( fp_out,"  The equation is to be solved for\n" );
    fprintf ( fp_out,"  X greater than XL = %f\n", xl );
    fprintf ( fp_out,"  and less than XR = %f\n", xr );
    fprintf ( fp_out,"\n" );
    fprintf (fp_out, "  The boundary conditions are:\n" );
    fprintf (fp_out, "\n" );

    if ( ibc == 1 || ibc == 3 )
    {
      fprintf (fp_out, "  At X = XL, U = %f\n", ul );
    }
    else
    {
      fprintf ( fp_out,"  At X = XL, U' = %f\n", ul );
    }

    if ( ibc == 2 || ibc == 3 )
    {
      fprintf (fp_out, "  At X = XR, U = %f\n", ur );
    }
    else
    {
      fprintf (fp_out, "  At X = XR, U' = %f\n", ur );
    }

    fprintf (fp_out, "\n" );
    fprintf (fp_out, "  Number of quadrature points per element is %d\n", nquad );

    return;
  }
/******************************************************************************/

/*
  Purpose:

    OUTPUT prints out the computed solution.

  Discussion:

    We simply print out the solution vector F, except that, for
    certain boundary conditions, we are going to have to get the
    value of the solution at XL or XR by using the specified
    boundary value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Input, double F(NU).
    ASSEMBLE stores into F the right hand side of the linear
    equations.
    SOLVE replaces those values of F by the solution of the
    linear equations.

    Input, int IBC.
    IBC declares what the boundary conditions are.
    1, at the left endpoint, U has the value UL,
       at the right endpoint, U' has the value UR.
    2, at the left endpoint, U' has the value UL,
       at the right endpoint, U has the value UR.
    3, at the left endpoint, U has the value UL,
       and at the right endpoint, U has the value UR.
    4, at the left endpoint, U' has the value UL,
       at the right endpoint U' has the value UR.

    Input, int INDX[NSUB+1].
    For a node I, INDX(I) is the index of the unknown
    associated with node I.
    If INDX(I) is equal to -1, then no unknown is associated
    with the node, because a boundary condition fixing the
    value of U has been applied at the node instead.
    Unknowns are numbered beginning with 1.
    If IBC is 2 or 4, then there is an unknown value of U
    at node 0, which will be unknown number 1.  Otherwise,
    unknown number 1 will be associated with node 1.
    If IBC is 1 or 4, then there is an unknown value of U
    at node NSUB, which will be unknown NSUB or NSUB+1,
    depending on whether there was an unknown at node 0.

    Input, int NSUB.
    The number of subintervals into which the interval [XL,XR] is broken.

    Input, int NU.
    NU is the number of unknowns in the linear system.
    Depending on the value of IBC, there will be NSUB-1,
    NSUB, or NSUB+1 unknown values, which are the coefficients
    of basis functions.

    Input, double UL.
    If IBC is 1 or 3, UL is the value that U is required
    to have at X = XL.
    If IBC is 2 or 4, UL is the value that U' is required
    to have at X = XL.

    Input, double UR.
    If IBC is 2 or 3, UR is the value that U is required
    to have at X = XR.
    If IBC is 1 or 4, UR is the value that U' is required
    to have at X = XR.

    Input, double XN(0:NSUB).
    XN(I) is the location of the I-th node.  XN(0) is XL,
    and XN(NSUB) is XR.
*/
  void output (){

	 int i;

	double u;

	fprintf (fp_sol,"\n" );
	fprintf (fp_sol,"  Computed solution coefficients:\n" );
	fprintf (fp_sol, "\n" );
	fprintf (fp_sol,"  Node    X(I)        U(X(I))\n" );
	fprintf (fp_sol,"\n" );


	for ( i = 0; i <= NSUB; i++ )
	{
	/*
	If we're at the first node, check the boundary condition.
	*/
	  if ( i == 0 )
	  {
	    if ( ibc == 1 || ibc == 3 )
	    {
	      u = ul;
	    }
	    else
	    {
	      u = f[indx[i]-1];
	    }
	  }
	/*
	If we're at the last node, check the boundary condition.
	*/
	  else if ( i == NSUB )
	  {
	    if ( ibc == 2 || ibc == 3 )
	    {
	      u = ur;
	    }
	    else
	    {
	      u = f[indx[i]-1];
	    }
	  }
	/*
	Any other node, we're sure the value is stored in F.
	*/
	  else
	  {
	    u = f[indx[i]-1];
	  }

	  fprintf ( fp_sol,"  %8d  %8f  %14f\n", i, xn[i], u );
	}

	return;
}
/******************************************************************************/

/*
  Purpose:

    PHI evaluates a linear basis function and its derivative.

  Discussion:

    The evaluation is done at a point X in an interval [XLEFT,XRITE].

    In this interval, there are just two nonzero basis functions.
    The first basis function is a line which is 1 at the left
    endpoint and 0 at the right.  The second basis function is 0 at
    the left endpoint and 1 at the right.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Input, int IL, the index of the basis function.
    1, the function which is 1 at XLEFT and 0 at XRITE.
    2, the function which is 0 at XLEFT and 1 at XRITE.

    Input, double X, the evaluation point.

    Output, double *PHII, *PHIIX, the value of the
    basis function and its derivative at X.

    Input, double XLEFT, XRITE, the left and right
    endpoints of the interval.
*/
  void phi ( int il, double x, double *phii, double *phiix, double xleft, 
    double xrite ){

    if ( xleft <= x && x <= xrite )
    {
      if ( il == 1 )
      {
        *phii = ( xrite - x ) / ( xrite - xleft );
        *phiix =         -1.0 / ( xrite - xleft );
      }
      else
      {
        *phii = ( x - xleft ) / ( xrite - xleft );
        *phiix = 1.0          / ( xrite - xleft );
      }
    }
  /*
    If X is outside of the interval, just set everything to 0.
  */
    else
    {
      *phii  = 0.0;
      *phiix = 0.0;
    }

    return;
  }
/******************************************************************************/

/*
  Purpose:

    PP evaluates the function P in the differential equation.

  Discussion:

    The function P appears in the differential equation as;

      - d/dx (p du/dx) + q u  =  f

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double PP, the value of the function.
*/
  double pp ( double x ){
    double value;

    value = 1.0;

    return value;
  }
/******************************************************************************/

/*
  Purpose:

    PRSYS prints out the tridiagonal linear system.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameter:

    Input, double ADIAG(NU), the "diagonal" coefficients.  That is, 
    ADIAG(I) is the coefficient of the I-th unknown in the I-th equation.

    Input, double ALEFT(NU), the "left hand" coefficients.  That is, ALEFT(I) 
    is the coefficient of the (I-1)-th unknown in the I-th equation.
    There is no value in ALEFT(1), since the first equation
    does not refer to a "0-th" unknown.

    Input, double ARITE(NU).
    ARITE(I) is the "right hand" coefficient of the I-th
    equation in the linear system.  ARITE(I) is the coefficient
    of the (I+1)-th unknown in the I-th equation.  There is
    no value in ARITE(NU) because the NU-th equation does not
    refer to an "NU+1"-th unknown.

    Input, double F(NU).
    ASSEMBLE stores into F the right hand side of the linear
    equations.
    SOLVE replaces those values of F by the solution of the
    linear equations.

    Input, int NU.
    NU is the number of unknowns in the linear system.
    Depending on the value of IBC, there will be NSUB-1,
    NSUB, or NSUB+1 unknown values, which are the coefficients
    of basis functions.
*/
  void prsys (){

    int i;

    fprintf (fp_out, "\n" );
    fprintf (fp_out,"Printout of tridiagonal linear system:\n" );
    fprintf (fp_out,"\n" );
    fprintf (fp_out,"Equation  ALEFT  ADIAG  ARITE  RHS\n" );
    fprintf (fp_out,"\n" );

    /* print statments left unparallelized for speed up */
    for ( i = 0; i < nu; i++ )
    {
      fprintf (fp_out, "  %8d  %14f  %14f  %14f  %14f\n",
        i + 1, aleft[i], adiag[i], arite[i], f[i] );
    }

    return;
  }
/******************************************************************************/

/*
  Purpose: 

    QQ evaluates the function Q in the differential equation.

  Discussion:

    The function Q appears in the differential equation as:

      - d/dx (p du/dx) + q u  =  f

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, double X, the argument of the function.

    Output, double QQ, the value of the function.
*/
  double qq ( double x ){
    double value;

    value = 0.0;

    return value;
  }
/******************************************************************************/

/*
  Purpose: 

    SOLVE solves a tridiagonal matrix system of the form A*x = b.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    29 May 2009

  Author:

    C version by John Burkardt

  Parameters:

    Input/output, double ADIAG(NU), ALEFT(NU), ARITE(NU).
    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
    left and right entries of the equations.
    On output, ADIAG and ARITE have been changed in order
    to compute the solution.
    Note that for the first equation, there is no ALEFT
    coefficient, and for the last, there is no ARITE.
    So there is no need to store a value in ALEFT(1), nor
    in ARITE(NU).

    Input/output, double F(NU).
    On input, F contains the right hand side of the linear
    system to be solved.
    On output, F contains the solution of the linear system.

    Input, int NU, the number of equations to be solved.
*/
  void solve (){

    int i;
  /*
    Carry out Gauss elimination on the matrix, saving information
    needed for the backsolve.
  */
    arite[0] = arite[0] / adiag[0];

    for ( i = 1; i < nu - 1; i++ )
    {
      adiag[i] = adiag[i] - aleft[i] * arite[i-1];
      arite[i] = arite[i] / adiag[i];
    }
    adiag[nu-1] = adiag[nu-1] - aleft[nu-1] * arite[nu-2];
  /*
    Carry out the same elimination steps on F that were done to the
    matrix.
  */
    f[0] = f[0] / adiag[0];

    for ( i = 1; i < nu; i++ )
    {
    	    //printf("%f\n",f[i]);

      f[i] = ( f[i] - aleft[i] * f[i-1] ) / adiag[i];
    }
    
  /*
    And now carry out the steps of "back substitution".
  */
    for ( i = nu - 2; 0 <= i; i-- )
    {
      f[i] = f[i] - arite[i] * f[i+1];
    }

    return;
  }
/******************************************************************************/

/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
  void timestamp ( void ){
    
  # define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
  //  size_t len;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    /*len =*/ strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    fprintf ( fp_out,"%s\n", time_buffer );

    return;
  # undef TIME_SIZE
  }