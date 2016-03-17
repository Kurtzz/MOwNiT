#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define N 20
long double lagrange_coefficients [N];
double newton_coefficients [N];
FILE * result_file;
FILE * coeff_file;

double double_rand_range(double min_n, double max_n)
{
    return (double)rand()/RAND_MAX * (max_n - min_n) + min_n;
}

void rand_nodes(double x [] ,double y [] )
{
    int i;
    srand(time(NULL));
    x[0]=0;
    y[0]=double_rand_range(-100,100);
    for ( i=0;i<N-1;i++)
    {
        x[i+1]=x[i]+4;
        y[i+1] = y[i] + double_rand_range(-10,10);
    }

}

void gsl_polym_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_polynomial, N);

        gsl_interp_init(workspace, x, y, N);
        h = (x[1]-x[0])/step;
        for(i = 0; i <= (step*(N-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_interp_eval(workspace, x, y, xtmp, acc);
            printf("GSL,%lf,%lf \n", xtmp,ytmp);
            fprintf(result_file,"GSL,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }
}

void gsl_spline_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline * workspace = gsl_spline_alloc(gsl_interp_cspline, N);

        gsl_spline_init(workspace, x, y, N);
        h = (x[1]-x[0])/step;
        for(i = 0; i <= (step*(N-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_spline_eval(workspace, xtmp, acc);
            printf("Spline,%lf,%lf \n", xtmp,ytmp);
            fprintf(result_file,"Spline,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_spline_free (workspace);
        gsl_interp_accel_free (acc);
    }
}

void gsl_akima_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_akima, N);

        gsl_interp_init(workspace, x, y, N);
        h = (x[1]-x[0])/step;
        for(i = 0; i <= (step*(N-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_interp_eval(workspace, x, y, xtmp, acc);
            printf("Akima,%lf,%lf \n", xtmp,ytmp);
            fprintf(result_file,"Akima,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }

}
// LaGrange init
void getLagrangeCoeffs(double x[],double y[])  {
    int i,j,k;
    // c[] are the coefficients of P(x) = (x-x[0])(x-x[1])...(x-x[n-1])
    long double c[N+1];
    for (i=0;i<N;i++)
    {
        lagrange_coefficients[i] = 0.0;
    }
    c[0] = 1.0;
    for ( i = 0; i < N; i++)
    {
        for ( j = i; j > 0; j--)
        {
            c[j] = c[j-1] - c[j] * x[i];
        }
        c[0] *= -x[i];
        c[i+1] = 1;
    }
    long double tc[N+1];
    for ( i = 0; i < N; i++) {
        // d = (x[i]-x[0])...(x[i]-x[i-1])(x[i]-x[i+1])...(x[i]-x[n-1])
        long double d = 1.0;
        for ( j = 0; j < N; j++) {
            if (i != j) {
                d *= (x[i] - x[j]);
            }
        }
        if (d == 0.0) {
            // This happens only when two abscissas are identical.
            for ( k = 0; k < N; ++k) {
                if ((i != k) && (x[i] == x[k])) {
                    printf("LocalizedFormats.IDENTICAL_ABSCISSAS_DIVISION_BY_ZERO,%d %d %lf",i, k, x[i]);
                    exit(-1);
                }

            }
        }

        long double t = y[i] / d;
        // Lagrange polynomial is the sum of n terms, each of which is a
        // polynomial of degree n-1. tc[] are the coefficients of the i-th
        // numerator Pi(x) = (x-x[0])...(x-x[i-1])(x-x[i+1])...(x-x[n-1]).
        tc[N-1] = c[N];     // actually c[n] = 1
        lagrange_coefficients[N-1] += t * tc[N-1];
        for ( j = N-2; j >= 0; j--) {
            tc[j] = c[j+1] + tc[j+1] * x[i];
            lagrange_coefficients[j] += (t * tc[j]);
        }
    }
    printf("lagrange_coefficients \n");
    fprintf(coeff_file,"lagrange_coefficients \n");
    for (i=0;i<N;i++)
    {
        printf("  %Lf x^%d",lagrange_coefficients[i],i );
        fprintf(coeff_file,"%Lf*x^%d +",lagrange_coefficients[i],i );
    }
    printf("\n");
    fprintf(coeff_file,"\n");
    return;
}
void print_newton(double x[])
{
    int i,j;
    printf("Newton_coefficients \n");
    fprintf(coeff_file,"Newton_coefficients \n");
    printf("%lf *1 +",newton_coefficients[0]);
    fprintf(coeff_file,"%lf *1 +",newton_coefficients[0]);
    for(i=1;i<N;i++)
    {
        printf("%lf",newton_coefficients[i] );
        fprintf(coeff_file,"%lf",newton_coefficients[i] );
        for (j=1;j<=i;j++)
        {
            printf("(x-%lf)",x[j-1] );
            fprintf(coeff_file,"(x-%lf)",x[j-1] );
        }
        if (i!=N-1) printf(" + ");
        if (i!=N-1) fprintf(coeff_file," + ");
    }
    printf("\n");
}

void getNetwonCoeffs(double x[],double y[])  {
    //  ilorazy roznicowe nie cala jest potrzebna, tylko  polowa
    double newton_coeffs [N-1][N-1];
    int i,j;

    // zero table
    for (i=0;i<(N-1);i++)
        for (j=0;j<(N-1);j++)
            newton_coeffs[i][j]=0;

    // first columns
    for (i=0;i<(N-1);i++)
        newton_coeffs[0][i]=(y[i+1]-y[i])/(x[i+1]-x[i]);


    for (i=1;i<(N-1);i++)
    {
        int counter=0;
        for(j=(N-2-i);j>=0;j--)
        {
            newton_coeffs[i][j]=(newton_coeffs[i-1][j+1]-newton_coeffs[i-1][j])/(x[(N-1-counter)]-x[(N-2-i-counter)]  );
            counter++;
        }
    }

    newton_coefficients[0]=y[0];
    for (i=1;i<N;i++)
    {
        newton_coefficients[i]=newton_coeffs[i-1][0];
    }
    print_newton(x);
}

double horner(long double *coeffs, int s,long double x)
{
    int i;
    long double res = 0.0;

    for(i=s-1; i >= 0; i--)
    {
        res = res * x + coeffs[i];
    }
    return res;
}

void eval_lagrange_interpol(double x[],double y[],int step)
{
    long double ytmp;
    long double xtmp,h;
    int i;
    h = (x[1]-x[0])/step;
    for(i = 0; i <= (step*(N-1)); i++)
    {
        xtmp = x[0] + h *i;
        ytmp = horner(lagrange_coefficients,N,xtmp);
        printf("Lagrange,%Lf,%Lf \n",xtmp,ytmp);
        fprintf(result_file,"Lagrange,%Lf,%Lf \n",xtmp,ytmp);
    }
}

double count_newton(double *coeffs,double *x, int s, double xtmp)
{
    double p = 1;
    int i;
    double result=0;
    for (i=0;i<s;i++)
    {
        result=result + (p*coeffs[i]);
        p = p * (xtmp-x[i]);
    }
    return result;

}
void eval_newton_interpol(double x[],double y[],int step)
{
    double ytmp;
    double xtmp,h;
    int i;
    h = (x[1]-x[0])/step;
    for(i = 0; i <= (step*(N-1)); i++)
    {
        xtmp = x[0] + h *i;
        ytmp = count_newton(newton_coefficients,x,N,xtmp);
        printf("Newton,%lf,%lf \n", xtmp,ytmp);
        fprintf(result_file,"Newton,%lf, %lf \n", xtmp,ytmp);
    }
}

int main ()
{
    int i;
    if  ((result_file = fopen ("result.txt", "w+")) == NULL )
    {
        printf("Result file open error \n");
        exit(-1);
    }
    if( (coeff_file = fopen ("coeff.txt", "w+")) == NULL)
    {
        printf("Coeff file open error \n");
        exit(-1);
    }
    double x[N];
    double y[N];
    rand_nodes(x,y);
    getLagrangeCoeffs(x,y);
    getNetwonCoeffs(x,y);
    printf("alg,x,y \n");
    fprintf(result_file,"alg,x,y \n");
    for ( i=0;i<N;i++)
    {
        printf("Gen_points,%lf,%lf \n", x[i],y[i]);
        fprintf(result_file,"Gen_points,%lf,%lf \n", x[i],y[i]);
    }

    eval_newton_interpol(x,y,4);
    eval_lagrange_interpol(x,y,4);
    gsl_polym_interpol(x,y,4);
    gsl_akima_interpol(x,y,4);
    gsl_spline_interpol(x,y,4);

    return 0;
}
