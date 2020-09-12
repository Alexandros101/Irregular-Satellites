/*
 * lessnasty.c
 */
#include <stdio.h>
#include <math.h>

void Euler(double omega, double N)
{
  double t;
  double dt;
  double x;
  double dx;
  int j;
  //declare all variables

  t = 0.0;
  x = 0.0;
  dx = 0.0;
  printf("%15.8f %15.8f %15.8f\n", t, x, dx);

  //initialize and print initial state
  dt = 6*M_PI / N;
  //the final time is 6pi so the time step size will
  //be 6pi over the number of steps
  for (j = 0; j < N; ++j)
  {
    double xold = x;
    double dxold = dx;
    
    //The Euler integration method: position moves
    //forward, then velocity, then time and the final
    //states are printed
    x += dt * dxold;
    dx += dt * (cos(omega*t) - xold - omega*dxold);
    t += dt;
    printf("%15.8f %15.8f %15.8f\n", t, x, dx);
  }
  printf("\n");
}

int main(double argc, int arg)
{
  double omega;
  int N;
  //Takes in two arguments, no error calls
  //as requested
  printf("Please input the value of frequency: ", omega);
  scanf("%lf", &omega);
  printf("Please input the number of steps: ", N);
  scanf("%d", &N);
  //Euler is executed and prints are in it
  Euler(omega,N);
  return 0;
}