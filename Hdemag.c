#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Temporary for code evaluation
#include <time.h>

void clear(double** arr, int n, int m){
	for(int i = 0; i<n; i++){
		free(arr[i]);
	}
	free(arr);
}

void update(double** u_old, double** u, int n, int m){
	for(int i = 0; i<n; i++){
		for(int j = 0; j<m; j++){
			u_old[i][j] = u[i][j];
		}
	}
}

double** alloc(int n, int m, double value){
	double** arr = (double**)malloc(n*sizeof(double*));
	for(int i = 0; i<n; i++){
		arr[i] = (double*)malloc(m*sizeof(double));
		for(int j = 0; j<m; j++){
			arr[i][j] = value;
		}
	}

	return arr;
}

// Lx, Ly, dx, max_iter 
double** run(double Lx, double Ly, double dx, int max_iter){

	// Convergence tolerance
	double tol = 1e-6;

	// Calculate the number of grid points
	int nx = Lx/dx + 1;
	int ny = Ly/dx + 1;

	// Initialize the potential
	double** u = alloc(ny,nx,0);
	double** u_old = alloc(ny,nx,0);

	// Define the magnetization
	double** Mx = alloc(ny,nx,1);
	double** My = alloc(ny,nx,0);

	// Preform iterations to solve the interior points
	for(int iter = 0; iter<max_iter; iter++){
		update(u_old,u,nx,ny);

		// Update the interior points
		for(int i = 1; i<nx-1; i++){
			for(int j = 1; j<ny-1; j++){

				double F = (Mx[j][i+1] - Mx[j][i] + \
							My[j+1][i] - My[j][i])/dx;

				u[j][i] = 0.25 * (u_old[j+1][i] + u_old[j-1][i] + \
								  u_old[j][i+1] + u_old[j][i-1] - \
								  dx*dx*F);

			} // End of y loop
		} // End of x loop

		// Apply Neumann boundary conditions
		for(int j = 0; j<ny; j++){
			u[j][0] = u[j][1] - dx*Mx[j][0];			// Left boundary (x = 0)
			u[j][nx-1] = u[j][nx-2] + dx*Mx[j][nx-1];	// Right boundary (x = Lx)
		}

		for(int i= 0; i<nx; i++){
			u[0][i] = u[1][i] - dx*My[0][i];			// Bottom boundary (y = 0)
			u[ny-1][i] = u[ny-2][i] + dx*My[ny-1][i];	// Top boundary (y = Ly)
		}

		// Check for convergence
		double max = 0;
		for(int i = 0; i<nx; i++){
			for(int j = 0; j<ny; j++){
				if(fabs(u[j][i]-u_old[j][i])>max){
					max = fabs(u[j][i]-u_old[j][i]);
				}
			}
		}

		if(max<tol)
		{
			printf("Converged after %d itearations.\n", iter);
			break;
		}

	} // End of attemp iteration loop
	clear(u_old,nx,ny);
	clear(Mx,nx,ny);
	clear(My,nx,ny);

	return u;
}

int main(void)
{

	// Ready clock for computing time
	clock_t start_time, end_time;
	double cpu_time_used;

	// Dimensions and spacing
	double Lx = 1.0;
	double Ly = 1.0;

	double dx = 0.01;

	// Calculate the number of grid points
	int nx = Lx/dx + 1;
	int ny = Ly/dx + 1;

	// Iterative solver parameters
	double max_iter = 1e5;

	
	start_time = clock(); // tic
	
	// Run
	double** V = run(Lx,Ly,dx,max_iter);
	
	end_time = clock(); // toc
	cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC * 1000.0;

	// Print the CPU time used
	printf("CPU time used: %f milliseconds\n", cpu_time_used);

	clear(V,nx,ny);

	printf("%s\n", "Memory cleared!");
	return 0;
}