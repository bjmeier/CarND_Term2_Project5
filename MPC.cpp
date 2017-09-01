#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 11;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
const double mph_to_mps = 0.44704;
const double max_delta = 25 * M_PI / 180; // max allowable steering angle is 25 degrees

size_t x_start     = 0;
size_t y_start     = x_start + N;
size_t psi_start   = y_start + N;
size_t v_start     = psi_start + N;
size_t cte_start   = v_start + N;
size_t epsi_start  = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start     = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {	  
	fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int i = 2; i < N; i++) {
      //  Reference value.  Exponent of 4 was selected to make small errors small, and large errors large
      fg[0] +=  1 * CppAD::pow(vars[cte_start + i], 4);
      
      //  At 50 m/s a one degree error will creat a cross track error of 0.85 m in 0.1s with a value of 10
      fg[0] +=  10 * CppAD::pow(vars[epsi_start + i], 2);
    }
    
    // Minimize the value gap between sequential actuations.
   for (int i = 1; i < N-1; i++) {
      fg[0] += 1000*CppAD::pow(vars[delta_start + i] - vars[delta_start + i - 1], 2); //
      fg[0] += 100*CppAD::pow(vars[cte_start  + i] - vars[cte_start  + i - 1], 2); // 5
    }
    
    // Using lateral acceleration to set desired speed
    AD<double> Rmax = 72;
    for (int i = 2; i < N; i++) {
      AD<double> xfit = vars[x_start + i];
      AD<double> dydxfit  = 3 * coeffs[3] * xfit * xfit + 2 * coeffs[2] * xfit + coeffs[1];
      AD<double> dydx2fit = CppAD::fabs(6 * coeffs[3] * xfit + 2 * coeffs[2]);
      AD<double> nearzero(0.0001);
      AD<double> Rfit = CppAD::fabs(CppAD::pow(1 + dydxfit * dydxfit, 1.5) / CppAD::CondExpGt(dydx2fit, nearzero, dydx2fit, nearzero));  
      Rfit = CppAD::CondExpGt(Rfit, Rmax, Rmax, Rfit);
      AD<double> vdes = CppAD::pow(40*Rfit,0.5);
      fg[0] +=  0.2*CppAD::pow(vars[v_start + i]  - vdes, 2);
    }

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + x_start]    = vars[x_start];
    fg[1 + y_start]    = vars[y_start];
    fg[1 + psi_start]  = vars[psi_start];
    fg[1 + v_start]    = vars[v_start];
    fg[1 + cte_start]  = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];
    
    //std::cout << "set costs and init values \n\n";

     // The rest of the constraints
     for (int i = 0; i < N - 1; i++) {
    	 // t + 1
          AD<double> x1   = vars[x_start    + i + 1];
          AD<double> y1   = vars[y_start    + i + 1];
          AD<double> psi1 = vars[psi_start  + i + 1];
          AD<double> v1   = vars[v_start    + i + 1];
          AD<double> cte1 = vars[cte_start  + i + 1];
          AD<double> epsi1= vars[epsi_start + i + 1];

          // t
          AD<double> x0    = vars[x_start    + i];
          AD<double> y0    = vars[y_start    + i];
          AD<double> psi0  = vars[psi_start  + i];
          AD<double> v0    = vars[v_start    + i];
          AD<double> cte0  = vars[cte_start  + i];
          AD<double> epsi0 = vars[epsi_start + i];

          AD<double> delta0 = vars[delta_start + i];
          AD<double> a0 = vars[a_start + i];
          
          AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * x0 * x0 + coeffs[3] * x0 * x0 * x0;
          AD<double> psides0 = CppAD::atan(3*coeffs[3]*x0*x0 + 2*coeffs[2]*x0 + coeffs[1]);
          
          AD<double> zero(0.);
          AD<double> vavg = v0 + 1 / 2 * a0;
          
          AD<double> d0x = x1 - (x0 + vavg * CppAD::cos(psi0) * dt);
          AD<double> d0y = y1 - (y0 + vavg * CppAD::sin(psi0) * dt);
          AD<double> psi1c = psi0 + vavg * delta0 / Lf * dt;
          
          AD<double> adelta0 = CppAD::fabs(delta0);
          AD<double> deladj;
          
          AD<double> small(0.000001);
          deladj = CppAD::CondExpEq(adelta0, zero, small, zero);
          
          AD<double> dn0x = x1 - (x0 + Lf / (delta0+deladj) * (CppAD::sin(psi1c) - CppAD::sin(psi0)));
          AD<double> dn0y = y1 - (y0 + Lf / (delta0+deladj) * (CppAD::cos(psi0)  - CppAD::cos(psi1c)));
          
          AD<double> dy0  = vavg * CppAD::sin(psi0) * dt;
          AD<double> dyn0 = Lf / (delta0+deladj) * (CppAD::cos(psi0)  - CppAD::cos(psi1c));
          
          
          fg[2 + x_start    + i] = CppAD::CondExpGt(adelta0, zero, dn0x, d0x);
          fg[2 + y_start    + i] = CppAD::CondExpGt(adelta0, zero, dn0y, d0y);
          fg[2 + psi_start  + i] = psi1  - psi1c;
          fg[2 + v_start    + i] = v1    - (v0 + a0 * dt);
          fg[2 + cte_start  + i] = cte1  - ((f0 - y0) + CppAD::CondExpGt(adelta0, zero, dyn0, dy0));
          fg[2 + epsi_start + i] = epsi1 - ((psi0 - psides0) + vavg * delta0 / Lf * dt);

     }  
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3] * mph_to_mps; // converts mph to m/s
  double cte  = state[4];
  double epsi = state[5];
  double sa   = - state[6]; // steering angle in radians
  double th   = state[7]; // throttle
  size_t n_vars = N * 6 + ( N - 1) * 2;
  // TODO: Set the number of constraints
  size_t n_constraints = N * 6;
  
  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  
  // Set the initial variable values
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;
  vars[delta_start] = sa;
  vars[a_start]    = th;

  //std::cout << "init values set \n\n";
  
  // Lower and upper limits for x
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  //std::cout << "Dvector set \n\n";
  
  // Set all non-actuators upper and lowerlimits
  // to the max negative and positive values.
  for (int i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  //std::cout << "non act bounds set \n\n";
  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  for (int i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -max_delta;
    vars_upperbound[i] =  max_delta;
  }

  // Acceleration/decceleration upper and lower limits.
  for (int i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] =  1.0;
  }
  
  vars_lowerbound[delta_start]= sa;
  vars_lowerbound[a_start]= th;
  vars_upperbound[delta_start]= sa;
  vars_upperbound[a_start]= th;
  
  // Lower and upper limits for constraints
  // All of these should be 0 except the initial
  // state indices.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  
  constraints_lowerbound[x_start]    = x;
  constraints_lowerbound[y_start]    = y;
  constraints_lowerbound[psi_start]  = psi;
  constraints_lowerbound[v_start]    = v;
  constraints_lowerbound[cte_start]  = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start]    = x;
  constraints_upperbound[y_start]    = y;
  constraints_upperbound[psi_start]  = psi;
  constraints_upperbound[v_start]    = v;
  constraints_upperbound[cte_start]  = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  //std::cout << "ipopt solved \n\n";
  // Check some of the solution values
  ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;

  //By setting dt to the delay and fixing the actuator values at t=0, returning the actuation at dt will help deal with lag.
  
  vector<double> result;
  
  result.clear();
  
  // To deal with delay, sent back NEXT time step to main.
  // Next time step and delay are both 100 ms.
  
  result.push_back(solution.x[delta_start + 1]);
  result.push_back(solution.x[a_start + 1]);
  
  for (int i = 0; i < N - 1; i++)
  {
	  result.push_back(solution.x[x_start + i + 1]);
	  result.push_back(solution.x[y_start + i + 1]);
  }
  
  return result;
}

