#include <iostream>
#include <math.h>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}

VectorXd Tools::ToPolar(const VectorXd &x)
{
  VectorXd ret = VectorXd(3);
  double rho = sqrt(x[0] * x[0] + x[1] * x[1]);
  double phi = atan(x[1] / x[0]);
  double rho_dot = (x[0] * x[2] + x[1] * x[3]) / rho;

  ret << rho, phi, rho_dot;

  return ret;
}

VectorXd Tools::ToCartesian(const VectorXd &x)
{
  VectorXd ret = VectorXd(4);
  ret << x[0] * cos(x[1]),
    x[0] * sin(x[1]),
    x[2] * cos(x[1]),
    x[2] * sin(x[1]);

  return ret;
}

VectorXd Tools::CalculateRMSE(
  const vector<VectorXd> &estimations,
  const vector<VectorXd> &ground_truth)
{
  /**
     TODO:
     * Calculate the RMSE here.
     */

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() != ground_truth.size()
     || estimations.size() == 0)
  {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i)
  {

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
  /**
     TODO:
     * Calculate a Jacobian here.
     */

  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if(fabs(c1) < 0.0001)
  {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    throw std::invalid_argument( "received negative value" );
    return Hj;
  }

  //compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
    -(py/c1), (px/c1), 0, 0,
    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
