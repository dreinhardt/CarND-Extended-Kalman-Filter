#include "tools.h"
#include <iostream>

using namespace std;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {		// very similar to the workshop example
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd e(4);
    e << 0,0,0,0;
  	unsigned int esize = estimations.size();

    if(ground_truth.size() == 0){
      cout << "ERR: Ground-truth is empty" << endl;
      return e;
    }
  
    if(esize == 0){
      cout << "ERR: Estimation vector is empty" << endl;
      return e;
    }

    if(esize != ground_truth.size()){
      cout << "ERR: Ground-truth and estimations vectors not equal." << endl;
      return e;
    }

    for(unsigned int i=0; i < esize; ++i){
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array()*diff.array();
      e += diff;
    }

    e = e / esize;
    e = e.array().sqrt();
    return e;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {    // Identical to our tutorials
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
      -(py/c1), (px/c1), 0, 0,
      py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
