#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;  // initial state (location and velocity) - Ort und Geschwindigkeit
  P_ = P_in;  // initial uncertainty 		- Anfangsunschaerfe
  F_ = F_in;  // next state function		- Zustandsübergangsmatrix
  H_ = H_in;  // measurement function		- Messmatrix
  R_ = R_in;  // measurement uncertainty	- Noise, Unschaerfe, measurement covariance matrix
  Q_ = Q_in;  // process covariance matrix	- Kovarianzmatrix, Prozessrauschen, 
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_ ;		// Identical to our tutorials
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  double rho     = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta   = atan2(x_(1), x_(0));
  double rho_dot;
  rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;  // formular from the script

  VectorXd h = VectorXd(3);
  h << rho, theta, rho_dot;		// form new h small
  VectorXd y = z - h;			// go into last function of Radar operation to result the common y
  
  // Normalize the angles - out of tips section
  while (y(1) > M_PI)
    y(1) -= 2 * M_PI;	// subtract 2PI
  
  while (y(1) <= -M_PI)
    y(1) += 2 * M_PI;	// add 2PI
  
  // Proceed with classic update
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
