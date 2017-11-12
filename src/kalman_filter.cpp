#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;

  cout << "Kalman Filter Predict" << endl;
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  cout << "Kalman Filter Update Begins" << endl;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd S_inv = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * S_inv;

  x_ = x_ +(K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);

  P_ = (I - K * H_) * P_;

  cout << "Kalman Filter Update" << endl;
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

  VectorXd z_pred = hx_;
  VectorXd y = z - z_pred;

  // if value of phi is not with in the range of   (-pi, pi), then put it in that range
  bool in_pi = false;
  while (in_pi == false) {
    if (y(1) > 3.14159) {
      y(1) = y(1) - 6.2831;
    }
    else if (y(1) < -3.14159) {
      y(1) = y(1) + 6.2831;
    } else {
      in_pi = true;
    }
  }

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd S_inv = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * S_inv;

  x_ = x_ + (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  P_ = (I - K * H_) * P_;

  cout << "Extended Kalman Filter Update" << endl;
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}
