#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd S = H_ * PHt + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = PHt * Si;

  //Update the posterior state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
	TODO:
	* update the state by using Extended Kalman Filter equations
	*/
	VectorXd z_pred(3);
  
	double px = x_[0];
	double py = x_[1];
	double vx = x_[2];
	double vy = x_[3];
  
	double p2     = px*px + py*py;
	double sqrt_p = sqrt(p2); // pow(p2, 0.5);
  
	double ro_hat = sqrt_p;
	double theta_hat = atan2(py,px);
	double ro_dot_hat = (px*vx + py*vy) / sqrt_p;

	//check division by zero
#if 0
	if (p2 == 0) {
		ro_dot_hat = 0;
	}else {
		ro_dot_hat = (px*vx + py*vy) / sqrt_p;
	}
#endif
	//Calculate the measurement based on the state prediction
	z_pred << ro_hat, theta_hat, ro_dot_hat;

	//Update
	VectorXd y = z - z_pred;	//VectorXd z_pred = H_ * x_;

	while(y[1] > M_PI) y[1] -= 2 * M_PI;
	while(y[1] <-M_PI) y[1] += 2 + M_PI;

	MatrixXd Hjt = H_.transpose();
	MatrixXd S = H_ * P_ * Hjt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Hjt;
	MatrixXd K = PHt * Si;

	//Update the posterior state 
	x_ = x_ + (K * y);
	double x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
