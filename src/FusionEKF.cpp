#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  // Initialise H_laser
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			  0, 1, 0, 0,
			  0, 0, 1000, 0,
			  0, 0, 0, 1000;

  //the initial transition matrix F_ with dt = 0
	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 0, 0,
			       0, 1, 0, 0,
			       0, 0, 1, 0,
			       0, 0, 0, 1;

  // Initialise H
  ekf_.H_ = MatrixXd(4, 4);
  ekf_.H_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    std::cout << "EKF: " << std::endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    ekf_.Q_ = Eigen::MatrixXd(4,4);

    float px;
    float py;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      std::cout << "Initialize Radar" << std::endl;
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      px = rho*cos(phi);
      py = rho*sin(phi);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      std::cout << "Initialize Laser" << std::endl;
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];

    }
    // To handel the case when the value of px or py takes on a very small value
    if(fabs(px) < 0.0001){
        std::cout << "The initial value of px is very small." << std::endl;
        px = 0.1;
    }

    if(fabs(py) < 0.0001){
        std::cout << "The initial value of py is very small." << std::endl;
        py = 0.1;
    }

    ekf_.x_ << px, py, 0, 0;
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = measurement_pack.timestamp_;

 	float dt_2 = dt * dt;
 	float dt_3 = dt_2 * dt;
 	float dt_4 = dt_3 * dt;

 	//Modify the F matrix so that the time is integrated
 	ekf_.F_(0, 2) = dt;
 	ekf_.F_(1, 3) = dt;

  // set noise values
  float noise_ax = 9;
  float noise_ay = 9;

 	//set the process covariance matrix Q
 	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
 			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
 			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
 			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
     ekf_.hx_ = VectorXd(3);

     float px = ekf_.x_[0];
     float py = ekf_.x_[1];
     float vx = ekf_.x_[2];
     float vy = ekf_.x_[3];

     float rho;
     float phi;
     float rhodot;
     if(fabs(px) < 0.0001 or fabs(py) < 0.0001){
        if(fabs(px) < 0.0001){
          px = 0.0001;
        }
        if(fabs(py) < 0.0001){
          py = 0.0001;
        }
        rho = sqrt(px*px + py*py);
        phi = 0;
        rhodot = 0;

      } else {
        rho = sqrt(px*px + py*py);
        phi = atan2(py,px);
        rhodot = (px*vx + py*vy) /rho;
      }

      ekf_.hx_ << rho, phi, rhodot;

      // Set the value of H_ to Hj with a radar measurement
      Hj_ = tools.CalculateJacobian(ekf_.x_);

      // The update step is not performed when the value of the Jacobian matrix is 0
      if (Hj_.isZero(0)){
        std::cout << "The value of Hj is 0" << std::endl;
        return;
      }

      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
