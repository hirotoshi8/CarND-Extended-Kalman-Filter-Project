
#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define IS_RADER_ARIVE true

//Helper tools
Tools tool;


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
  R_laser_ << 0.0225,      0,
                   0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09,      0,    0,
                 0, 0.0009,    0,
                 0,      0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  //state covariance matrix P
  ekf_.P_ = MatrixXd::Identity(4, 4);

  //the initial transition matrix F_
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
    
  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
#if 0
  ekf_.Q_ << 0, 0,  0,    0,
             0, 0,  0,    0,
	         0, 0,  0,    0,
	         0, 0,  0,    0;
#endif
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
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
    cout << "********** EKF: Initialization ********** " << endl;
    ekf_.x_ = VectorXd(4);
    //ekf_.x_ << 1, 1, 1, 1;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		/**
		Convert radar from polar to cartesian coordinates and initialize state.
		*/
		cout << "EKF with first Radar: " << endl;
		//get the variables
		double ro = measurement_pack.raw_measurements_[0];
		double theta = measurement_pack.raw_measurements_[1];
		double ro_dot = measurement_pack.raw_measurements_[2];

		double px = ro * cos(theta);
		double py = ro * sin(theta);
		double vx = ro_dot * cos(theta);
		double vy = ro_dot * sin(theta);

		//set the state with the initial location (and zero velocity ?)
  		ekf_.x_ << px, py, vx, vy;
	  	
		//Update the Covariance Matrix "R" of Radar
		ekf_.R_ = MatrixXd(3,3);
		ekf_.R_ = R_radar_;
	  //ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
		/**
		Initialize state.
		*/
		cout << "EKF with first Laser: " << endl;
		//cout << measurement_pack.raw_measurements_ << endl;
		//set the state with the initial location (and zero velocity ?)
		double px = measurement_pack.raw_measurements_[0];
		double py = measurement_pack.raw_measurements_[1];
		double vx = 0;
		double vy = 0;
      
		//cout << measurement_pack.raw_measurements_ << endl;
		ekf_.x_ << px, py, vx, vy;
            
		//Update the Covariance Matrix "R" of Radar
		ekf_.R_ = MatrixXd(2,2);
		ekf_.R_ = R_laser_;
	  //ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
	}
    
	previous_timestamp_ = measurement_pack.timestamp_;
    
	// done initializing, no need to predict or update
	is_initialized_ = true;

	return;
  }

  /*****************************************************************************
  *  After Initialization
  ****************************************************************************/

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
  //debug
  cout << "********** Prediction: ********** " << endl;

  double noise_ax = 9;
  double noise_ay = 9;

  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;//dt - expressed in seconds
  //Update the previous time
  previous_timestamp_ = measurement_pack.timestamp_;
  
  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  //set the process covariance matrix Q (Process noise)
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax,               0, dt_3/2*noise_ax,               0,
                            0, dt_4/4*noise_ay,               0, dt_3/2*noise_ay,
	          dt_3/2*noise_ax,               0,   dt_2*noise_ax,               0,
	                        0, dt_3/2*noise_ay,               0,   dt_2*noise_ay;

  // Predict
  cout << " Predicte Execute: " << endl;
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
   //debug
  cout << "********** Update: ********** " << endl;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	  //Radar updates
#if IS_RADER_ARIVE
	  //Update the Transfer Matrix "H" with Jacobian
	  ekf_.H_ = MatrixXd(3,4);
	  ekf_.H_ = tool.CalculateJacobian(ekf_.x_);
    
	  //Update the Covariance Matrix "R" of Radar
	  ekf_.R_ = MatrixXd(3,3);
	  ekf_.R_ = R_radar_;
    
	  //Measurement of the Radar
	  VectorXd z_rader(3);
	  z_rader(0) = measurement_pack.raw_measurements_[0]; //ro
	  z_rader(1) = measurement_pack.raw_measurements_[1]; //theta
	  z_rader(2) = measurement_pack.raw_measurements_[2]; //ro_dot
    
	  cout << "Radar Update: EKF" << endl;
	  ekf_.UpdateEKF(z_rader);
#endif

  } else {
	  //Update the Transfer Matrix "H"
	  ekf_.H_ = MatrixXd(2,4);
	  ekf_.H_ = H_laser_;
	  //Update the Covariance Matrix "R" of Lider
	  ekf_.R_ = MatrixXd(2,2);
	  ekf_.R_ = R_laser_;

	  cout << "Lider Update: EKF" << endl;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
