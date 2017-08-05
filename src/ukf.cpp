#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // State covariance matrix
  P_ << 1, 0, 0, 0, 0,
       	0, 1, 0, 0, 0,
       	0, 0, 1, 0, 0,
       	0, 0, 0, 1, 0,
       	0, 0, 0, 0, 1;  

  //set state dimension
  n_x = 5;

  // set aug state dimension
  n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z = 3;

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(5,15);

  //create vector for weights
  weights_ = VectorXd(2*n_aug+1);
  
  // define spreading parameter
  lambda = 3 - n_aug;  

  // set weights
  double weight_0 = lambda/(lambda+n_aug);
  weights_(0) = weight_0;
  double weight = 0.5/(lambda+n_aug);
  
  for (int i = 1; i < 2*n_aug+1; i++) {  //2n+1 weights_    
    weights_(i) = weight;
  }
  
  // Used in UpdateLidarEKF function
  // Init matrices
  R_ = MatrixXd(2, 2);
  H_ = MatrixXd(2, 5);
  //measurement covariance matrix - laser
  R_ << 0.0225, 0,
  		0, 0.0225;

  // measurement matrix - laser
  H_ << 1, 0, 0, 0, 0,
  		0, 1, 0, 0, 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ukf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "UKF: " << endl;   
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
       */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[1];

      float x = rho*cos(phi);
      float y = rho*sin(phi);
      
      // set the state vector
      x_ << x,y,0,0,0;      

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;        
    }

  	time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
	
    return;
	
	}

	//compute the time elapsed between the current and previous measurements

	float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0; //delta_t - expressed in seconds
	time_us_ = meas_package.timestamp_;

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
	Prediction(delta_t);
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_ == true) {		
		//cout << "RADAR" << endl;
		UpdateRadar(meas_package);
	}
	else if(meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_ == true){
		//cout << "LIDAR" << endl;
		UpdateLidar(meas_package);			// UKF based calculations
		//UpdateLidarEKF(meas_package);		// EKF based calculations : Since lidar calc are
											//	linear, EKF equations can be used here too.
	}  
	// print the output
	//cout << "x_ = " << x_ << endl;
	//cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /***************************** GENERATING AUGMENTED SIGMA POINTS ************************/
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  x_aug.fill(0.0);

  //create augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);

  //create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug.fill(0.0);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix  
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
  }

  /************************ PREDICTING SIGMA POINTS **************************/  
  //predict sigma points
  for (int i = 0; i < 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
        
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
	
    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;	
    
    //write predicted sigma point into right column           
    Xsig_pred_(0,i) = px_p;    
    Xsig_pred_(1,i) = py_p;    
    Xsig_pred_(2,i) = v_p;    
    Xsig_pred_(3,i) = yaw_p;    
    Xsig_pred_(4,i) = yawd_p;        
  }

  /****************************** PREDICTED MEAN X AND P **************************/

  //create mean predict state vector
  VectorXd x_pred = VectorXd(5);
  
  //predicted state mean
  x_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }

  //create predicted state covariance matrix
  MatrixXd P_pred = MatrixXd(5, 5);
  //predicted state covariance matrix
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose() ;
	}

  x_ = x_pred;
  P_ = P_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  /******************************** PREDICTING LIDAR MEASUREMENTS **********************************/   
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z-1, 2 * n_aug + 1); // nz-1 = 2 since the lidar measures only px, py
  Zsig.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    
    // measurement model
    Zsig(0,i) = p_x;			// px
    Zsig(1,i) = p_y;    		// py
  }
    
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z-1);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z-1,n_z-1);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z-1,n_z-1); 
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;          
  S = S + R;

  /***************************** UPDATE STATE  ************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z-1);
  Tc.fill(0.0);

  //create vector for incoming lidar measurement
  VectorXd z_ = VectorXd(n_z-1);

  //calculate cross correlation matrix  
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();  

  // sensor measurement
  z_ = meas_package.raw_measurements_;

  //residual
  VectorXd z_diff = z_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  /******************************* NIS Calculations *************************/
}

void UKF::UpdateLidarEKF(MeasurementPackage meas_package) {

	// This function is the same as the EKF equations of updatelidar function

	//create vector for incoming lidar measurement
  	VectorXd z_ = VectorXd(n_z-1);

	// sensor measurement
  	z_ = meas_package.raw_measurements_;
  
	VectorXd z_pred = H_ * x_;
	VectorXd y = z_ - z_pred;
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

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /******************************** PREDICTING RADAR MEASUREMENTS **********************************/  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  /***************************** UPDATE STATE  ************************************/
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  Tc.fill(0.0);

  //create vector for incoming radar measurement
  VectorXd z_ = VectorXd(n_z);

  //calculate cross correlation matrix  
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // sensor measurement
  z_ = meas_package.raw_measurements_;
  //residual
  VectorXd z_diff = z_ - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  /******************************* NIS Calculations *************************/
}
