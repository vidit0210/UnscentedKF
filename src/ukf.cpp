#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initialize Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored  (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = .75;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/8.0;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  is_initialized_ = false;

  // state dimension
  n_x_ = 5;

  // augmented state dimension
  n_aug_ = 7;

  // spreading parameter
  lambda_ = 3 - n_aug_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_ + 1);

  // weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }

  // R matrix for Radar calculations
  R_rdr_ = MatrixXd(3,3);
  R_rdr_.fill(0);
  R_rdr_(0,0) = std_radr_*std_radr_;
  R_rdr_(1,1) = std_radphi_*std_radphi_;
  R_rdr_(2,2) = std_radrd_*std_radrd_;

  // R matrix for Laser calculations
  R_lsr_ = MatrixXd(2,2);
  R_lsr_.fill(0);
  R_lsr_(0,0) = std_laspx_*std_laspx_;
  R_lsr_(1,1) = std_laspy_*std_laspy_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if(!is_initialized_) {
    // initialize ukf with first measurement
    cout << "Initializing UKF" << endl;
    x_ = VectorXd(n_x_);
    x_ << 1,1,0,0,0;
    P_ = MatrixXd::Identity(n_x_,n_x_);

    if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // convert from polar to cartesian coordinates and initialize state
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      x_(0) = rho*cos(theta);
      x_(1) = rho*sin(theta);
      R_ = R_rdr_;
    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      R_ = R_lsr_;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    R_ = R_rdr_;
    UpdateRadar(meas_package);
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
    R_ = R_lsr_;
    UpdateLidar(meas_package);
  }
    
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */

void UKF::Prediction(double delta_t) {
  cout << "Predicting" << endl;
  // Complete this function! Estimate the object's location. Modify the state
  // vector, x_. Predict sigma points, the state, and the state covariance matrix.

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;
  x_aug(n_x_+1) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_x_,n_x_) = std_a_*std_a_;
  P_aug(n_x_+1,n_x_+1) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++) {
    //extract values for better readability
    double x,y,v,yaw,yaw_rate,noise_a,noise_yr;
    x = Xsig_aug(0,i);
    y = Xsig_aug(1,i);
    v = Xsig_aug(2,i);
    yaw = Xsig_aug(3,i);
    yaw_rate = Xsig_aug(4,i);
    noise_a = Xsig_aug(5,i);
    noise_yr = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yaw_rate) > 0.001) {
      px_p = x + v/yaw_rate * ( sin(yaw + yaw_rate*delta_t) - sin(yaw));
      py_p = y + v/yaw_rate * ( cos(yaw) - cos(yaw+yaw_rate*delta_t) );
    }
    else {
      px_p = x + v*delta_t*cos(yaw);
      py_p = y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yaw_rate*delta_t;
    double yawd_p = yaw_rate;

    //add noise
    px_p = px_p + 0.5*noise_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*noise_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + noise_a*delta_t;

    yaw_p = yaw_p + 0.5*noise_yr*delta_t*delta_t;
    yawd_p = yawd_p + noise_yr*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  cout << "Updating Lidar" << endl;
  // Complete this function! Use lidar data to update the belief about the object's
  // position. Modify the state vector, x_, and covariance, P_.

  // You'll also need to calculate the lidar NIS.

  //measurement matrix
  int n_z = 2;
  MatrixXd H_ = MatrixXd(n_z, n_x_);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  VectorXd z = meas_package.raw_measurements_;

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(n_x_,n_x_);
  P_ = (I - K * H_) * P_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  cout << "Updating Radar" << endl;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    // extract values for better readability
    double x,y,v,yaw,v1,v2;
    x = Xsig_pred_(0,i);
    y = Xsig_pred_(1,i);
    v  = Xsig_pred_(2,i);
    yaw = Xsig_pred_(3,i);
    
    v1 = cos(yaw)*v;
    v2 = sin(yaw)*v;
    
    // measurement model
    Zsig(0,i) = sqrt(x*x + y*y);                    //r
    Zsig(1,i) = atan2(y,x);                         //phi
    Zsig(2,i) = (x*v1 + y*v2 ) / sqrt(x*x + y*y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S += R_;

  // Create vector for incomming measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 sigma points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
}
