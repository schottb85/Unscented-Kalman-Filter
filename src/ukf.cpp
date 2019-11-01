#include "ukf.h"
#include "Eigen/Dense"
#include "iostream"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8; // (!) optimize based on NIS

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.55; // (!) optimize based on NIS
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;
  
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  double weight = 0.5/(lambda_+n_aug_);
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = weight_0;

  for (int i=1; i<2*n_aug_+1; ++i) {  
    weights_(i) = weight;
  }

  Xsig_pred_ = MatrixXd(5, 15);
  Xsig_pred_.fill(0.0);

  double NIS_lidar_ = 0.;
  double NIS_radar_ = 0.;
}

UKF::~UKF() {
  //std::cout << "!!!!! DELETE INSTANCE !!!!!" << std::endl;
}

void UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {
  if(meas_package.sensor_type_ == meas_package.LASER && use_laser_){
     ProcessLidarMeasurement(meas_package);
  }
  else if(meas_package.sensor_type_ == meas_package.RADAR && use_radar_){
     ProcessRadarMeasurement(meas_package);
  }
}

void UKF::ProcessLidarMeasurement(const MeasurementPackage& meas_package){
  if(!is_initialized_){
    InitLidarMeasurement(meas_package);
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);
  UpdateLidar(meas_package);
}

void UKF::ProcessRadarMeasurement(const MeasurementPackage& meas_package){
  if(!is_initialized_){
    InitRadarMeasurement(meas_package);
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(dt);
  UpdateRadar(meas_package);
}

void UKF::InitLidarMeasurement(const MeasurementPackage& meas_package){
  x_.fill(0.);
  x_[0] = meas_package.raw_measurements_[0]; // px
  x_[1] = meas_package.raw_measurements_[1]; // py

  // initialize covariance matrix with diagonal matrix
  P_.fill(0.);

  P_(0, 0) = std_laspx_*std_laspx_;
  P_(1, 1) = std_laspy_*std_laspy_;
  P_(2, 2) = std_radr_*std_radr_;
  P_(3, 3) = std_radphi_*std_radphi_;
  P_(4, 4) = std_radrd_*std_radrd_;
}

void UKF::InitRadarMeasurement(const MeasurementPackage& meas_package){
  x_.fill(0.);

  double rho = meas_package.raw_measurements_[0];      // range: radial distance from origin
  double phi = meas_package.raw_measurements_[1];      // bearing: angle between rho and x axis
  //double rho_dot = meas_package.raw_measurements_[2];  // radial velocity: change of rho
  double px = cos(phi) * rho;
  double py = sin(phi) * rho;

  x_[0] = px; // px
  x_[1] = py; // py

  // initialize covariance matrix with diagonal matrix
  P_.fill(0.);

  P_(0, 0) = std_laspx_*std_laspx_;
  P_(1, 1) = std_laspy_*std_laspy_;
  P_(2, 2) = std_radr_*std_radr_;
  P_(3, 3) = std_radphi_*std_radphi_;
  P_(4, 4) = std_radrd_*std_radrd_;
}


void UKF::GenerateAugmentedSigmaPoints(MatrixXd* Xsig_out){
  // create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_; // mean of a and yaw_dd is both 0

  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_*std_a_; // add variance of longitudinal acc
  P_aug(n_x_+1, n_x_+1) = std_yawdd_*std_yawdd_; // add variance of yaw rate derivative

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i=0; i<n_aug_; ++i){
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
      Xsig_aug.col(n_aug_+i+1) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i); 
  }

  // write result
  *Xsig_out = Xsig_aug;
}

void UKF::PredictSigmaPoints(const MatrixXd& Xsig_aug, MatrixXd* Xsig_out, double dt) {
  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // predict sigma points
  // avoid division by zero
  // write predicted sigma points into right column
  for(int i=0; i<2*n_aug_+1; ++i) //loop sigma points
  {
      double px = Xsig_aug(0,i);
      double py = Xsig_aug(1,i);
      double v = Xsig_aug(2,i);
      double psi = Xsig_aug(3,i);
      double psid = Xsig_aug(4,i);
      double noise_a = Xsig_aug(5,i);
      double noise_psidd = Xsig_aug(6,i);

      if(fabs(psid)<1e-04)
      {
          // linear simplified model
          Xsig_pred.col(i)(0) = px + v*cos(psi)*dt;
          Xsig_pred.col(i)(1) = py + v*sin(psi)*dt;
      }
      else{
          Xsig_pred.col(i)(0) = px + v/psid*(sin(psi+psid*dt)-sin(psi));
          Xsig_pred.col(i)(1) = py + v/psid*(-cos(psi+psid*dt)+cos(psi));
      }
      Xsig_pred.col(i)(0) += 0.5*dt*dt*cos(psi)*noise_a;
      Xsig_pred.col(i)(1) += 0.5*dt*dt*sin(psi)*noise_a;
      Xsig_pred.col(i)(2) = v + dt*noise_a;
      Xsig_pred.col(i)(3) = psi + psid*dt + 0.5*dt*dt*noise_psidd;
      Xsig_pred.col(i)(4) = psid + dt*noise_psidd;
  }

  // write result
  *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(const MatrixXd& Xsig_pred, VectorXd* x_out, MatrixXd* P_out) {
  // create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  // create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  // predicted state mean
  x.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    x = x + weights_(i) * Xsig_pred.col(i);
  }

  // predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  // write result
  *x_out = x;
  *P_out = P;
}


void UKF::PredictMeasurement(int n_z, const MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, const MatrixXd &R){
  // mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // innovation covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  S = S + R;
}

void UKF::UpdateState(const VectorXd& z, const VectorXd& z_pred, const MatrixXd& S, const MatrixXd& Zsig) {
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, z.size());

  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred; // z_pred_ = (r, phi, rdot)

  // angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
}



void UKF::Prediction(double delta_t) {
  // 1. Setup the sigma points to prepare the prediction step (augmented by noise components)
  // 2. Predict Sigma points
  // 3. Predict mean and covariance
  MatrixXd Xsig_aug = MatrixXd(7, 15);
  GenerateAugmentedSigmaPoints(&Xsig_aug);

  PredictSigmaPoints(Xsig_aug, &Xsig_pred_, delta_t);

  PredictMeanAndCovariance(Xsig_pred_, &x_, &P_);
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  const VectorXd& z = meas_package.raw_measurements_; // marker.rho, marker.phi, marker.rho_dot

 // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;                       // px
    Zsig(1,i) = p_y;                       // py
  }

  MatrixXd R = MatrixXd(n_z,n_z);
      // add measurement noise covariance matrix
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.);

  // fill S and z_pred
  PredictMeasurement(n_z, Zsig, z_pred, S, R);

  // update x, P based on lidar measurement data
  UpdateState( z, z_pred, S, Zsig);

  NIS_lidar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  std::cout << "NIS level of 7.8: Lidar: " << NIS_radar_ << std::endl;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  const VectorXd& z = meas_package.raw_measurements_; // marker.rho, marker.phi, marker.rho_dot

 // set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_radar = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar, 2 * n_aug_ + 1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }

  MatrixXd R = MatrixXd(n_z_radar,n_z_radar);
      // add measurement noise covariance matrix
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0,std_radrd_*std_radrd_;

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar);
  
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
  S.fill(0.);

  // fill S and z_pred
  PredictMeasurement(n_z_radar, Zsig, z_pred, S, R);

  // update based on radar measurement data
  UpdateState(z, z_pred, S, Zsig);

  NIS_radar_ = (z - z_pred).transpose() * S.inverse() * (z - z_pred);

  std::cout << "NIS level of 7.8: Radar: " << NIS_radar_ << std::endl;
}