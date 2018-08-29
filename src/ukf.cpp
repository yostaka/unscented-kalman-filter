#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  P_ <<
    1, 0, 0, 0, 0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2 * n_aug_ + 1);

  // set weights
  int na = 2 * n_aug_ + 1;
  float w1 = lambda_ / (lambda_ + n_aug_);
  float wn = 1.0 / (2.0 * (lambda_ +  n_aug_));

  for(int i=0; i<na; i++) {
      if(i==0) {
          weights_(i) = w1;
      } else {
          weights_(i) = wn;
      }
  }

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

  // Initialization
  if (!is_initialized_) {

      // first measurement
      x_ = VectorXd(5);

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          cout << "Initialize state vector x with Radar measurement" << endl;
          double ro = meas_package.raw_measurements_[0];
          double theta = meas_package.raw_measurements_[1];
          double ro_dot = meas_package.raw_measurements_[2];

          double px = ro * cos(theta);
          double py = ro * sin(theta);
          double v = ro_dot;
          x_ << px, py, v, 0, 0;
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          cout << "Initialize state vector x with Laser measurement" << endl;
          double px = meas_package.raw_measurements_[0];
          double py = meas_package.raw_measurements_[1];
          x_ << px, py, 0, 0, 0;
      }

      // previous timestamp
      time_us_ = meas_package.timestamp_;

      // done initializing, no need to predict or update
      is_initialized_ = true;
      return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  long long previous_timestamp = time_us_;
  long long current_timestamp = meas_package.timestamp_;
  double dt = (current_timestamp - previous_timestamp) / 1000000.0;

  Prediction(dt);

  time_us_ = current_timestamp;



  /*****************************************************************************
   *  Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
  }
  else {
      UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << endl << x_ << endl;
  cout << "P_ = " << endl << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // Generate sigma points
  // Xsig: 7 x 15 Matrix
  MatrixXd Xsig = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  GenerateAugmentedSigmaPoints(&Xsig);

  // Predict sigma points
  // Xsig_pred_: 5 x 15 Matrix
  Xsig_pred_.fill(0.0);
  SigmaPointPrediction(Xsig, &Xsig_pred_, delta_t);

  // Calculate predicted mean and covariance
  // x_pred: 5 length vector
  // P_pred: 5 x 5 Matrix
  VectorXd x_pred = VectorXd(n_x_);
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  PredictMeanAndCovariance(&x_pred, &P_pred);
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

    int n_z = 2;

    MatrixXd H = MatrixXd(n_z, 5);
    H << 1, 0, 0, 0, 0,
         0, 1, 0, 0, 0;

    VectorXd z = VectorXd(n_z);
    z <<
        meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1];

    MatrixXd R = MatrixXd(2, 2);
    R <<
        std_laspx_ * std_laspx_, 0,
        0, std_laspy_ * std_laspy_;

    VectorXd z_pred = H * x_;
    VectorXd y = z - z_pred;
    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P_ * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // new estimate
    // x_ is new state matrix expressing most likely state (2d position and 2d velocity)
    // P_ is new covarience matrix expressing uncertainty of the state
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H) * P_;
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
    int n_z = 3;

    // Calculate mean predicted measurement z_pred and predicted measurement covariance matrix S_pred
    // in radar measurement space [rho, phi, rho_dot]
    VectorXd z_pred = VectorXd(n_z);
    MatrixXd S_pred = MatrixXd(n_z, n_z);
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    PredictRadarMeasurement(&z_pred, &S_pred, &Tc);

    // Update state vector x_ and covariance matrix P_ with radar measurement
    VectorXd z = VectorXd(n_z);
    z <<
        meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1],
        meas_package.raw_measurements_[2];

    UpdateStateWithRadarMeasurement(z, z_pred, S_pred, Tc);

}

void UKF::GenerateAugmentedSigmaPoints(MatrixXd *Xsig_out) {

    // create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.fill(0.0);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
    Xsig_aug.fill(0.0);

    // create augmented mean state
    x_aug << x_, 0, 0;
    Xsig_aug.col(0) << x_aug;

    // create augmented covariance matrix
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    // calculate square root of P
    MatrixXd A = P_aug.llt().matrixL();

    // create augmented sigma points
    for (int i=0; i<n_aug_; i++) {
        Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
        Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
    }

    *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig, MatrixXd* Xsig_out, double delta_t) {

    MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
    Xsig_pred.fill(0.0);

    // predict sigma points
    for (int i=0; i<2*n_aug_+1; i++) {
        VectorXd Statek = VectorXd(n_aug_);
        Statek << Xsig.col(i);

        float px = Statek(0);
        float py = Statek(1);
        float v  = Statek(2);
        float psi = Statek(3);
        float psi_dot = Statek(4);
        float nu_ak = Statek(5);
        float nu_psi_dotdot = Statek(6);

        if(psi_dot > 0.001) {
            VectorXd firstTerm = VectorXd(5);
            firstTerm <<
                      px, py, v, psi, psi_dot;

            VectorXd secondTerm = VectorXd(5);
            secondTerm <<
                       v/psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi)),
                    v/psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi)),
                    0,
                    psi_dot * delta_t,
                    0;

            VectorXd thirdTerm = VectorXd(5);
            thirdTerm <<
                      1.0/2.0 * delta_t * delta_t * cos(psi) * nu_ak,
                    1.0/2.0 * delta_t * delta_t * sin(psi) * nu_ak,
                    delta_t * nu_ak,
                    1.0/2.0 * delta_t * delta_t * nu_psi_dotdot,
                    delta_t * nu_psi_dotdot;

            Xsig_pred.col(i) << firstTerm + secondTerm + thirdTerm;
        }
        else {
            VectorXd firstTerm = VectorXd(5);
            firstTerm <<
                      px, py, v, psi, psi_dot;

            VectorXd secondTerm = VectorXd(5);
            secondTerm <<
                       v * cos(psi) * delta_t,
                    v * sin(psi) * delta_t,
                    0,
                    psi_dot * delta_t,
                    0;

            VectorXd thirdTerm = VectorXd(5);
            thirdTerm <<
                      1.0/2.0 * delta_t * delta_t * cos(psi) * nu_ak,
                    1.0/2.0 * delta_t * delta_t * sin(psi) * nu_ak,
                    delta_t * nu_ak,
                    1.0/2.0 * delta_t * delta_t * nu_psi_dotdot,
                    delta_t * nu_psi_dotdot;

            Xsig_pred.col(i) << firstTerm + secondTerm + thirdTerm;
        }

    }

    *Xsig_out = Xsig_pred;
}

void UKF::PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred) {

    VectorXd x = VectorXd(n_x_);
    MatrixXd P = MatrixXd(n_x_, n_x_);

    int na = 2 * n_aug_ + 1;

    // clear x and P
    x.fill(0.0);
    P.fill(0.0);

    // predict state mean
    for(int i=0; i<na; i++) {
        x += weights_(i) * Xsig_pred_.col(i);
    }

    // predict state covariance matrix
    for(int i=0; i<na; i++) {
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x;

        // angle normalization
        while (x_diff(3)> M_PI) { x_diff(3)-=2.*M_PI; }
        while (x_diff(3)<-M_PI) { x_diff(3)+=2.*M_PI; }

        P = P + weights_(i) * x_diff * x_diff.transpose();
    }

    *x_pred = x;
    *P_pred = P;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd *Tc_out) {

    int n_z = 3;

    // create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);

    // transform sigma points into measurement space
    for(int i=0; i<2*n_aug_+1; i++) {
        float px = Xsig_pred_(0, i);
        float py = Xsig_pred_(1, i);
        float v  = Xsig_pred_(2, i);
        float psi = Xsig_pred_(3, i);
//        float psi_dot = Xsig_pred_(4, i);

        Zsig.col(i) <<
                sqrt(px*px + py*py),
                atan2(py,px),
                (px*cos(psi)*v + py*sin(psi)*v) / sqrt(px*px + py*py);
    }

    //calculate mean predicted measurement
    for(int i=0; i<2*n_aug_+1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //calculate innovation covariance matrix S
    MatrixXd R = MatrixXd(n_z, n_z);
    R.fill(0.0);
    R <<
            std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    S.fill(0.0);
    for(int i=0; i<2*n_aug_+1; i++) {
        S = S + weights_(i) * (Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose();
    }
    S = S + R;

    // create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);

    for(int i=0; i<2*n_aug_+1; i++) {
        VectorXd x_diff = VectorXd(n_x_);
        x_diff = Xsig_pred_.col(i) - x_;

        VectorXd z_diff = VectorXd(n_z);
        z_diff = Zsig.col(i) - z_pred;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    *z_out  = z_pred;
    *S_out  = S;
    *Tc_out = Tc;
}

void UKF::UpdateStateWithRadarMeasurement(VectorXd z, VectorXd z_pred, MatrixXd S_pred, MatrixXd Tc) {

    int n_z = 3;

    // calculate Kalmain gain K
    MatrixXd K = MatrixXd(n_x_, n_z);
    K.fill(0.0);
    K = Tc * S_pred.inverse();

    // update state mean and covariance matrix
    x_ = x_ + K * (z - z_pred);
    P_ = P_ - K * S_pred * K.transpose();
}

