#include <iostream>
#include "ukf.h"

UKF::UKF() {
    //TODO Auto-generated constructor stub
    Init();
}

UKF::~UKF() {
    //TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/

void UKF::SigmaPointPrediction(MatrixXd *Xsig_out) {

    //set state dimension
    int n_x = 5;

    //set augmented dimension
    int n_aug = 7;

    //create example sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
    Xsig_aug <<
             5.7441, 5.85768, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.63052, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441,
            1.38, 1.34566, 1.52806, 1.38, 1.38, 1.38, 1.38, 1.38, 1.41434, 1.23194, 1.38, 1.38, 1.38, 1.38, 1.38,
            2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.2049, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.2049, 2.2049,
            0.5015, 0.44339, 0.631886, 0.516923, 0.595227, 0.5015, 0.5015, 0.5015, 0.55961, 0.371114, 0.486077, 0.407773, 0.5015, 0.5015, 0.5015,
            0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.3528, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.3528,
            0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641, 0,
            0, 0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641;

    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

    double delta_t = 0.1; //time diff in sec
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

    int n_points = 15;

    for (size_t p = 0; p < n_points; p++) {
        auto sig = Xsig_aug.col(p);
        double px = sig(0);
        double py = sig(1);
        double v = sig(2);
        double yaw = sig(3);
        double yawd = sig(4);
        double nu_a = sig(5);
        double nu_yawdd = sig(6);

        VectorXd x(5);
        x << px, py, v, yaw, yawd;

        VectorXd a(5);

        if (abs(yawd) > 0.001) {
            a.row(0) << v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            a.row(1) << v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
        } else {
            a.row(0) << v * cos(yaw) * delta_t;
            a.row(1) << v * sin(yaw) * delta_t;
        }
        a.row(2) << 0;
        a.row(3) << yawd * delta_t;
        a.row(4) << 0;

        VectorXd b(5);
        b.row(0) << .5 * delta_t * delta_t * cos(yaw) * nu_a;
        b.row(1) << .5 * delta_t * delta_t * sin(yaw) * nu_a;
        b.row(2) << delta_t * nu_a;
        b.row(3) << .5 * delta_t * delta_t * nu_yawdd;
        b.row(4) << delta_t * nu_yawdd;

        Xsig_pred.col(p) << x + a + b;
    }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

    //print result
    std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

    //write result
    *Xsig_out = Xsig_pred;

}
