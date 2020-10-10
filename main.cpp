#include <iostream>
#include <Eigen/Geometry>
#include "UnderwaterSFM.h"

int main() {
    double n[3] = {1.1, 1.49, 1.0};
    double d[3] = {0, 40, 100};
    Vector3d eularAngle(-0.25*M_PI, -0.05*M_PI, 0.01*M_PI);
    AngleAxisd rollAngle(AngleAxisd(eularAngle(2), Vector3d::UnitX()));
    AngleAxisd pitchAngle(AngleAxisd(eularAngle(1), Vector3d::UnitY()));
    AngleAxisd yawAngle(AngleAxisd(eularAngle(0), Vector3d::UnitZ()));
    Matrix3d R;
    R = yawAngle*pitchAngle*rollAngle;
    Vector3d t = Vector3d(0,15,60);
    Matrix3d K;
    K<< 535.4, 0, 320.1,
        0, 539.2, 247.6,
        0, 0, 1;
    bool isScaleMeter = false;
    bool is_nw_error = true;
    int pointscale = 1000;
    Mat img = imread("./3.png");
    Mat img_d = imread("./3.pgm");
    vector<Vector3d> vRealPoint;
    vector<Vector2d> vImg1Point, vImg2Point;
    fstream errornw;
    String txtfilename = to_string(n[0]) +to_string(d[1])+ "underwater_nw_error_p.txt";
    errornw.open(txtfilename, ios::binary | ios::app | ios::in | ios::out);
    for(int i = 0; i<100; i++){
        double error_nw = -0.1+i/500.0;
        UnderwaterSFM myunderwater(n, d, R, t, K, isScaleMeter, is_nw_error, pointscale , error_nw);
//    myunderwater.getUnderwaterImg(img, img_d, vRealPoint);
//    myunderwater.backProjRay(p_vec1, p_vec2, vrw1, vrw2, vpout1, vpout2);
        String test = "test.txt";
        myunderwater.runUnderwaterSFM(test, img, img_d);
        double used_nw = n[0] + error_nw;
        double error_p_mean = myunderwater._error_p_mean;
        double error_pig_mean = myunderwater._error_pigrefrax;
        double error_t = myunderwater._error_t;
        double error_t_norm = myunderwater._error_t_norm;
        errornw<<fixed<<setprecision(8)<<used_nw<<"\t" <<error_p_mean<<"\t"<<error_pig_mean<<"\t"<<error_t<<"\t"<<error_t_norm<<endl;
    }
    errornw.close();
    return 0;
}