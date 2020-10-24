//
// Created by ouc on 2020/4/21.
//

#ifndef NEWUNDERWATERSFM_UNDERWATERSFM_H
#define NEWUNDERWATERSFM_UNDERWATERSFM_H

#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;
using namespace Eigen;
/*
 * underwaterSFM类，
 * getUnderwaterImg()水下成像、
 * backProjRay()求解rw
 * findMatrixG()估计g
 * recoverPose()分解r/t
 * getImgPoint() 取点
 *
 * vector<Vector3d> vRealPoint
 * vector<Vector3d> vEstiPoint
 * vector<Vector2d> vImg1Point, vImg2Point
 * real_nw, ng, na, backProj_nw; dw, dg, da;
 * Matrix3d _R, Vector3d _t
 */
class UnderwaterSFM {
public:
    UnderwaterSFM(){};
    ~UnderwaterSFM(){}

    //需要初始化的参数只能在构造函数里赋值吗
    /*
     * 构造函数的输入参数：
     * double n[3]:空气、玻璃、水的折射率，na，ng，nw
     * double d[3]：相机到玻璃板距离、玻璃板厚度、待求点到玻璃板的距离
     * Matrix& R: 图像对应的相机位姿（wc）
     * Vector3d& t：
     * Matrix3d K: 相机内参
     * bool isScaleMeter: 位姿和空间点坐标的尺度是否才哦嗯米制
     * bool is_nw_error：折射率nw是否包含误差
     * int pointscale：取点间隔
     */
    UnderwaterSFM(double n[3], double d[3], Matrix3d& R, Vector3d& t, Matrix3d K, bool isScaleMeter, bool is_nw_error, int pointscale, double nw_error = 0):
            _real_nw(n[0]), _ng(n[1]), _na(n[2]), _dg(d[1]), _da(d[2]),
            _real_R(R), _real_t(t), _K(K),
            _isScaleMeter(isScaleMeter), _is_nw_error(is_nw_error), _pointScale(pointscale),
            _nw_error(nw_error)
            {
        cout<<"construction succedd!"<<endl;
            };

    void runUnderwaterSFM(String outputfilenametxt, Mat colorImg, Mat depthImg);

    /*
     * 获取水下图像
     * 输入：原彩色图及深度图，3d点，图1中的2d点，图2中的2d点?没有写好
     * 输出：
     */
    void getUnderwaterImg(Mat colorImg, Mat depthImg);
    void Img2match(Mat img1, Mat img2);


    void underwaterSFM();

    /*
     * 反向投影
     */
    void backProjRay();

    void mytriangulation();

    void myrecoverPose(Matrix3d R, Vector3d t);

    void reconstructPoint(vector<Vector3d> vEstiPoint);

    void saveUnderwaterImg(Mat colorImg, vector<Vector2d> point_origin2,
            bool isShow, String filename1, String filename2);

    /*
     * 输出点云
     * 参数1：彩色图像; 参数2：图像点3D坐标vecotr；参数3：图像点2D坐标vector；参数4：保存的ply文件名(带后缀)
     */
    void savePointCloud(Mat colorImg, vector<Vector3d> p3d, vector<Vector2d> p2d, String filename);
    double _error_p_mean;
    double _error_pigrefrax;
    double _error_t_norm;
    double _error_t;

private:
    /*
     * 求解点Pc的theta值并得到对应水下成像点projectpoint
     */
    void ceres_solve_theta(Vector3d& Pc, Vector2d& projectpoint);
    void pose_estimation_2d2d(vector<Vector2d> vpoints1, vector<Vector2d> vpoints2, Mat& R_ignore, Mat& t_ignore);
    vector<Vector3d> _vRealPoint, _vEstiPoint, _vIgnoreRefraxPoint, _vrw1, _vrw2, _vpout1, _vpout2;
    vector<Vector2d> _vImg1Point, _vImg2Point;
    double _nw_error;
    double _real_nw = 1.33, _backProj_nw = _real_nw+_nw_error;
    double _ng, _na;
    double _dg, _da;
    double _theta_initial = 0.15;
    int _pointScale;    //取点间隔
    bool _isScaleMeter;
    bool _is_nw_error;
    Matrix3d _real_R, _esti_R;    //Rwc, from camera coordinate to world coordinate
    Vector3d _real_t, _esti_t;    //twc
    Matrix3d _K;

    Mat _img_proj1, _img_proj2;



};


#endif //NEWUNDERWATERSFM_UNDERWATERSFM_H
