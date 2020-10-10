//
// Created by ouc on 2020/4/21.
//

#include "UnderwaterSFM.h"
#include <pcl/point_types.h>
#include <pcl/io/ply_io.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <opencv2/calib3d.hpp>

struct costFunctor
{
    double _theta;
    double _px;

    double _d[3];   //order:[dw, dg, ga]
    double _n[3];

    //use constructor to apply un-optimized variables
    costFunctor(double d[3], double n[3], double px)
    {
        _px = abs(px);
        for(int i = 0; i<3; i++)
        {
            _d[i] = d[i];
            _n[i] = n[i];
        }
    }

    template<typename T>
    bool operator()(const T* const theta, T* residual) const {      //is here the minimum of residual = 0? if it's bigger than 0, the solution of theta is wrong!
        residual[0] = _px -_d[0]*_n[2]*sin(theta[0])/(_n[0]*sqrt(T(1)-pow((_n[2]*sin(theta[0])/_n[0]),2)))
                      - _d[1]*_n[2]*sin(theta[0])/(_n[1]*sqrt(T(1)-pow((_n[2]*sin(theta[0])/_n[1]),2)))
                      -_d[2]*sin(theta[0])/sqrt(T(1)-(sin(theta[0]))*(sin(theta[0])));
        return true;
    }
};


struct Solveg
{
    double _p1[3];
    double _p2[3];
    double _d1;
    double _d2;
    Solveg(Vector3d p1, Vector3d p2, double d1, double d2):_d1(d1), _d2(d2)
    {
        for(int i = 0; i<3; i++)
        {
            _p1[i] = p1(i);
            _p2[i] = p2(i);
        }
    }

    template<typename T>
    bool operator()(const T* const g, T* residual) const
    {
        residual[0] =  T(_p1[0])*T(_p2[0])*g[0] + T(_p1[1])*T(_p2[0])*g[1] + T(_p1[2])*T(_p2[0])*g[2]
                       + T(_p1[0])*T(_p2[1])*g[3] + T(_p1[1])*T(_p2[1])*g[4] + T(_p1[2])*T(_p2[1])*g[5]
                       + T(_p1[0])*T(_p2[2])*g[6] + T(_p1[1])*T(_p2[2])*g[7] + T(_p1[2])*T(_p2[2])*g[8]
                       + (_d1*T(_p1[1])*T(_p2[0]) + _d2*T(_p1[0])*T(_p2[1]))*g[9] + (-_d1*T(_p1[0])*T(_p2[0]) + _d2*T(_p1[1])*T(_p2[1]))*g[10] + (_d2*T(_p1[2])*T(_p2[1]))*g[11]
                       + (_d1*T(_p1[1])*T(_p2[1]) - _d2*T(_p1[0])*T(_p2[0]))*g[12] + (-_d1*T(_p1[0])*T(_p2[1]) - _d2*T(_p1[1])*T(_p2[0]))*g[13] + (-_d2*T(_p1[2])*T(_p2[0]))*g[14]
                       +_d1*T(_p1[1])*T(_p2[2])*g[15] -_d1*T(_p1[0])*T(_p2[2])*g[16];
        residual[1] = g[9]*g[9] + g[10]*g[10] + g[11]*g[11] - T(1.0);
        residual[2] = g[12]*g[12] + g[13]*g[13] + g[14]*g[14] - T(1.0);
        return true;
    }
};


void UnderwaterSFM::getUnderwaterImg(Mat colorImg, Mat depthImg){
    //读取图像中的有深度信息的点，
    Mat img = colorImg;
    Mat img_depth = depthImg;
//    String outputnum_string =
    String outputnum_string = to_string(int(_da))+ to_string(int(_dg))+to_string(_real_nw); //玻璃厚度及距离做输出文档名称
    double depthscale;
    if(_isScaleMeter)
        depthscale = 5*1000;
    else
        depthscale = 5*1.000;
    int rows = img.rows, cols = img.cols;
//    Mat img_projection, img_projection2;
    _img_proj1 = Mat(rows, cols,CV_8UC3,Scalar(255,255,255));
    _img_proj2 = Mat(rows, cols,CV_8UC3,Scalar(255,255,255));
    vector<Vector3d> vptc, vptw;    //point_cvec单位位姿/第一帧;point_wvecR,t位姿/第二帧
    vector<Vector2d> point_origin;  //

    //read the 3D points
    //这里获取深度存在的水下真实3d点（2类坐标系下的）及对应2d坐标
    for(int u = 0; u<cols; u++)      //pay attention to the (u,v) = (640,480)(u:colums, v:rows)傻逼
        for(int v = 0; v<rows; v++)
        {
            unsigned int d = img_depth.ptr<unsigned short>(v)[u];
            if( d == 0)
                continue;
            Vector3d point, point_w, point_c, p_homo;
            point[2] = double(d)/depthscale;
//            point[0] = point[2]*(u-cx)/fx;
//            point[1] = point[2]*(v-cy)/fy;
            p_homo<<u, v, 1;
            point = point[2]*_K.inverse()*p_homo;   //与注释行等价，但速度可能更慢一点
            vptc.push_back(point);  //相机坐标系下点坐标(单位位姿)
            point_w = _real_R*point + _real_t;
            vptw.push_back(point_w);    //世界系下点坐标（给定R，T位姿）
            point_origin.push_back(Vector2d(u,v));   //原始图像的二维点point(u,v)->u, |v
        }
    cout<<"size of all point with non-zero depth is "<<vptw.size()<<endl;

    //按一定规模取点
    vector<Vector2d> point_origin2; //相对原始二维点而言，去掉了在水下相机中不能成像的点
    vector<Vector3d> Pvec;  //3d坐标真值
    cout<<"original 3D points are: "<<endl;
    for(int i = 0; i<(vptw.size()/_pointScale); i++)
    {

        Vector2d p, p2;
        Vector3d P2 = vptc[_pointScale*i];      //P2的输入3d点是相机坐标系下的
        Vector3d P = vptw[_pointScale*i];       //P的输入3d点是世界坐标系下的，相当于给了相机一个Rwc
        //size equals to that of later estimated points?
        if(P[2]<_dg) p = Vector2d(-1,-1);
        else if(P2[2]<_dg) p2 = Vector2d(-1,-1);
        else{
            ceres_solve_theta(P, p);   //3D coordinate P is with camera coordinate
            ceres_solve_theta(P2, p2); //2d点的颜色丢失了吗,这里的点是坐标，不是图像，所以本身就不带颜色
        }
        //这样删点就很愉快了！（不是删，而是不加
        if (p(0) < 0 || p(1) < 0 || p2(0) < 0 || p2(1) < 0 ||
        p(0) > cols || p(1) > rows || p2(0) > cols || p2(1) > rows)
        {
            //cout<<"Point is out of view"<<"p_c: "<<p2.x<<", "<<p2.y<<"\tp_w: "<<p.x<<", "<<p.y<<endl;
            continue;
        }
        else {
            Pvec.push_back(P);
            _vRealPoint.push_back(P);
            point_origin2.push_back(point_origin[_pointScale*i]);
            _vImg1Point.push_back(p);     //p_worldv含旋转平移
            _vImg2Point.push_back(p2);   //p_camv
            //cout<<P2.transpose()<<endl;
        }
    }
    String real_pc = "real_pc"+outputnum_string;
    savePointCloud(img, Pvec, point_origin2, real_pc);  //这里传递的参数是否对应？
    cout<<"size of p_vec = "<<_vImg1Point.size()<<endl;

    cout<<"imax = "<<vptw.size()/_pointScale<<endl;
    ofstream outfile3d;
    String outputname1 = "real3dpoint_200"+outputnum_string+".txt";
    outfile3d.open(outputname1, ios::binary | ios::trunc | ios::in | ios::out);

    //输出两张图像 (save and display)
    String imgname = "underwater_w"+outputnum_string+".png";
    String imgname2 = "underwater_c"+outputnum_string+".png";
    //这里保存的是颜色信息，同一i，原图中point_origin2[i]（2是取了部分点的1）、对应第一帧中的p_vec[i]（3D点是世界坐标）、第二帧中的p_vec2[i]（3D点是Rcw之后的坐标）
    saveUnderwaterImg(img, point_origin2, false, imgname, imgname2);
    /*
    //save the image points, but why only save image2?
    ofstream fout;
    fout.open("pvec2.txt", ios::binary | ios::trunc | ios::in | ios::out);
    for(int i = 0; i<p_vec.size(); i++)
    {
        if(runoutput)
            fout<<fixed << setprecision(8)<<p_vec2[i].x<<" "<<p_vec2[i].y<<endl;
    }
    fout.close();
*/
}


void UnderwaterSFM::underwaterSFM(){
    backProjRay();


}

void UnderwaterSFM::ceres_solve_theta(Vector3d &Pc, Vector2d &projectpoint)
{

    double px = sqrt(Pc(0)*Pc(0) + Pc(1)*Pc(1));
    double pz = Pc(2);
    double dw = pz - _dg - _da;
    if(Pc(2)<_dg){
        cout<<"3D coordinates: "<<Pc.transpose()<<endl;
        cout<<"dw= "<<dw<<endl;
        cout<<"dw<0, 无法水下成像"<<endl;
    }
    double _px = abs(px);
    double theta = _theta_initial;
    double d[3] = {dw, _dg, _da};
    double n[3] = {_real_nw, _ng, _na};
    ceres::Problem problem;
    costFunctor* c_functor = new costFunctor(d, n, _px);
    ceres::CostFunction* c_function = new ceres::AutoDiffCostFunction<costFunctor, 1, 1>(c_functor);

    problem.AddResidualBlock(c_function, NULL, &theta);

    ceres::Solver::Options option;
//    option.linear_solver_type = ceres::DENSE_QR;  //why dont use?
//    option.minimizer_progress_to_stdout = true;

    problem.SetParameterLowerBound(&theta, 0, 0.0);
    problem.SetParameterUpperBound(&theta, 0, M_PI/2-0.000001);
    ceres::Solver::Summary summary;
    ceres::Solve(option, &problem, &summary);
    /*
    double error = d[0]*n[2]*sin(theta)/(n[0]*sqrt(1-pow(n[2]*sin(theta)/n[0],2)))
                   + d[1]*n[2]*sin(theta)/(n[1]*sqrt(1-pow(n[2]*sin(theta)/n[1],2)))
                   + d[2]*sin(theta)/sqrt(1-pow(sin(theta),2))
                   - px;
    if(debug_theta_error) cout<<"error = "<<fixed << setprecision(8)<<error<<endl;  //nearly equals 0!
    if(debug_theta_error) cout<<"theta = "<<fixed << setprecision(8)<<theta<<endl;
    */
    cout<<summary.BriefReport()<<endl;
    Vector3d qin;
    //这里的qin[0]没有乘cosfai，是否需要乘？如何乘？//已经乘进去了,Pc[0]/px
    qin(0) = d[2]*tan(theta)*Pc(0)/px; //here inverted can make the image positive, and use Pc[0]/px is right, cause px =Pc[0]
    qin(1) = Pc(1)*qin(0)/Pc[0];
    qin(2) = d[2];
    //solve the project point coordinate u,v
    Vector3d x = _K*qin;  //why inverted?
    projectpoint<<x(0)/x(2), x(1)/x(2);
    //cout<<"project point: "<<projectpoint.x<<projectpoint.y<<endl<<endl;
}


void UnderwaterSFM::savePointCloud(Mat colorImg, vector<Vector3d> p3d, vector<Vector2d> p2d, String filename){
    //这里同时保存水下真实的点云
    typedef pcl::PointXYZRGB PointT;
    typedef pcl::PointCloud<PointT> PointCloud;
    PointCloud::Ptr pointCloud_gt(new PointCloud);
    for(int i = 0; i<p3d.size(); i++)      //pay attention to the (u,v) = (640,480)(u:colums, v:rows)
    {
        Vector3d P2 = p3d[i];
        Vector2d p2 = p2d[i];
        PointT pc;
        pc.x = P2[0];
        pc.y = P2[1];
        pc.z = P2[2];
        pc.b = colorImg.at<Vec3b>(p2(0), p2(1))[0];
        pc.g = colorImg.at<Vec3b>(p2(0), p2(1))[1];
        pc.r = colorImg.at<Vec3b>(p2(0), p2(1))[2]; //这里是u,v还是v,y?
        pointCloud_gt->points.push_back(pc);
        //cout<<P2.transpose()<<endl;
    }
    pointCloud_gt->is_dense = false;
    cout<<"there are "<<pointCloud_gt->size()<<" points in the groundtruth cloud"<<endl;
//    pcl::io::savePLYFile(to_string(outputnum)+to_string(pointscale)+string("_groundtruth.ply"),*pointCloud_gt);
    pcl::io::savePLYFile(filename+".ply",*pointCloud_gt);
    cout<<"size of all point with non-zero depth in point cloud is "<<p3d.size()<<endl;
}


void UnderwaterSFM::saveUnderwaterImg(Mat colorImg, vector<Vector2d> point_origin2,
        bool isShow, String filename1, String filename2){
    for(int i = 0; i<_vImg1Point.size(); i++)
    {
        Point2d pveci = Point2d(_vImg1Point[i](0), _vImg1Point[i](1));
        Point2d pveci2 = Point2d(_vImg2Point[i](0), _vImg2Point[i](1));
        Point2d porigin2i = Point2d(point_origin2[i](0), point_origin2[i](1));
        _img_proj1.at<Vec3b>(pveci)[0] = colorImg.at<Vec3b>(porigin2i)[0];
        _img_proj1.at<Vec3b>(pveci)[1] = colorImg.at<Vec3b>(porigin2i)[1];
        _img_proj1.at<Vec3b>(pveci)[2] = colorImg.at<Vec3b>(porigin2i)[2];

        _img_proj2.at<Vec3b>(pveci2)[0] = colorImg.at<Vec3b>(porigin2i)[0];
        _img_proj2.at<Vec3b>(pveci2)[1] = colorImg.at<Vec3b>(porigin2i)[1];
        _img_proj2.at<Vec3b>(pveci2)[2] = colorImg.at<Vec3b>(porigin2i)[2];
    }
    vector<int> compression_params;
    compression_params.push_back(CV_IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    imwrite(filename1+".png",_img_proj1, compression_params);
    imwrite(filename2+".png",_img_proj2, compression_params);
    if(isShow){
        imshow(filename1, _img_proj1);
        waitKey(0);
        imshow(filename2, _img_proj2);
        waitKey(0);
    }
}


Point2f pixel2cam(const Point2d& p, const Matrix3d& K){  //this is on normalized plane
    return Point2f((p.x - K(0,2))/K(0,0),
                   (p.y - K(1,2))/K(1,1));
}


void solvedp(const Point2d& p2d, const Matrix3d& K, const double n[3], const double d[3],
             Vector3d& r_out, double& dd, Vector3d& p)
{
    double n1 = n[2];
    double n2 = n[1];
    double n3 = n[0];
    double l = d[2];
    double w = d[1];
    Point2f p_norm(pixel2cam(p2d, K));
    Vector3d xi;
    xi<<p_norm.x, p_norm.y, 1;
    Vector3d N(0,0,1);                              //undefine
    Vector3d r_in = xi.normalized();    //仅需单位化rin，r_mid 和r_out都直接是单位向量
    double cos_theta1 = r_in(2);
    double sin_theta1 = sqrt(r_in(0)*r_in(0)+r_in(1)*r_in(1));

    double p1 = n1/n2;
    double q1 = -p1*cos_theta1+sqrt(1-pow(p1*sin_theta1,2));
    Vector3d r_mid = p1*r_in + q1*N;
//    cout<<"norm of r_mid is "<<sqrt(r_mid(0)*r_mid(0) + r_mid(1)*r_mid(1) + r_mid(2)*r_mid(2))<<endl;//whether it's unit vector?
    double cos_theta2 = p1*cos_theta1 + q1; //(13)
//    double cos_theta2_1 = sqrt(1-pow(p1*sin_theta1,2));   //(14)
    double sin_theta2 = sqrt( pow(p1*r_in(0),2) + pow(p1*r_in(1),2) );//(11)
//    cout<<"costheta2^2 + sintheta2^2  = "<<cos_theta2*cos_theta2 + sin_theta2*sin_theta2<<endl;//whether ()=1?

    double p2 = n2/n3;  //这里出问题了
    double q2 = -p2*cos_theta2+sqrt(1-pow(p2*sin_theta2,2));
    double sin_theta3 = p2*sin_theta2;

    r_out = p2*r_mid + q2*N;
//    cout<<"rout = "<<r_out.transpose()<<"\t norm = "<<(r_out(0)*r_out(0)+r_out(1)*r_out(1)+r_out(2)*r_out(2))<<endl;
    dd = l+w-(l*sin_theta1/cos_theta1+w*sin_theta2/cos_theta2)*sqrt(1-sin_theta3*sin_theta3)/sin_theta3;
    p = l*r_in/(N.transpose()*r_in) + w*r_mid/(N.transpose()*r_mid);    //world or camera? we need camera, but it seems like world
}


void recoverRt(const double g[17], Matrix3d& R, Vector3d& t)
{
    Matrix3d RR;
    Vector3d r1, r2, r3, r2rec;
    r1<<g[9], g[10], g[11];
    r2<<g[12], g[13], g[14];    //in this way whether r1*r2=0?

    double r12 = r1.transpose()*r2; //not exactly equals 0
    r2rec = r2 - r12*r1;
    r2 = r2rec.normalized();
    r3 = r1.cross(r2);  //need to verify whether it's the right usage of r1Xr2 A:yes!
    r3 = r3.normalized();
    RR<<r1, r2, r3;
    R = RR.transpose();
    cout<<"see whether R*R^T = I: "<<fixed << setprecision(8)<<R*R.transpose()<<endl;
    Matrix3d E;
    E<<g[0], g[1], g[2],
            g[3], g[4], g[5],
            g[6], g[7], g[8];
    Matrix3d T = R.transpose()*E;    //error is too large, the diagonal is not 0,0,0
//    Vector3d t2 = Sophus::SO3::vee(T);
    t(0) = T(2,1)*(abs(T(1,2)) + abs(T(2,1)))/(2.0*abs(T(2,1)));
    t(1) = T(0,2)*(abs(T(0,2)) + abs(T(2,0)))/(2.0*abs(T(0,2)));
    t(2) = T(1,0)*(abs(T(0,1)) + abs(T(1,0)))/(2.0*abs(T(1,0)));
//    t = Sophus::SO3::vee(T);
//    t<<-T(1,2), -T(2,0), T(1,0);  //which element should I recover t?
}



void UnderwaterSFM::pose_estimation_2d2d(vector<Vector2d> vpoints1, vector<Vector2d> vpoints2, Mat& R_ignore, Mat& t_ignore){
    Mat K = ( Mat_<double> ( 3,3 ) << _K(0,0), 0, _K(0,2), 0, _K(1,1), _K(1,2), 0, 0, 1 );
    Mat essential_matrix;
//    essential_matrix = findEssentialMat(points1, points2,focal_length,principle_pt);
//    recoverPose(essential_matrix, points1, points2, R_ignore, t_ignore, focal_length, principle_pt);
    vector<Point2d> points1, points2;
    for(int i = 0; i<vpoints1.size(); i++){
        Point2d p1 = Point2d(vpoints1[i](0), vpoints1[i](1));
        Point2d p2 = Point2d(vpoints2[i](0), vpoints2[i](1));
        points1.push_back(p1);
        points2.push_back(p2);
    }
    essential_matrix = findEssentialMat(points1, points2,K);
    recoverPose(essential_matrix, points1, points2, K, R_ignore, t_ignore);
    cout<<"R_ignore is "<<R_ignore<<endl<<"t_ignore is "<<t_ignore<<endl;

    //start to triangulate
    Mat T1 = (Mat_<double>(3,4)<<1,0,0,0,
            0,1,0,0,
            0,0,1,0);
    Mat T2 = (Mat_<double>(3,4)<<R_ignore.at<double>(0,0),R_ignore.at<double>(0,1),R_ignore.at<double>(0,2),t_ignore.at<double>(0,0),
            R_ignore.at<double>(1,0),R_ignore.at<double>(1,1),R_ignore.at<double>(1,2),t_ignore.at<double>(1,0),
            R_ignore.at<double>(2,0),R_ignore.at<double>(2,1),R_ignore.at<double>(2,2),t_ignore.at<double>(2,0));
    Mat pts_4d;
    //像素坐标转为相机坐标
    vector<Point2d> points_cam_1, points_cam_2;
    for(int i = 0; i<points1.size(); i++){
        Point2d cam_1 (( points1[i].x-K.at<double>(0,2) ) / K.at<double>(0,0),
                       ( points1[i].y-K.at<double>(1,2) ) / K.at<double>(1,1));
        Point2d cam_2 (( points2[i].x-K.at<double>(0,2) ) / K.at<double>(0,0),
                       ( points2[i].y-K.at<double>(1,2) ) / K.at<double>(1,1));
        points_cam_1.push_back(cam_1);
        points_cam_2.push_back(cam_2);
    }

    triangulatePoints(T1, T2, points_cam_1, points_cam_2, pts_4d);

    //transform the homogenous coordinate to non-homogenous coordinate
    for(int i = 0; i<pts_4d.cols; i++){
        Mat x = pts_4d.col(i);
        x = x/x.at<double>(3,0);
        Vector3d p3d = Vector3d(x.at<double>(0,0),
                    x.at<double>(1,0),
                    x.at<double>(2,0));
        _vIgnoreRefraxPoint.push_back(p3d);
    }
}

/*
 * 反向投影
 */
void UnderwaterSFM::backProjRay(){
    ceres::Problem problem;
    double g[17] = {0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0};
//    vector<Point2d> p2d1v, p2d2v;
    for(int i = 0; i<_vImg1Point.size(); i++){
        ////////for 2d point we estimate the 3d rout1 & d1, rout2 & d2
        Point2d p2d1 = Point2d(_vImg1Point[i](0), _vImg1Point[i](1));
        Point2d p2d2 = Point2d(_vImg2Point[i](0), _vImg2Point[i](1));
        Vector3d rout1, rout2, p_out1, p_out2;
        double d1, d2;
        double n[3] = {_backProj_nw, _ng, _na};
        double d[3] = {0, _dg, _da};
        //here must be trouble!
        //改变实际用的折射率情况，看错一点会增加多少误差
        solvedp(p2d1, _K, n, d, rout1, d1, p_out1);
        solvedp(p2d2, _K, n, d, rout2, d2, p_out2);
        ////////for 3d rout1,2 we estimate the g[17]
        Solveg* s_g = new Solveg(rout1, rout2, d1, d2);
        ceres::CostFunction* c_f = new ceres::AutoDiffCostFunction<Solveg,3,17>(s_g);
        problem.AddResidualBlock(c_f, nullptr, g);
        _vrw1.push_back(rout1); //这里把rout1,rout2, pout1, pout2入栈是要做啥？——————为了三角化方便
        _vrw2.push_back(rout2);
        _vpout1.push_back(p_out1);
        _vpout2.push_back(p_out2);
    }
    //here we need all points to calculate g;
    ceres::Solver::Options option;
    ceres::Solver::Summary summary;
    option.function_tolerance = 1e-16;
    option.gradient_tolerance = 1e-20;
    option.parameter_tolerance = 1e-18;
//    option.function_tolerance = 1.0e-16;
    option.max_num_iterations = 500;
    option.linear_solver_type = ceres::DENSE_QR;    //the problem is too large for DENSE_QR?
    option.minimizer_progress_to_stdout = true;
    ceres::Solve(option, &problem, &summary);
    cout<<summary.BriefReport()<<endl;
    cout<<"optimized g = ";
    for(auto k:g)
        cout<<fixed << setprecision(8)<<k<<" ";
    cout<<endl;

    Matrix3d Rcw;
    recoverRt(g, Rcw, _esti_t);//here the t is twc

    Matrix3d _esti_R = Rcw.transpose();  //Rwc   //it is for comparision
    Vector3d tcw = -_esti_R*_esti_t;

    cout<<"estimated Rcw = "<<fixed << setprecision(8)<<Rcw<<endl<<"estimated tcw = "<<tcw.transpose()<<endl; //it is for triangulation input
    Quaterniond rr(Rcw);
    cout<<"estimated R(quaternion) = "<<fixed << setprecision(8)<<rr.coeffs().transpose()<<endl;

    Quaterniond r_compare(_esti_R);
    cout<<"this R suit to campare with the groundtruth:"<<endl<<endl<<endl;
    cout<<"estimated Rwc = "<<fixed << setprecision(8)<<_esti_R<<endl<<"estimated twc = "<<_esti_t.transpose()<<endl;
    cout<<"estimater Rwc(quaternion) = "<<fixed << setprecision(8)<<r_compare.coeffs().transpose()<<endl<<endl<<endl;

    Vector3d t_error = _esti_t-_real_t;
    double error_t = sqrt(t_error.transpose()*t_error)/3.0;     //平移向量样除以3,而坐标误差不用除以3？
    cout<<fixed << setprecision(8)<<error_t<<endl;    //这是算的什么？

    //无折射时位姿估计结果，这里是否测试过？
    Mat R_ignore, t_ignore;
    pose_estimation_2d2d(_vImg1Point, _vImg2Point, R_ignore, t_ignore);    //R,t are from p_vec to p_vec2
    //pose_estimation_2d2d(p_vec, p_vec, R_ignore, t_ignore);
}


void UnderwaterSFM::mytriangulation(){
    ////////for R,t, 2d point we triangulate the depth
    /// a new triangulation method
//    vector<Vector3d> pv2;   //pv是考虑折射时，求解得到的3d坐标（第一帧下）
    ceres::Problem problemlamda;
    for(int i = 0; i<_vrw1.size(); i++){
        //firstly we check whether rout1& rout2 lie in common plane
        //translate the pout2& rout2 into world coordinate, pout is a point thus need translate, rout is a vector thus dont need to translate, only need to rotate
        Vector3d rout2_w = _real_R*_vrw2[i];
        Vector3d pout2_w = _real_R*_vpout2[i] + _real_t;
        Vector3d p12_w = -pout2_w + _vpout1[i];
        Vector3d plane_norm = _vrw1[i].cross(rout2_w);  //might because rout1, rout2 are small so their cross is small(not because vertical)
        Vector3d common_vertical = plane_norm;
        Vector3d normbeta = _vrw1[i].cross(common_vertical); //(N,-O,P)
        double  molecular = p12_w.transpose()*normbeta;
        double denominator = rout2_w.transpose()*normbeta;
        double kt ;
        if(molecular == 0) kt = 0;
        else kt = (double)molecular/denominator;                             //parameter k
        Vector3d estimated_p3d = pout2_w+kt*rout2_w;
        _vEstiPoint.push_back(estimated_p3d);

        //here to check the distance
        Vector3d aaa = estimated_p3d-_vpout1[i];
        Vector3d abb = aaa.cross(_vrw1[i]);
        double distancea = sqrt(abb(0)*abb(0)+abb(1)*abb(1)+abb(2)*abb(2));
        Vector3d aaa2 = estimated_p3d-pout2_w;
        Vector3d abb2 = aaa2.cross(_vrw2[i]);
        double distancea2 = sqrt(abb2(0)*abb2(0)+abb2(1)*abb2(1)+abb2(2)*abb2(2));  //should be 0
//        cout<<"distance should close to 0 "<<distancea+distancea2<<"\t";
    }
    //here we compute the mean error
    double error_sum1 = 0, error_sum2 = 0;
    for(int i = 0; i<_vEstiPoint.size(); i++){
        Vector3d point_error = _vEstiPoint[i]- _vRealPoint[i]; //enlarge the error!
        double error_point = sqrt(point_error.transpose()*point_error);
        error_sum1 += error_point;
//        cout<<"point_error[1]"<<point_error.transpose()<<endl;
        Vector3d point_error_igrefrax = _vIgnoreRefraxPoint[i] - _vRealPoint[i];
        double error_point_igrefrax = sqrt(point_error_igrefrax.transpose()*point_error_igrefrax);
        error_sum2 += error_point_igrefrax;
    }
    _error_p_mean = error_sum1/_vEstiPoint.size();
    _error_pigrefrax = error_sum2/_vEstiPoint.size();
    Vector3d t_error = _real_t - _esti_t;
    _error_t = sqrt(t_error.transpose()*t_error)/3.0;
    _error_t_norm = sqrt(t_error.transpose()*t_error)/sqrt(_real_t.transpose()*_real_t);

    cout<<"in the estimation of d[3]={da, dg}, nw = "<<endl<<fixed << setprecision(8)<<_da<<", "<<_dg<<", "<<_backProj_nw<<", mean error of estimated points = "<<_error_p_mean<<endl;
    cout<<fixed << setprecision(8)<<"error of t_ = "<<_error_t<<",normalized error = "<<_error_t_norm<<endl;
    cout<<fixed << setprecision(8)<<", mean error of ignore_refrax points = "<<_error_pigrefrax<<endl;
    cout<<"estimated 3d data has been saved to estimated3d.txt"<<endl;
}

void UnderwaterSFM::runUnderwaterSFM(String outputfilenametxt, Mat colorImg, Mat depthImg){
    getUnderwaterImg(colorImg, depthImg);
    backProjRay();
    mytriangulation();
//保存不同的3d点
//输出最终的误差文件
}