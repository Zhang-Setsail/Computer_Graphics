#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <cmath>
// add some other header files you need

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;
    // std::clog << "view" << std::endl << view << std::endl;  // check data

    return view;
}


Eigen::Matrix4f get_model_matrix(float rotation_angle, Eigen::Vector3f T, Eigen::Vector3f S, Eigen::Vector3f P0, Eigen::Vector3f P1)
{

    //Step 1: Build the Translation Matrix M_trans:
    Eigen::Matrix4f M_trans;
    M_trans << 1, 0, 0, T(0),
               0, 1, 0, T(1),
               0, 0, 1, T(2),
               0, 0, 0, 1;

    //Step 2: Build the Scale Matrix S_trans:
    Eigen::Matrix4f S_trans;
    S_trans << S(0), 0, 0, 0,
               0, S(1), 0, 0,
               0, 0, S(2), 0,
               0, 0, 0, 1;

    //Step 3: Implement Rodrigues' Rotation Formular, rotation by angle theta around axix u, then get the model matrix
	// The axis u is determined by two points, u = P1-P0: Eigen::Vector3f P0 ,Eigen::Vector3f P1  
    // Create the model matrix for rotating the triangle around a given axis. // Hint: normalize axis first

    Eigen::Matrix4f R_trans;

    Eigen::Vector3f u;
    u = P1- P0;
    u = u.normalized();
    float ux = u(0);
    float uy = u(1);
    float uz = u(2);

    float angle = rotation_angle*MY_PI/180;
    float cos_m = cos(angle);
    float sin_m = sin(angle);

    R_trans << ux*ux*(1-cos_m)+cos_m, ux*uy*(1-cos_m)-uz*sin_m, ux*uz*(1-cos_m)+uy*sin_m, 0,
               ux*uy*(1-cos_m)+uz*sin_m, uy*uy*(1-cos_m)+cos_m, uy*uz*(1-cos_m)-ux*sin_m, 0,
               ux*uz*(1-cos_m)-uy*sin_m, uy*uz*(1-cos_m)+ux*sin_m, uz*uz*(1-cos_m)+cos_m, 0,
               0, 0, 0, 1;

	//Step 4: Use Eigen's "AngleAxisf" to verify your Rotation
	//Eigen::AngleAxisf rotation_vector(radian, Vector3f(axis[0], axis[1], axis[2]));  
	//Eigen::Matrix3f rotation_matrix;
	//rotation_m = rotation_vector.toRotationMatrix();
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f translate;

    // float angle = rotation_angle*MY_PI/180;
    // float cos_m = cos(angle);
    // float sin_m = sin(angle);

    translate << cos_m, -sin_m, 0, 0,
                 sin_m, cos_m, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;
    
    // 0, 0, 1;

    model = M_trans * S_trans * R_trans * model;

    return model;
}



Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Implement this function
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    // frustum -> cubic

    // orthographic projection

    // squash all transformations

    // std::clog << "projection" << std::endl << projection << std::endl; //check
    float n = zNear;
    float f = zFar;
    float fov_radian = eye_fov*MY_PI/180;
    float wc = 1/(aspect_ratio*tan(fov_radian/2));
    float hc = 1/tan(fov_radian/2);
    float m33 = -(n+f)/(n-f);
    float m34 = 2*n*f/(n-f);

    Eigen::Matrix4f translate;
    translate << wc, 0, 0, 0, 0, hc, 0, 0, 0, 0, m33, m34, 0, 0, 1, 0;
    
    Eigen::Matrix4f minus;
    minus<< -1,0,0,0, 0,-1,0,0, 0,0,-1,0, 0,0,0,1;
    projection =  minus * translate * projection;

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "result.png";

    if (argc >= 3) {
        command_line = true;
        angle = std::stof(argv[2]); // -r by default
        if (argc == 4) {
            filename = std::string(argv[3]);
        }
        else
            return 0;
    }

    rst::rasterizer r(1024, 1024);
    rst::rasterizer r2(1024, 1024);
    // define your eye position "eye_pos" to a proper position

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};
    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};
    std::vector<Eigen::Vector3f> pos2{{1, 0, -1}, {0, 1, -1}, {-1, 0, -1}};
    std::vector<Eigen::Vector3i> ind2{{0, 1, 2}};
    // define a triangle named by "pos" and "ind"


    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto pos_id2 = r2.load_positions(pos2);
    auto ind_id2 = r2.load_indices(ind2);

    int key = 0;
    int frame_count = 0;

    // added parameters for get_projection_matrix(float eye_fov, float aspect_ratio,float zNear, float zFar)

    // Eigen::Vector3f axis(0, 0, 1);
    Eigen::Vector3f axis(1, 0, 0);
    float eye_fov, aspect_ratio, zNear, zFar;
    eye_fov = 45;
    aspect_ratio = 1;
    zNear = 0.1;
    zFar = 50;

    Eigen::Vector3f T, S, P0, P1, P2;
    T << 0, 0, 0;
    S << 1, 1, 1;
    P0 << 0, 0, 0;
    P1 << 0, 1, 0;

    // P1 << 0, 0, 1;

    int shot_time = 0;



    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, T, S, P0, P1));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(eye_fov, aspect_ratio, zNear, zFar));

        std::cout<<"end"<<std::endl;

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        r2.draw(pos_id2, ind_id2, rst::Primitive::Triangle);
        cv::Mat image(1024, 1024, CV_32FC3, r.frame_buffer().data());
        cv::Mat image2(1024, 1024, CV_32FC3, r2.frame_buffer().data());

        image = image + image2;

        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle, T, S, P0, P1));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(eye_fov, aspect_ratio, zNear, zFar));

        r2.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r2.set_model(get_model_matrix(angle, T, S, P0, P1));
        r2.set_view(get_view_matrix(eye_pos));
        r2.set_projection(get_projection_matrix(eye_fov, aspect_ratio, zNear, zFar));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        r2.draw(pos_id2, ind_id2, rst::Primitive::Triangle);

        cv::Mat image(1024, 1024, CV_32FC3, r.frame_buffer().data());
        cv::Mat image2(1024, 1024, CV_32FC3, r2.frame_buffer().data());

        image = image + image2;

        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);

        if (shot_time == 0)
        {
            /* code */
            cv::imwrite("./Result_1.png", image);
            angle = angle + 20;
            shot_time = 1;
        }
        else if (shot_time == 1)
        {
            /* code */
            cv::imwrite("./Result_2.png", image);
            angle = angle + 20;
            shot_time = 2;
        }
        else if (shot_time == 2)
        {
            /* code */
            cv::imwrite("./Result_3.png", image);
            angle = 0;
            shot_time = 3;
        }

        // image2.convertTo(image2, CV_8UC3, 1.0f);
        // cv::imshow("image", image2);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
        std::clog << "angle: " << angle << std::endl;
    

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
