#include <iostream>
#include <curses.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <eigen3/Eigen/Dense>

class Cubic {

    private:
        Eigen::Matrix3d rot_mat;  // Actually not influence the coordinate; only work when shedding
        Eigen::Matrix3d inv_rot_mat;
        Eigen::Vector3d _view_norm_vec;  // normal vector of the viewing planet
        int light_on_surface[6];

        bool is_close(double a, double b);

        void get_light_order();

    public:

        Cubic();

        Cubic(unsigned int theta, unsigned int phi, unsigned int psi);
        
        /**
         * Find which surface the point x, y, z is on.
         * The return value would be the indexes of the four cornor of the surface.
         * @param[in] x, y, z The coordinate of the point
         * @return index of the surface, ordering x=0, y=0, z=0, x=1, y=1, z=1
        */
        int which_surface(Eigen::Vector3d vec);

        char light_of(Eigen::Vector2d coord);

        void rotate_horizonal(int theta);

        void rotate_vertical(int theta);

};