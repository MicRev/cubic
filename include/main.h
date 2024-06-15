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

        /**
         * Find which surface the point x, y, z is on.
         * The return value would be the indexes of the four cornor of the surface.
         * @param[in] x, y, z The coordinate of the point
         * @return index of the surface, ordering x=0, y=0, z=0, x=1, y=1, z=1
        */
        int which_surface(Eigen::Vector3d vec);

        /**
         * Find a point on surface first met tracing from certain point through certain direction
         * @param[in] point the point to begin tracing
         * @param[in] direction the direction the tracing obey
         * @return the point on surface. It is ensured that the returned point is the first point met on surface. When no point found, the function return (0, 0, 0)
        */
        Eigen::Vector3d point_through_vec_to_surface(Eigen::Vector3d point, Eigen::Vector3d direction);

    public:

        Cubic();

        Cubic(unsigned int theta, unsigned int phi, unsigned int psi);

        char light_of(Eigen::Vector2d coord);

        char shadow_of(Eigen::Vector2d coord);

        void rotate_horizonal(int theta);

        void rotate_vertical(int theta);

        void rotate_view_vertical(int theta);

        void rotate_view_horizonal(int theta);

        void reset();

};