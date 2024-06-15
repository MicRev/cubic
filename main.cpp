#include "main.h"

// #define DEBUG

#ifdef DEBUG
#include<fstream>
std::fstream log_file;
bool printed = false;
#endif

const char lightness[6] = {'.', '-', ':', '=', '#', '@'};
const char shadow = '\\';
const double ground_z = -sqrt(3.) / 2;
const double length = 1.;

const int D = 100;
const double CERTAINTY = 1e-6;
Eigen::Vector3d light(1e-5, -0.5, -0.86602540378443864676372317075294); // (0, -1/2, -\sqrt{3}/2)


bool Cubic::is_close(double a, double b) {
    return fabs(a - b) < CERTAINTY;
}

void Cubic::get_light_order() {
    
    Eigen::Vector3d _light = inv_rot_mat * light;
    double dot_val[6];
    for (int i = 0; i < 6; ++i) {
        int idx = i % 3;
        int val = i > 2 ? 1 : -1;
        Eigen::Vector3d _norm(0, 0, 0);
        _norm(idx) = val;
        dot_val[i] = _light.dot(_norm);
    }

    for (int i = 0; i < 6; ++i) {
        light_on_surface[i] = 0;
        for (int j = 0; j < 6; ++j) {
            if (dot_val[i] > dot_val[j]) light_on_surface[i]++;
        }
    }
}

Eigen::Vector3d Cubic::point_through_vec_to_surface(Eigen::Vector3d point, Eigen::Vector3d direction) {
    Eigen::Vector3d unrot_vec       = inv_rot_mat * point;
    Eigen::Vector3d unrot_view_norm = inv_rot_mat * direction;

    for (int i = 0; i < 6; ++i) {  // find meeting point of eyeshot through this point with 6 surfaces

        int idx = i % 3;  // x, y, or z
        double val = double(i / 3) - 0.5; // -0.5 for 0, 1, 2 and 0.5 for 3, 4, 5

        if (val * unrot_view_norm(idx) > 0) continue;  // make sure the point is on the surface can be seen

        Eigen::Vector3d p;
        p(idx) = val;
        bool is_valid_point = true;
        
        for (int j = 0; j < 3; ++j) {
            if (j == idx) continue;
            p(j) = unrot_vec(j) + (val - unrot_vec(idx)) * unrot_view_norm(j) / unrot_view_norm(idx);
            if (p(j) < -0.5 || p(j) > 0.5) {
                is_valid_point = false; // not the point
            }
        }

        if (!is_valid_point) continue;

        return p;
        break;
    }
    Eigen::Vector3d p(0, 0, 0);
    return p;
}

Cubic::Cubic() {
    rot_mat     << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    inv_rot_mat << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    get_light_order();
    _view_norm_vec << -0.8660254037844*0.5, -0.8660254037844*-0.8660254037844, -0.5;
}

Cubic::Cubic(unsigned int theta, unsigned int phi, unsigned int psi) : Cubic() {

    theta = theta % 360;
    phi = phi % 360;
    psi = psi % 180;

    double d_theta, d_phi, d_psi;
    d_theta = theta / 180. * M_PI;
    d_phi = phi / 180. * M_PI;

    Eigen::Matrix3d theta_rot_mat, phi_rot_mat, psi_rot_mat;
    theta_rot_mat << cos(d_theta), sin(d_theta), 0., 
                    -sin(d_theta), cos(d_theta), 0., 
                     0.            , 0.             , 1.;
    phi_rot_mat   << 1.,  0.          , 0.           ,
                     0.,  cos(d_phi), sin(d_phi),
                     0., -sin(d_phi), cos(d_phi);

    #ifdef DEBUG
        // log_file << "Theta" << theta << "\n" << theta_rot_mat << std::endl;
        // log_file << "Phi"   << phi   << "\n" << phi_rot_mat   << std::endl;
    #endif

    rot_mat = phi_rot_mat * theta_rot_mat;
    inv_rot_mat = rot_mat.inverse();
    get_light_order();
    _view_norm_vec << -0.8660254037844*0.5, -0.8660254037844*-0.8660254037844, -0.5;
}

int Cubic::which_surface(Eigen::Vector3d vec) {
    Eigen::Vector3d unrot_vec = vec;
    if (is_close(unrot_vec(0), -0.5)) return 0;
    else if (is_close(unrot_vec(1), -0.5)) return 1;
    else if (is_close(unrot_vec(2), -0.5)) return 2;
    else if (is_close(unrot_vec(0), 0.5)) return 3;
    else if (is_close(unrot_vec(1), 0.5)) return 4;
    else if (is_close(unrot_vec(2), 0.5)) return 5;
    else {
        std::cout << "This point is not on any surface!" << std::endl;
        throw "This point is not on any surface!";
    }
}

char Cubic::light_of(Eigen::Vector2d coord) {

    #ifdef DEBUG
        // log_file << "coord\n" << coord << std::endl;
    #endif

    Eigen::Vector3d vec;

    #ifdef DEBUG
        if (!printed) {
            log_file << "rot_mat\n" << rot_mat << std::endl;
            log_file << "inv_rot_mat\n" << inv_rot_mat << std::endl;
            log_file << "unitary test\n" << rot_mat * rot_mat.transpose() << std::endl;
            printed = true;
        }
    #endif

    double nx = _view_norm_vec(0);
    double ny = _view_norm_vec(1);
    double nz = _view_norm_vec(2);

    Eigen::Matrix<double, 3, 2> base_mat;
    double norm_coord_x = sqrt(ny*ny + nx*nx);
    double norm_coord_y = norm_coord_x;
    base_mat <<  ny / norm_coord_x, -nx*nz        / norm_coord_y, 
                -nx / norm_coord_x, -ny*nz        / norm_coord_y, 
                 0.               , (nx*nx+ny*ny) / norm_coord_y; 
    // col 1: \vec{b}_1 = \vec{n} \times \hat{k}; col 2: \vec{b}_2 = \vec{b}_1 \times \vec{n}, then normalize

    #ifdef DEBUG
        // log_file << "Base matrix\n" << base_mat << std::endl; 
    #endif

    vec = base_mat * coord + -D * _view_norm_vec;

    #ifdef DEBUG
        log_file << "Vec in origin rotated space\n" << vec << std::endl; 
    #endif

    Eigen::Vector3d point_on_surface = point_through_vec_to_surface(vec, _view_norm_vec);

    if (point_on_surface(0) == 0 && point_on_surface(1) == 0 && point_on_surface(2) == 0) {
        return ' ';
    }

    int surface = which_surface(point_on_surface);

    #ifdef DEBUG
        log_file << "coord\n" << coord << std::endl;
        log_file << "Found at point\n" << point_on_surface << std::endl;
        log_file << "On surface " << surface << std::endl;
    #endif
    get_light_order();
    return lightness[light_on_surface[surface]];
}

char Cubic::shadow_of(Eigen::Vector2d coord) {
    Eigen::Vector3d vec;
    double nx = _view_norm_vec(0);
    double ny = _view_norm_vec(1);
    double nz = _view_norm_vec(2);

    Eigen::Matrix<double, 3, 2> base_mat;
    double norm_coord_x = sqrt(ny*ny + nx*nx);
    double norm_coord_y = norm_coord_x;
    base_mat <<  ny / norm_coord_x, -nx*nz        / norm_coord_y, 
                -nx / norm_coord_x, -ny*nz        / norm_coord_y, 
                 0.               , (nx*nx+ny*ny) / norm_coord_y;
    
    vec = base_mat * coord + -D * _view_norm_vec;

    Eigen::Vector3d ground_point(0, 0, ground_z);
    ground_point(0) = nx / nz * (ground_z - vec(2)) + vec(0);
    ground_point(1) = ny / nz * (ground_z - vec(2)) + vec(1);

    Eigen::Vector3d point_on_surface = point_through_vec_to_surface(ground_point, -light);
    if (point_on_surface(0) == 0 && point_on_surface(1) == 0 && point_on_surface(2) == 0) {
        return ' ';
    } else {
        return shadow;
    }
}

void Cubic::rotate_horizonal(int theta) {
    double d_theta = theta / 180. * M_PI;

    Eigen::Matrix3d cur_rot_mat;
    cur_rot_mat <<  cos(d_theta), sin(d_theta), 0,
                    -sin(d_theta), cos(d_theta), 0,
                    0             , 0              , 1;

    double nx = _view_norm_vec(0);
    double ny = _view_norm_vec(1);
    double nz = _view_norm_vec(2);
    double norm_coord_view = sqrt(nx*nx+ny*ny);

    Eigen::Matrix3d base_transform_to_view_mat;
    base_transform_to_view_mat << ny/norm_coord_view, -nx/norm_coord_view, 0,
                                  nx/norm_coord_view,  ny/norm_coord_view, 0,
                                  0                 ,  0                 , 1;  // unitary

    rot_mat = base_transform_to_view_mat.transpose() * cur_rot_mat * base_transform_to_view_mat * rot_mat;
    inv_rot_mat = rot_mat.inverse();
}

void Cubic::rotate_vertical(int theta) {
    double d_theta = theta / 180. * M_PI;

    double nx = _view_norm_vec(0);
    double ny = _view_norm_vec(1);
    double nz = _view_norm_vec(2);
    double norm_coord_view = sqrt(nx*nx+ny*ny);

    Eigen::Matrix3d base_transform_to_view_mat;
    base_transform_to_view_mat << ny/norm_coord_view, -nx/norm_coord_view, 0,
                                  nx/norm_coord_view,  ny/norm_coord_view, 0,
                                  0                 ,  0                 , 1;  // unitary
    Eigen::Matrix3d cur_rot_mat;
    cur_rot_mat << 1, 0             ,  0              ,
                   0, cos(d_theta), -sin(d_theta),
                   0, sin(d_theta),  cos(d_theta);

    rot_mat = base_transform_to_view_mat.transpose() * cur_rot_mat * base_transform_to_view_mat * rot_mat;
    inv_rot_mat = rot_mat.inverse();
}

void Cubic::rotate_view_horizonal(int theta) {
    double d_theta = theta / 180. * M_PI;
    Eigen::Matrix3d cur_rot_mat;
    cur_rot_mat <<  cos(d_theta), sin(d_theta), 0,
                    -sin(d_theta), cos(d_theta), 0,
                    0             , 0              , 1;
    _view_norm_vec = cur_rot_mat.transpose() * _view_norm_vec;
}

void Cubic::rotate_view_vertical(int theta) { // !BUG here
    double d_theta = theta / 180. * M_PI;
    Eigen::Matrix3d cur_rot_mat;
    cur_rot_mat << 1, 0             ,  0              ,
                   0, cos(d_theta), -sin(d_theta),
                   0, sin(d_theta),  cos(d_theta);
    
    double nx = _view_norm_vec(0);
    double ny = _view_norm_vec(1);
    double nz = _view_norm_vec(2);
    double norm_coord_view = sqrt(nx*nx+ny*ny);

    Eigen::Matrix3d base_transform_to_view_mat;
    base_transform_to_view_mat << ny/norm_coord_view, -nx/norm_coord_view, 0,
                                  nx/norm_coord_view,  ny/norm_coord_view, 0,
                                  0                 ,  0                 , 1;  // unitary

    _view_norm_vec = base_transform_to_view_mat * cur_rot_mat.transpose() * base_transform_to_view_mat.transpose() * _view_norm_vec;

}

void Cubic::reset() {
    rot_mat << 1., 0., 0., 0., 1., 0., 0., 0., 1.;
    inv_rot_mat = rot_mat;
    _view_norm_vec << -0.8660254037844*0.5, -0.8660254037844*-0.8660254037844, -0.5;
}

int main() {

    #ifdef DEBUG
    log_file.open("../log.log", std::ios::out);
    #endif

    initscr();
    noecho();  // do not show input chars
    curs_set(0);  // hide cursor
    // nodelay(stdscr, true);
    // leaveok(stdscr, true);
    clear();
    refresh();

    winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);  // get current windowsize 

    LINES = w.ws_row;
    COLS  = w.ws_col;

    #define min(a, b) ((a) < (b) ? (a) : (b))
    int size = min(LINES, COLS/2);
    #undef min

    #ifdef DEBUG
    log_file << "Lines\tCols\tsize\ti begin\ti end\tj begin\tj end\n" 
             << LINES << "\t" << COLS << "\t" << size <<"\t"
             << i_begin << "\t" << i_end << "\t" << j_begin << "\t" << j_end << std::endl;
    #endif

    Cubic cubic(0, 0, 0);
    while (1) {

        #ifdef DEBUG
            printed = false;
        #endif

        for (int i = 0; i < COLS; ++i) {
            for (int j = 0; j < LINES; ++j) {
                double x = (i - COLS/2) * 2. / 2.5 / size + 1e-6;
                double y = (LINES/2 - j) * 2. / size + 1e-6;  // 1e-6 for numerical stability
                #ifdef DEBUG
                    // log_file << "i " << i << " j " << j << " x " << x << " y " << y << std::endl;
                #endif
                Eigen::Vector2d coord(x, y);
                char cur_char;
                cur_char = cubic.shadow_of(coord);
                if (cur_char != ' ') mvaddch(j, i, cur_char);
                cur_char = cubic.light_of(coord);
                if (cur_char != ' ') mvaddch(j, i, cur_char);
            }
        }
        refresh();
        int k = getch();

        #ifdef DEBUG
            log_file << "Pressed " << k << std::endl;   
        #endif

        bool to_break = false;
        switch (k) {
            case 'w':
                cubic.rotate_vertical(-10);
                clear();
                break;
            case 's':
                cubic.rotate_vertical( 10);
                clear();
                break;
            case 'a':
                cubic.rotate_horizonal( 10);
                clear();
                break;
            case 'd':
                cubic.rotate_horizonal(-10);
                clear();
                break;
            // case 'i':
            //     cubic.rotate_view_vertical(-10);
            //     clear();
            //     break;
            // case 'k':
            //     cubic.rotate_view_vertical( 10);
            //     clear();
            //     break;
            case 'j':
                cubic.rotate_view_horizonal( -10);
                clear();
                break;
            case 'l':
                cubic.rotate_view_horizonal(10);
                clear();
                break;
            case 'r':
                cubic.reset();
                clear();
                break;
            case 'q':
                to_break = true;
                break;
            default:
                break;
        }
        if (to_break) break;
    }

    #ifdef DEBUG
        log_file.close();
    #endif
    endwin();
}