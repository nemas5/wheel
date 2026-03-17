#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <functional>
 
using namespace std;

// Для расширения либо создавать набор экземпляров для разных поверхностей,
// либо устанавливать зависимость от времени для коэффициентов.
class Ground {
    private:
        double c_gr;  // Коэффициент деформации грунта
        double nu;  // Показатель плотности грунта
        double phi_x, phi_y;  // Продольный и поперечный коэффициенты анизотропного трения
        double a1, a2;  // Волшебные коэффициенты грунта из работы Кацыгина (a1 связан с рыхлостью грунта)
    public:
        Ground(double a, double b, double c, double d, double e, double f) 
            : c_gr(a), nu(b), phi_x(c), phi_y(d), a1(e), a2(f) {};
        double z_gr_nd(double x, double y) {return 0.;};
        double get_c_gr() {return c_gr;};
        double get_nu() {return nu;};
        double get_a1() {return a1;};
        double get_a2() {return a2;};
        double get_phi_x() {return phi_x;};
        double get_phi_y() {return phi_y;};
};
 
class Wheel {
    private:
        struct WheelStepState {
            double x_0, y_0, z_0;  // Координаты центра колеса
            double v_x, v_y, v_z;  // Скорости центра колеса
            double w;
            WheelStepState(double x_, double y_, double z_, double vx, double vy, double vz, double ww)
                : x_0(x_), y_0(y_), z_0(z_), v_x(vx), v_y(vy), v_z(vz), w(ww) {};
        };
        double hg_max_prev;  // Хранится глубина колеи (наибольшая) только за предыдущий шаг (необходимо для вычисления L1)
        vector<double> hg_prev;  // А тут всё распределение глубины колеи вдоль ширины колеса на предыдущем шагу
        vector<vector<double>> z_prev;  // Для вычисления скорости деформации профиля точки
        vector<WheelStepState> states;  // Изменяющееся состояние колеса по ходу движения
        double r, b;  // Радиус и ширина колеса
        double j;  // Момент инерции колеса
        double m;  // Масса колеса
        double steps_x, steps_y;  // Кол-во точек на пятне контакта
        double dx, dy;  // Расстояние между точками на пятне контакта
        double r_elastic, r_damping;  // Упругие и демпфирующие свойства колеса
        Ground gr;  // Поверхность движения
    public:
        Wheel(double wheel_radius, double wheel_width, double x, double y, 
            Ground ground, double wheel_j, double m, double delta_x, double delta_y);
        void modeling(double mk, double f_side, double t_start, double t_end, double dt);
        pair<double, double> calc_c(double vx, double vy, double w);  // Вычисление моментального центра скольжения колеса
        double calc_alfa(double x);  // Посчитать угол альфа для точки контакта
        double calc_hg(double q);  // Глубина колеи
        double calc_xg(double hg);  // Длина колеи ?)
        double calc_zgri(double x, double y, double xg, double hg);  // Деформированный профиль грунта в точке
        double calc_l1();  // Длина пятна контакта после центра колеса
        double calc_l2(const vector<double>& z, const vector<double>& x);  // Длина пятна контакта до центра колеса
        double calc_f(double l1, double l2);  // Площадь пятна контакта
        double calc_qmax(double f);  // Максимальное нормальное давление в пятне контакта
        double calc_q(double qmax, double y);  // Нормальное давление в точке пятна контакта
        double calc_dri(double z, double cos_alfa);  // Деформация шины в радиальном направлении в точке пятна контакта
        double calc_vri(double dr, double alfa, double w, double vx);  // Скорость токи пятна контакта в радиальном направлении
        double calc_dR(double dr, double vr, double cos_alfa, double z, double z_prev, double dt);  // Реакция в точке пятна контакта (в радиальном направлении)
        double calc_phi_x(double vx, double w, double dri);  // Вычисление коэффициента трения в точке пятна контакта для продольного направления
        double calc_phi_y(double slip);  // Для поперечного (одно значение для всех точек)
        vector<double> calc_contact_reactions(double vx_add, double vy_add, double w_add, double dt);
        vector<double> solve_integral(
            double xc, double yc, const vector<double>& x, 
            const vector<double>& y, const vector<vector<double>>& phi_x, double phi_y,
            const vector<vector<double>>& dR,
            vector<vector<double>> dr  // TODO: это можно сделать ссылкой?
        );  // Интегрирование по пятну контакта
        void out();
};

void Wheel::out()
{
    ofstream file("file.csv");
    if (!file.is_open()) return;

    file << "x_0,y_0,v_x,v_y,w\n";

    for (const auto& s : states) {
        file << s.x_0 << ","
             << s.y_0 << ","
             // << s.z_0 << ","
             << s.v_x << ","
             << s.v_y << ","
             // << s.v_z << ","
             << s.w   << "\n";
    }

    file.close();
}

// Составная формула трапеций для интегрирования по площади контакта
vector<double> Wheel::solve_integral(  // TODO: повторения кода вырезать
    double xc, double yc, const vector<double>& x, const vector<double>& y,
    const vector<vector<double>>& phi_x, double phi_y, const vector<vector<double>>& dR,
    vector<vector<double>> dr
)
{
    double p_x = 0.;
    double p_y = 0.;
    double m_x_cos = 0., m_x_sin = 0.;
    double p_x_cos = 0., p_x_sin = 0.;
    double p_y_cos = 0.;
    for (int i = 1; i < steps_x - 1; i++) {
        for (int j = 1; j < steps_y - 1; j++) {
            dr[i][j] = r - dr[i][j];
        }
    }

    // cout << xc << " " << yc << endl;
    // for (int i = 0; i < x.size(); i++) {
    //     cout << x[i] << " ";
    // }
    // cout << endl;
    // for (int i = 0; i < y.size(); i++) {
    //     cout << y[i] << " ";
    // }
    // for (int i = 0; i < phi_x.size(); i++) {
    //     for (int j = 0; j < phi_x[i].size(); j++) {
    //         cout << phi_x[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // for (int i = 0; i < dR.size(); i++) {
    //     for (int j = 0; j < dR[i].size(); j++) {
    //         cout << dR[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // for (int i = 0; i < dr.size(); i++) {
    //     for (int j = 0; j < dr[i].size(); j++) {
    //         cout << dr[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;


    // Внутренние точки
    for (int i = 1; i < steps_x - 2; i++) {
        for (int j = 1; j < steps_y - 2; j++) {
            if (abs(xc - x[i]) < dx && abs(yc - y[j]) < dy) {
                m_x_cos += (phi_x[i][j] * dR[i][j] * cos(calc_alfa(x[i])) * dr[i][j]);
                m_x_sin += (phi_x[i][j] * dR[i][j] * sin(calc_alfa(x[i])) * dr[i][j]);
                p_x_cos += (phi_x[i][j] * dR[i][j] * cos(calc_alfa(x[i])));
                p_x_sin += (phi_x[i][j] * dR[i][j] * sin(calc_alfa(x[i])));
                p_y_cos += (-phi_y * dR[i][j] * cos(calc_alfa(x[i])));
            }
            else {
                m_x_cos += (phi_x[i][j] * (yc - y[j]) * dR[i][j] * cos(calc_alfa(x[i])) * dr[i][j]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                m_x_sin += (phi_x[i][j] * (yc - y[j]) * dR[i][j] * sin(calc_alfa(x[i])) * dr[i][j]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_x_cos += (phi_x[i][j] * (yc - y[j]) * dR[i][j] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_x_sin += (phi_x[i][j] * (yc - y[j]) * dR[i][j] * sin(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_y_cos += (-phi_y * (xc - x[i]) * dR[i][j] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
            }
        }
    }
    
    // Рёбра без углов по x
    for (int i = 1; i < steps_x - 2; i++) {
        if (abs(xc - x[i]) < dx && abs(yc - y[0]) < dy) {
            m_x_cos += 0.5 * (phi_x[i][0] * dR[i][0] * cos(calc_alfa(x[i])) * dr[i][0]);
            m_x_sin += 0.5 * (phi_x[i][0] * dR[i][0] * sin(calc_alfa(x[i])) * dr[i][0]);
            p_x_cos += 0.5 * (phi_x[i][0] * dR[i][0] * cos(calc_alfa(x[i])));
            p_x_sin += 0.5 * (phi_x[i][0] * dR[i][0] * sin(calc_alfa(x[i])));
            p_y_cos += 0.5 * (-phi_y * dR[i][0] * cos(calc_alfa(x[i])));
        }
        else {
            m_x_cos += 0.5 * (phi_x[i][0] * (yc - y[0]) * dR[i][0] * cos(calc_alfa(x[i])) * dr[i][0]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[0], 2));
            m_x_sin += 0.5 * (phi_x[i][0] * (yc - y[0]) * dR[i][0] * sin(calc_alfa(x[i])) * dr[i][0]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[0], 2));
            p_x_cos += 0.5 * (phi_x[i][0] * (yc - y[0]) * dR[i][0] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[0], 2));
            p_x_sin += 0.5 * (phi_x[i][0] * (yc - y[0]) * dR[i][0] * sin(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[0], 2));
            p_y_cos += 0.5 * (-phi_y * (xc - x[i]) * dR[i][0] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[0], 2));
        }

        if (abs(xc - x[i]) < dx && abs(yc - y[steps_y - 2]) < dy) {
            m_x_cos += 0.5 * (phi_x[i][steps_y - 2] * dR[i][steps_y - 2] * cos(calc_alfa(x[i])) * dr[i][steps_y - 2]);
            m_x_sin += 0.5 * (phi_x[i][steps_y - 2] * dR[i][steps_y - 2] * sin(calc_alfa(x[i])) * dr[i][steps_y - 2]);
            p_x_cos += 0.5 * (phi_x[i][steps_y - 2] * dR[i][steps_y - 2] * cos(calc_alfa(x[i])));
            p_x_sin += 0.5 * (phi_x[i][steps_y - 2] * dR[i][steps_y - 2] * sin(calc_alfa(x[i])));
            p_y_cos += 0.5 * (-phi_y * dR[i][steps_y - 2] * cos(calc_alfa(x[i])));
        }
        else {
            m_x_cos += 0.5 * (phi_x[i][steps_y - 2] * (yc - y[steps_y - 2]) * dR[i][steps_y - 2] * cos(calc_alfa(x[i])) * dr[i][steps_y - 2]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[steps_y - 2], 2));
            m_x_sin += 0.5 * (phi_x[i][steps_y - 2] * (yc - y[steps_y - 2]) * dR[i][steps_y - 2] * sin(calc_alfa(x[i])) * dr[i][steps_y - 2]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[steps_y - 2], 2));
            p_x_cos += 0.5 * (phi_x[i][steps_y - 2] * (yc - y[steps_y - 2]) * dR[i][steps_y - 2] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[steps_y - 2], 2));
            p_x_sin += 0.5 * (phi_x[i][steps_y - 2] * (yc - y[steps_y - 2]) * dR[i][steps_y - 2] * sin(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[steps_y - 2], 2));
            p_y_cos += 0.5 * (-phi_y * (xc - x[i]) * dR[i][steps_y - 2] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[steps_y - 2], 2));

        }
    }

    // Рёбра без углов по y
    for (int j = 1; j < steps_y - 2; j++) {
        
        if (abs(xc - x[0]) < dx && abs(yc - y[j]) < dy) {
            m_x_cos += 0.5 * (phi_x[0][j] * dR[0][j] * cos(calc_alfa(x[0])) * dr[0][j]);
            m_x_sin += 0.5 * (phi_x[0][j] * dR[0][j] * sin(calc_alfa(x[0])) * dr[0][j]);
            p_x_cos += 0.5 * (phi_x[0][j] * dR[0][j] * cos(calc_alfa(x[0])));
            p_x_sin += 0.5 * (phi_x[0][j] * dR[0][j] * sin(calc_alfa(x[0])));
            p_y_cos += 0.5 * (-phi_y * dR[0][j] * cos(calc_alfa(x[0])));
        }
        else {
            m_x_cos += 0.5 * (phi_x[0][j] * (yc - y[j]) * dR[0][j] * cos(calc_alfa(x[0])) * dr[0][j]) / sqrt(pow(xc - x[0], 2) + pow(yc - y[j], 2));
            m_x_sin += 0.5 * (phi_x[0][j] * (yc - y[j]) * dR[0][j] * sin(calc_alfa(x[0])) * dr[0][j]) / sqrt(pow(xc - x[0], 2) + pow(yc - y[j], 2));
            p_x_cos += 0.5 * (phi_x[0][j] * (yc - y[j]) * dR[0][j] * cos(calc_alfa(x[0]))) / sqrt(pow(xc - x[0], 2) + pow(yc - y[j], 2));
            p_x_sin += 0.5 * (phi_x[0][j] * (yc - y[j]) * dR[0][j] * sin(calc_alfa(x[0]))) / sqrt(pow(xc - x[0], 2) + pow(yc - y[j], 2));
            p_y_cos += 0.5 * (-phi_y * (xc - x[0]) * dR[0][j] * cos(calc_alfa(x[0]))) / sqrt(pow(xc - x[0], 2) + pow(yc - y[j], 2));
        }


        if (abs(xc - x[steps_x - 2]) < dx && abs(yc - y[j]) < dy) {
            m_x_cos += 0.5 * (phi_x[steps_x - 2][j] * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2])) * dr[steps_x - 2][j]);
            m_x_sin += 0.5 * (phi_x[steps_x - 2][j] * dR[steps_x - 2][j] * sin(calc_alfa(x[steps_x - 2])) * dr[steps_x - 2][j]);
            p_x_cos += 0.5 * (phi_x[steps_x - 2][j] * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2])));
            p_x_sin += 0.5 * (phi_x[steps_x - 2][j] * dR[steps_x - 2][j] * sin(calc_alfa(x[steps_x - 2])));
            p_y_cos += 0.5 * (-phi_y * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2])));
        }
        else {
            m_x_cos += 0.5 * (phi_x[steps_x - 2][j] * (yc - y[j]) * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2])) * dr[steps_x - 2][j]) / sqrt(pow(xc - x[steps_x - 2], 2) + pow(yc - y[j], 2));
            m_x_sin += 0.5 * (phi_x[steps_x - 2][j] * (yc - y[j]) * dR[steps_x - 2][j] * sin(calc_alfa(x[steps_x - 2])) * dr[steps_x - 2][j]) / sqrt(pow(xc - x[steps_x - 2], 2) + pow(yc - y[j], 2));
            p_x_cos += 0.5 * (phi_x[steps_x - 2][j] * (yc - y[j]) * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2]))) / sqrt(pow(xc - x[steps_x - 2], 2) + pow(yc - y[j], 2));
            p_x_sin += 0.5 * (phi_x[steps_x - 2][j] * (yc - y[j]) * dR[steps_x - 2][j] * sin(calc_alfa(x[steps_x - 2]))) / sqrt(pow(xc - x[steps_x - 2], 2) + pow(yc - y[j], 2));
            p_y_cos += 0.5 * (-phi_y * (xc - x[steps_x - 2]) * dR[steps_x - 2][j] * cos(calc_alfa(x[steps_x - 2]))) / sqrt(pow(xc - x[steps_x - 2], 2) + pow(yc - y[j], 2));
        }

    }

    // Углы
    for (int i = 0; i < steps_x - 1; i += steps_x - 2) {
        for (int j = 0; j < steps_y - 1; j += steps_y - 2) {

            if (abs(xc - x[i]) < dx && abs(yc - y[j]) < dy) {
                m_x_cos += 0.25 * (phi_x[i][j] * dR[i][j] * cos(calc_alfa(x[i])) * dr[i][j]);
                m_x_sin += 0.25 * (phi_x[i][j] * dR[i][j] * sin(calc_alfa(x[i])) * dr[i][j]);
                p_x_cos += 0.25 * (phi_x[i][j] * dR[i][j] * cos(calc_alfa(x[i])));
                p_x_sin += 0.25 * (phi_x[i][j] * dR[i][j] * sin(calc_alfa(x[i])));
                p_y_cos += 0.25 * (-phi_y * dR[i][j] * cos(calc_alfa(x[i])));
            }
            else {
                m_x_cos += 0.25 * (phi_x[i][j] * (yc - y[j]) * dR[i][j] * cos(calc_alfa(x[i])) * dr[i][j]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                m_x_sin += 0.25 * (phi_x[i][j] * (yc - y[j]) * dR[i][j] * sin(calc_alfa(x[i])) * dr[i][j]) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_x_cos += 0.25 * (phi_x[i][j] * (yc - y[j]) * dR[i][j] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_x_sin += 0.25 * (phi_x[i][j] * (yc - y[j]) * dR[i][j] * sin(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
                p_y_cos += 0.25 * (-phi_y * (xc - x[i]) * dR[i][j] * cos(calc_alfa(x[i]))) / sqrt(pow(xc - x[i], 2) + pow(yc - y[j], 2));
            }
        }
    }
    // cout << m_x_cos << " " << m_x_sin << " " << p_x_cos << " " << p_x_sin << " " << p_y_cos << endl;
    return {m_x_cos * dx * dy, m_x_sin * dx * dy, p_x_cos * dx * dy, p_x_sin * dx * dy, p_y_cos * dx * dy};
}

pair<double, double> Wheel::calc_c(double vx, double vy, double w) {
    if (vx == 0 && vy == 0 && w == 0) {
        return {0., 0.};
    }
    // В системе координат колеса
    if (w < 1e-8) {
        return {0., 0.};
    }
    double x_c = -vy / w;
    double y_c = vx / w;
    return {x_c, y_c};
}

double Wheel::calc_xg(double hg) {
    return sqrt(pow(r, 2) - pow((states.back().z_0 - hg), 2));
}

double Wheel::calc_l1() {
    // cout << "hg:" << hg_max_prev << " " << calc_xg(hg_max_prev) << endl;
    double xg = calc_xg(hg_max_prev);  // hg с предыдущего шага
    if (xg < 10e-8)
        return 0;
    double a = hg_max_prev / pow(xg, 2);
    if (a < 1e-15)
        return 0;
    return (1 / (2 * a)) * (a * xg * sqrt(a * pow(xg, 2) + 1) + sqrt(a) * pow(sinh(xg * sqrt(a)), -1));
}

double Wheel::calc_l2(const vector<double>& z, const vector<double>& x) {
    double sm = 0.;
    for (int i = 0; i < z.size() - 1; i++) {
        sm += sqrt(pow(z[i + 1] - z[i], 2) + pow(x[i + 1] - x[i], 2));
    }
    return sm;
}

double Wheel::calc_alfa(double x) {
    return asin(x / r);
}

double Wheel::calc_f(double l1, double l2) {
    return (l1 + l2) * b;
}

double Wheel::calc_qmax(double f) {
    return (4 * m * 9.8) / (M_PI * f); 
}

double Wheel::calc_q(double qmax, double y) {
    // cout << "inside q: " << qmax << " " << y << " " << b << sqrt(1 - pow((2 * y / b - 1), 2)) << endl;
    return qmax * sqrt(1 - pow((2 * y / b - 1), 2));
}

double Wheel::calc_hg(double q) {
    return pow(q / (gr.get_c_gr() * 1e6), 1 / gr.get_nu());
}

double Wheel::calc_zgri(double x, double y, double xg, double hg) {
    double z = gr.z_gr_nd(x, y);
    if (x >= xg) {
        return z;
    }
    else if (x <= 0) {
        return z - hg;
    }
    else {
        return z - (1 - (hg / pow(xg, 2)) * pow(x, 2));
    }
}

double Wheel::calc_dri(double z, double cos_alfa) {
    double z_r = states.back().z_0 - r * cos_alfa;
    if (z > z_r)
        return (z - z_r) * cos_alfa;
    return 0.;
}

double Wheel::calc_vri(double dr, double alfa, double w, double vx) {
    double vz = states.back().v_z;
    vx = w * (r - dr) * cos(alfa) + vx;
    vz = w * (r - dr) * sin(alfa) + vz;
    return vx * sin(alfa) + vz * cos(alfa);
}

double Wheel::calc_dR(double dr, double vr, double cos_alfa, double z, double z_prev, double dt) {
    // cout << dr << " " << vr << " " << cos_alfa << " " << z << z_prev << " " << dt << " ";
    double z_diff = (z - z_prev) / dt;
    return dr - vr + cos_alfa * z_diff;
}

Wheel::Wheel(double wheel_radius, double wheel_width, double x, double y, Ground ground, 
    double wheel_j, double mass, double delta_x, double delta_y)
    : r(wheel_radius), b(wheel_width), gr(ground), j(wheel_j), m(mass), dx(delta_x), dy(delta_y)
{
    steps_x = (r * 2.) / dx + 1;
    steps_y = b / dy + 1;
    hg_max_prev = 0.;
    hg_prev.resize(steps_y);
    z_prev.assign(steps_x - 1, std::vector<double>(steps_y - 1, 0.0));
    for (int i = 0; i < steps_x - 1; i++) {
        for (int j = 0; j < steps_y - 1; j++) {
            z_prev[i][j] = gr.z_gr_nd(i * dx, j * dy);
        }
    }
    states.emplace_back(x, y, gr.z_gr_nd(0, 0) + r, 0., 0., 0., 0.);
    //     WheelStepState(double x_, double y_, double z_, double vx, double vy, double vz, double ww)
    //         : x_0(x_), y_0(y_), z_0(z_), v_x(vx), v_y(vy), v_z(vz), w(ww) {};
}

double Wheel::calc_phi_y(double slip) {
    return gr.get_phi_y() * (1 + (gr.get_a1() / cosh(slip / gr.get_a2()))) * tanh(slip / gr.get_a2());
}

double Wheel::calc_phi_x(double vx, double w, double dri) {
    if (w < 1e-8)
        return gr.get_phi_x();
    double slip = (vx - w * (r - dri)) / (w * (r - dri));
    return gr.get_phi_x() * (1 + (gr.get_a1() / cosh(slip / gr.get_a2()))) * tanh(slip / gr.get_a2());
}

// hg снаружи пересчитываем
vector<double> Wheel::calc_contact_reactions(double vx_add, double vy_add, double w_add, double dt) {
    WheelStepState last_step_state = states.back();
    double vx = last_step_state.v_x + vx_add;
    double w = last_step_state.w + w_add;
    double vy = last_step_state.v_y + vy_add;
    pair<double, double> c = calc_c(vx, vy, w);
    double cur_alfa, cur_z;
    vector<double> x, y;
    vector<vector<double>> dr, dR;
    vector<vector<double>> phi_x;
    double phi_y;
    if (abs(vx) < 1e-8)
        phi_y = gr.get_phi_y();
    else
        phi_y = calc_phi_y(vy / abs(vx));
    int i, j;
    for (j = 0; j < steps_y - 1; j++)
        y.push_back((-b / 2.) + j * dy);
    for (i = 0; i < steps_x - 1; i++) {
        x.push_back(-r + i * dx);
        cur_alfa = calc_alfa(x[i]);
        dr.push_back({});
        dR.push_back({});
        phi_x.push_back({});
        for (j = 0; j < steps_y - 1; j++) {
            cur_z = calc_zgri(x[i], y[j], calc_xg(hg_prev[j]), hg_prev[j]);
            dr[i].push_back(calc_dri(cur_z, cos(cur_alfa)));
            dR[i].push_back(calc_dR(dr[i][j], calc_vri(dr[i][j], cur_alfa, w, vx), cos(cur_alfa), cur_z, z_prev[i][j], dt));
            phi_x[i].push_back(calc_phi_x(vx, w, dr[i][j]));
        }
    }
    return solve_integral(c.first, c.second, x, y, phi_x, phi_y, dR, dr);
}

// Метод Рунге-Кутты 4-го порядка
void Wheel::modeling(double mk, double f_side, double t_start, double t_end, double dt) {
    double steps_t = (t_end - t_start) / dt;
    double vx_k1, vx_k2, vx_k3, vx_k4, vx_ans;
    double vy_k1, vy_k2, vy_k3, vy_k4, vy_ans;
    double w_k1, w_k2, w_k3, w_k4, w_ans;
    vector<double> integral;  // return {m_x_cos, m_x_sin, p_x_cos, p_x_sin, p_y_cos};
    for (int step = 0; step < steps_t - 1; step++) {  // Считаем следующий шаг
        
        WheelStepState last_step_state = states.back();
        integral = calc_contact_reactions(0, 0, 0, dt);
        vx_k1 = dt * (integral[2] - integral[3]) / m;
        vy_k1 = dt * (integral[4] + f_side) / m;
        w_k1 = dt * (mk - integral[0] - integral[1]) / j;
        integral = calc_contact_reactions(vx_k1 * 0.5, vy_k1 * 0.5, w_k1 * 0.5, dt);
        vx_k2 = dt * (integral[2] - integral[3]) / m;
        vy_k2 = dt * (integral[4] + f_side) / m;
        w_k2 = dt * (mk - integral[0] - integral[1]) / j;
        integral = calc_contact_reactions(vx_k2 * 0.5, vy_k2 * 0.5, w_k2 * 0.5, dt);
        vx_k3 = dt * (integral[2] - integral[3]) / m;
        vy_k3 = dt * (integral[4] + f_side) / m;
        w_k3 = dt * (mk - integral[0] - integral[1]) / j;
        integral = calc_contact_reactions(vx_k3, vy_k3, w_k3, dt);
        vx_k4 = dt * (integral[2] - integral[3]) / m;
        vy_k4 = dt * (integral[4] + f_side) / m;
        w_k4 = dt * (mk - integral[0] - integral[1]) / j;

        
        vx_ans = last_step_state.v_x + (1. / 6.) * (vx_k1 + 2. * vx_k2 + 2. * vx_k3 + vx_k4);
        vy_ans = last_step_state.v_y + (1. / 6.) * (vy_k1 + 2. * vy_k2 + 2. * vy_k3 + vy_k4);
        w_ans = last_step_state.w + (1. / 6.) * (w_k1 + 2. * w_k2 + 2. * w_k3 + w_k4);
        
        // cout << vx_ans << " " << vy_ans << " " << w_ans << endl;
        vector<double> x, z;
        // Перезаписываем для следующего шага координаты z на этом шагу
        // И одновременно записываем x и z для y = 0, x <= 0, чтобы посчитать длину пятна контакта
        for (int i = 0; i < steps_x - 1; i++) {
            if (i * dx <= r) {
                x.push_back(-r + i * dx);
                z.push_back(calc_zgri(-r + i * dx, 0, calc_xg(hg_max_prev), hg_max_prev));
            }
            for (int j = 0; j < steps_y - 1; j++) {
                z_prev[i][j] = calc_zgri(-r + i * dx, -b + j * dy, calc_xg(hg_prev[j]), hg_prev[j]);
            }
        }
        
        // Пересчитываем hg для следующего шага
        double qmax = calc_qmax(calc_f(calc_l1(), calc_l2(x, z)));
        // cout << qmax << endl;
        hg_max_prev = calc_hg(calc_q(qmax, b / 2.));
        // cout << hg_max_prev << endl;
        for (int j = 0; j < steps_y - 1; j++) {
            hg_prev[j] = calc_hg(calc_q(qmax, j * dy));
            // cout << "hg_prev: " << hg_prev[j] << endl;
        }
        // cout << endl;
        // WheelStepState(double x_, double y_, double z_, double vx, double vy, double vz, double ww)
        states.emplace_back(
            last_step_state.x_0 + (vx_ans * dt), last_step_state.y_0 + (vy_ans * dt), r - hg_max_prev,
            vx_ans, vy_ans, 0, w_ans
        );
    }
    
}

int main() 
{
    // c_gr из статьи, nu из статьи, phi_x наугад, phi_y наугад, a1 из статьи для рыхлого грунта, a2 наугад
    Ground ground(0.85, 0.77, 1, 1, 1, 1);
    Wheel wheel(0.3, 0.2, 0.1, 0.1, ground, 0.045, 1, 0.01, 0.01);
    wheel.modeling(50, 0, 0, 5, 1e-3);
    wheel.out();

}

    // Подсказки по заполнению полей:

    // Wheel(double wheel_radius, double wheel_width, double x, double y, 
    //     Ground ground, double wheel_j, double m, double delta_x, double delta_y);

    // struct WheelStepState {
    //     double x_0, y_0, z_0;  // Координаты центра колеса
    //     double v_x, v_y, v_z;  // Скорости центра колеса
    //     double w;
    //     WheelStepState(double x_, double y_, double z_, double vx, double vy, double vz, double ww)
    //         : x_0(x_), y_0(y_), z_0(z_), v_x(vx), v_y(vy), v_z(vz), w(ww) {};
    // };
    
    // states.emplace_back(x, y, gr.z_gr_nd(0, 0) + r, 0.1, 0.1, 0, 0.1);