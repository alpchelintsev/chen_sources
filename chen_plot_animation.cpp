// Source code of the programs for the article "A New Reliable Numerical
// Method for Computing Chaotic Solutions of Dynamical Systems: The Chen Attractor Case"
// https://www.worldscientific.com/doi/abs/10.1142/S0218127415501874

// ************************************************************************
// * (c) 2015 Alexander Pchelintsev, pchelintsev.an@yandex.ru             *
// * The program for constructing of animation the arc of the trajectory  *
// * of the Chen system in a given time interval, animation in GIF format *
// * (with comments in Russian). The program generates the text file      *
// * anim_gnuplot.txt with commands for Gnuplot.                          *
// ************************************************************************
// * Программа построения анимации дуги траектории динамической системы   *
// * Чена на заданном отрезке времени (файл в формате GIF).               *
// * Программа генерирует текстовый файл anim_gnuplot.txt с командами     *
// * для Gnuplot.                                                         *
// ************************************************************************

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;

#include <mpreal.h>
using namespace mpfr;

#include <cstdio>

// Точность по степенному ряду
#define eps_c   "1e-80"

// Параметры системы уравнений Чена
#define a   35
#define b   3
#define c   28

// Количество бит под мантиссу вещественного числа
#define prec    300

// Как считать шаг: 1 - по оценке отрезка сходимости, 0 - заданная величина
#define FL_CALC 1

// Фиксированный шаг по времени
#define step_t  "0.02"

#define STEP_PLOT 3

// Функция возвращает длину отрезка сходимости степенного ряда
mpreal get_delta_t(const mpreal &alpha0, const mpreal &beta0, const mpreal &gamma0)
{
    mpreal h2 = (mpfr::max)((mpfr::max)(fabs(alpha0), fabs(beta0)), fabs(gamma0)),
           h1 = (mpfr::max)((mpfr::max)(2*fabs(a), fabs(c-a) + 2*h2 + fabs(c)), fabs(b) + 2*h2);
    mpreal h3 = h2 >= 1 ? h1 * h2 : (mpfr::max)((mpfr::max)(2*fabs(a), fabs(c-a) + fabs(c) + 1), fabs(b) + 1);
    return 1/(h3 + "1e-10");
}

// Функция вычисления значений фазовых координат в конечный момент времени
// x, y и z - координаты начальной точки; T - длина отрезка интегрирования;
// way - направление поиска решений: 1 - вперед по времени, -1 - назад по времени
void calc(mpreal &x, mpreal &y, mpreal &z, const mpreal &T, int way = 1)
{
    mpreal t = 0, delta_t, L, p, s1, s2;
    int iter = 0;
    ofstream f("anim_gnuplot.txt");
    f << "set term gif animate optimize delay 1 size 800, 800 crop font \"Times-Roman,13\"\nset output \"chen_traj_anim.gif\"\n";
    f << "set xlabel \"x\"\nset ylabel \"y\"\nset zlabel \"z\"\n";
    f << "set xrange [-40:40]\nset yrange [-40:40]\nset zrange [0:60]\n";
    ostringstream buf;
    buf << x.toString(3) << " " << y.toString(3) << " " << z.toString(3) << endl;
    bool fl_rp;
    do
    {
        if(!(iter % STEP_PLOT))
        {
            if(iter)
                f << "e\n";
            f << "splot '-' with line\n";
            f << buf.str();
        }

        if(FL_CALC)
            delta_t = get_delta_t(x, y, z);
        else
            delta_t = step_t;
        t += delta_t;
        if(t < T)
            fl_rp = true;
        else if(t > T)
        {
            delta_t -= t-T;
            fl_rp = false;
        }
        else
            fl_rp = false;

        vector<mpreal> alpha, beta, gamma;
        alpha.push_back(x);
        beta.push_back(y);
        gamma.push_back(z);

        int i = 0;
        L = sqrt(alpha[0]*alpha[0] + beta[0]*beta[0] + gamma[0]*gamma[0]);
        p = way * delta_t;
        while(L > eps_c)
        {
            // Вычисляем новые коэффициенты степенных рядов
            s1 = s2 = 0;
            for(int j = 0; j <= i; j++)
            {
                s1 += alpha[j] * gamma[i-j];
                s2 += alpha[j] * beta[i-j];
            }
            alpha.push_back(a*(beta[i] - alpha[i])/(i+1));
            beta.push_back(((c-a)*alpha[i] - s1 + c*beta[i])/(i+1));
            gamma.push_back((s2 - b*gamma[i])/(i+1));

            i++;

            x += alpha[i] * p;
            y += beta[i] * p;
            z += gamma[i] * p;
            L = fabs(p) * sqrt(alpha[i]*alpha[i] + beta[i]*beta[i] + gamma[i]*gamma[i]);
            p *= way * delta_t;
        }
        iter++;
        if(!(iter % STEP_PLOT))
            buf << x.toString(3) << " " << y.toString(3) << " " << z.toString(3) << endl;
    }
    while(fl_rp);
    f.close();
}

int main()
{
    mpreal::set_default_prec(prec);
    mpreal T;
    cout << "\nВведите длину отрезка времени > ";
    cin >> T;

    mpreal x, y, z;
    cout << "\nx0 > ";
    cin >> x;

    cout << "y0 > ";
    cin >> y;

    cout << "z0 > ";
    cin >> z;
    cout << endl;

    calc(x, y, z, T);
    return 0;
}
