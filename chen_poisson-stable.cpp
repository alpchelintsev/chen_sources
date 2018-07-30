// Source code of the programs for the article "A New Reliable Numerical
// Method for Computing Chaotic Solutions of Dynamical Systems: The Chen Attractor Case"
// https://www.worldscientific.com/doi/abs/10.1142/S0218127415501874

// *******************************************************************
// * (c) 2015 Alexander Pchelintsev, pchelintsev.an@yandex.ru        *
// * The program for finding the time of return of trajectory        *
// * in the ρ-neighborhood of the starting point for the Chen system *
// * (with comments in Russian)                                      *
// *******************************************************************
// * Программа поиска времени возврата траектории системы Чена в     *
// * окрестность начальной точки                                     *
// *******************************************************************

#include <iostream>
#include <vector>
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
    bool fl_rp;
    do
    {
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
    }
    while(fl_rp);
}

int main()
{
    mpreal::set_default_prec(prec);
    mpreal x = "-10.33913519761", y = "-11.10031188035", z = "23.84877914089";
    mpreal x0 = x, y0 = y, z0 = z, T = "0.1", delta_t = "0.001";
    calc(x, y, z, T);
    while(true)
    {
        calc(x, y, z, delta_t);
        T += delta_t;
        mpreal rho = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
        cout << "\nrho = " << rho.toString() << "\nT = " << T.toString() << "\nx = " << x.toString() << "\ny = " << y.toString() << "\nz = " << z.toString() << endl;
        if(rho <= "0.6")
           break;
        //getchar();
    }
    return 0;
}
