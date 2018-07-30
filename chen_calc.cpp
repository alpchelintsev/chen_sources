// Source code of the programs for the article "A New Reliable Numerical
// Method for Computing Chaotic Solutions of Dynamical Systems: The Chen Attractor Case"
// https://www.worldscientific.com/doi/abs/10.1142/S0218127415501874

// ******************************************************************
// * (c) 2015 Alexander Pchelintsev, pchelintsev.an@yandex.ru       *
// * The program implementing the algorithm shown in Figure 1       *
// * of this article for the Chen system (with comments in Russian) *
// ******************************************************************
// * Программа, реализующая алгоритм, показанный на рисунке 1       *
// * статьи, для системы Чена                                       *
// ******************************************************************

#include <iostream>
#include <vector>
using namespace std;

#include <mpreal.h>
using namespace mpfr;

// Точность по степенному ряду
#define eps_c   "1e-53"

// Параметры системы уравнений Чена
#define a   35
#define b   3
#define c   28

// Количество бит под мантиссу вещественного числа
#define prec    189

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
    cout << "\nВ начальный момент времени:\n\ndx/dt = " << a*(y-x) << "\ndy/dt = " <<
            (c-a)*x-x*z+c*y << "\ndz/dt = " << x*y-b*z << endl;

    mpreal t = 0, delta_t, L, p, s1, s2, tmin, tmax, dmin, dmax, tdmin, tdmax;
    int l = 0, imax, imin, lmin, lmax = 0, lmindt = 0, lmaxdt = 0;
    bool fl_rp;
    do
    {
        l++;
        if(FL_CALC)
            delta_t = get_delta_t(x, y, z);
        else
            delta_t = step_t;
        bool flagimin = false;
        if(t == 0)
        {
            imax = 0;
            flagimin = true;
            tdmax = tdmin = tmax = 0;
            dmax = dmin = delta_t;
        }
        else if(delta_t < dmin)
        {
             dmin = delta_t;
             tdmin = t;
             lmindt = l;
        }
        else if(delta_t > dmax)
        {
             dmax = delta_t;
             tdmax = t;
             lmaxdt = l;
        }

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

        if(delta_t < dmin)
        {
             dmin = delta_t;
             tdmin = T-delta_t;
             lmindt = l;
        }
        else if(delta_t > dmax)
        {
             dmax = delta_t;
             tdmax = T-delta_t;
             lmaxdt = l;
        }

        vector<mpreal> alpha, beta, gamma;
        alpha.push_back(x);
        beta.push_back(y);
        gamma.push_back(z);

        L = sqrt(alpha[0]*alpha[0] + beta[0]*beta[0] + gamma[0]*gamma[0]);
        int i = 0;
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
        if(i > imax)
        {
            imax = i;
            lmax = l;
            if(t >= T)
                tmax = T-delta_t;
            else
                tmax = t-delta_t;
        }
        if(flagimin)
        {
            imin = i;
            lmin = l;
            tmin = 0;
        }
        else if(i < imin)
        {
            imin = i;
            lmin = l;
            if(t >= T)
                tmin = T-delta_t;
            else
                tmin = t-delta_t;
        }
    }
    while(fl_rp);

    cout << "\nКоординаты в конечный момент времени:\nx = " << x.toString() << "\ny = " <<
            y.toString() << "\nz = " << z.toString() << endl;
    cout << "\nЗначения производных:\n\ndx/dt = " << a*(y-x) << "\ndy/dt = " <<
            (c-a)*x-x*z+c*y << "\ndz/dt = " << x*y-b*z << endl;
    cout << "\nmin Степень полиномов = " << imin << endl;
    cout << "Соответствующий момент времени = " << tmin << endl;
    cout << "Разница = " << T-tmin << endl;
    cout << "Номер момента времени = " << lmin << endl;
    cout << "\nmax Степень полиномов = " << imax << endl;
    cout << "Соответствующий момент времени = " << tmax << endl;
    cout << "Разница = " << T-tmax << endl;
    cout << "Номер момента времени = " << lmax << endl;
    cout << "\nmin delta_t = " << dmin << endl;
    cout << "Соответствующий момент времени = " << tdmin << endl;
    cout << "Разница = " << T-tdmin << endl;
    cout << "Номер момента времени = " << lmindt << endl;
    cout << "\nmax delta_t = " << dmax << endl;
    cout << "Соответствующий момент времени = " << tdmax << endl;
    cout << "Разница = " << T-tdmax << endl;
    cout << "Номер момента времени = " << lmaxdt << endl;
    
    cout << "\nN = " << l << endl;
}

int main()
{
    mpreal::set_default_prec(prec);
    cout << "Машинный эпсилон = " << machine_epsilon() << endl;

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
    cout << "\n\n*** Проход назад ***\n";
    calc(x, y, z, T, -1);

    return 0;
}
