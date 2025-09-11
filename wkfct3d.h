#include "olb3D.h"
#include "olb3D.hh"

// for Pedrizzetti 1996
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = FLOATING_POINT_TYPE;

inline T pedrizettiFlow(int iT, int iTperiod) {
   // Average flow: 0.4355
   // Max flow: ca. 1
   T pedrizzetiValue = 0.4355;
   pedrizzetiValue += 0.05 * cos(2. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue += 0.25 * sin(2. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue -= 0.13 * cos(4. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue += 0.13 * sin(4. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue -= 0.10 * cos(6. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue -= 0.02 * sin(6. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue -= 0.01 * cos(8. * M_PI * (T(iT) / T(iTperiod)));
   pedrizzetiValue -= 0.03 * sin(8. * M_PI * (T(iT) / T(iTperiod)));
   return pedrizzetiValue;
}

inline T wang2023Flow(int iT, int iTperiod) {
   // Average flow: 0.1
   // Max flow: ca. 0.35
   T Value = 0.1;
   Value += 0.01870 * cos(2. * M_PI * (T(iT) / T(iTperiod)));
   Value += 0.15871 * sin(2. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.09044 * cos(4. * M_PI * (T(iT) / T(iTperiod)));
   Value += 0.02530 * sin(4. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00407 * cos(6. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.02210 * sin(6. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00996 * cos(8. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00002 * sin(8. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00156 * cos(10. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00645 * sin(10. * M_PI * (T(iT) / T(iTperiod)));
   Value += 0.00161 * cos(12. * M_PI * (T(iT) / T(iTperiod)));
   Value -= 0.00274 * sin(12. * M_PI * (T(iT) / T(iTperiod)));

   return Value;
}

inline T corAflow(int iT, int iTperiod) {
   /*
   Needed average flow = 250 ml/min / Area
   = [(0.25 L/min)/(60 s/min * 1000 m^3/L)] / (2 * 3.141 * (0.0017647m)^2)
   = 0.21298502 m/s
   => * 1/0.6, because effective time with 0.21m/s is only 0.6 of the period
   */
   T flowVelocity = 0;
   if (T(iT % iTperiod) < (0.3 * T(iTperiod))) {
      flowVelocity = 0;
   } else if (T(iT % iTperiod) < (0.4 * T(iTperiod))) {
      flowVelocity = (T(iT % iTperiod) / T(iTperiod) - 0.3) * 2.12985 / 0.6;
   } else if (T(iT % iTperiod) < (0.9 * T(iTperiod))) {
      flowVelocity = 0.212985 / 0.6;
   } else {
      flowVelocity = (1. - T(iT % iTperiod) / T(iTperiod)) * 2.12985 / 0.6;
   }
   return flowVelocity;
}

inline T cosScale(int iT, int maxT) {
   return 0.5 * (1 - cos(M_PI * T(iT) / T(maxT)));
}

inline T linearScale(int iT, int maxT) { return T(iT) / T(maxT); }

inline T threeElmWindkessel(T pressure, T flow, T delta_flow, T alpha = 1,
                            T beta = 1, T gamma = 1) {
   return (alpha * flow + beta * delta_flow - gamma * pressure);
}

inline T rungeKutta(T pressure, T flow, T delta_flow, T delta_delta_flow,
                    T delta_t, T alpha = 1, T beta = 1, T gamma = 1) {
   T flow_half_future = flow + (delta_t * delta_flow / 2);
   T flow_future = flow + (delta_t * delta_flow);
   T delta_flow_half_future = delta_flow + (delta_t * delta_delta_flow / 2);
   T delta_flow_future = delta_flow + (delta_t * delta_delta_flow);

   T A = threeElmWindkessel(pressure, flow, delta_flow, alpha, beta, gamma);
   T pressure_A = pressure + (A * delta_t / 2);

   T B = threeElmWindkessel(pressure_A, flow_half_future,
                            delta_flow_half_future, alpha, beta, gamma);
   T pressure_B = pressure + (B * delta_t / 2);

   T C = threeElmWindkessel(pressure_B, flow_half_future,
                            delta_flow_half_future, alpha, beta, gamma);
   T pressure_C = pressure + (C * delta_t);

   T D = threeElmWindkessel(pressure_C, flow_future, delta_flow_future, alpha,
                            beta, gamma);

   return (pressure + ((A + 2 * B + 2 * C + D) / 6) * delta_t);
}

inline T backwardFiniteDifference(int derivative, int order, T deltaT, T a[],
                                  int size, bool zeroFirst = true) {
   // write values of a into a new array in the right order, depending on
   // zeroFirst
   T array[size];
   for (int i = 0; i < size; ++i) {
      if (!zeroFirst) {
         array[i] = a[(size - 1) - i];
      } else {
         array[i] = a[i];
      }
   }

   // third order approximations
   if (order == 3) {
      if ((derivative == 1) & (size > 2)) {
         return ((-11 / 16 * array[0]) + 3 * array[1] - (3 / 2 * array[2]) +
                 (1 / 3 * array[3])) /
                (deltaT);
      } else if ((derivative == 2) & (size > 3)) {
         return ((35 / 12 * array[0]) - (26 / 3 * array[1]) +
                 (19 / 2 * array[2]) - (14 / 3 * array[3]) +
                 (11 / 12 * array[4])) /
                (std::pow(deltaT, 2));
      } else {
         return NAN;
      }
   }
   // second order approximations
   else if (order == 2) {
      if ((derivative == 1) & (size > 1)) {
         return (1.5 * array[0] - 2 * array[1] + 0.5 * array[2]) / (deltaT);
      } else if ((derivative == 2) & (size > 2)) {
         return (2.0 * array[0] - 5 * array[1] + 4 * array[2] - array[3]) /
                (std::pow(deltaT, 2));
      } else {
         return NAN;
      }
   }
   // first order approximation
   else if (order == 1) {
      if ((derivative == 1) & (size > 0)) {
         return (array[0] - array[1]) / deltaT;
      } else if ((derivative == 2) & (size > 1)) {
         return (array[0] - 2 * array[1] + array[2]) / (std::pow(deltaT, 2));
      } else {
         return NAN;
      }
   } else {
      return NAN;
   }
}
