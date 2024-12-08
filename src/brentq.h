#pragma once

typedef double (*callback_type)(double, void *);

double brentinv(callback_type f, double y, void *func_data);
double brentq(callback_type f, double y, double xa, double xb, double xtol, double rtol, int iter, void *func_data);
