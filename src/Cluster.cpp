#include "Cluster.h"
#include <cstdlib>
#include <cmath>

double Cluster::dpd(double x1, double y1, double x2, double y2) const
{
    double dx = std::abs(x1-x2);
    double dy = std::abs(y1-y2);

    if (dx < dy)
        return (2 * dy - dx);
    else
        return (2 * dx - dy);
}

int Cluster::dpd(int x1, int y1, int x2, int y2)
{
    int dx = std::abs(x1-x2);
    int dy = std::abs(y1-y2);

    if (dx < dy)
        return (2 * dy - dx);
    else
        return (2 * dx - dy);
}

double Cluster::dpd() const
{
    return dpd(begin_x, begin_y, end_x, end_y);
}

int Cluster::kspd(const int& x1, const int& y1, const int &x2, const int& y2)
{
    int dx = std::abs(x1-x2);
    int dy = std::abs(y1-y2);

    return ((dx>dy) ? dx : dy);
}
