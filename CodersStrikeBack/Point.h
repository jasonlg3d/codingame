#if !defined POINT_H
#define POINT_H
#include <cmath>

class Point
{
public:
    Point()
        : x(0.0f)
        , y(0.0f)
    {}

    Point(float _x, float _y)
        : x(_x)
        , y(_y)
    {}

    float x;
    float y;

    float distance(const Point &p)
    {
        return std::sqrt((this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y));
    }
    float distance2(const Point &p)
    {
        return (this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y);
    }

    Point closest(const Point &a, const Point &b)
    {
        float da = b.y - a.y;
        float db = a.x - b.x;

        float c1 = da*a.x + db*a.y;
        float c2 = -db*this->x + da* this->y;
        float det = da*da + db*db;
        float cx = 0.0f;
        float cy = 0.0f;

        if (det != 0.0)
        {
            cx = (da*c1 - db*c2) / det;
            cy = (da*c2 - db*c1) / det;
        }
        else
        {
            cx = this->x;
            cy = this->y;
        }

        return Point(cx, cy);
    }
};

#endif

