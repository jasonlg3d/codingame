#if !defined UNIT_H
#define UNIT_H

#include "Point.h"
class Unit;

class Collision
{
public:
    Collision(const Unit *_a, const Unit *_b, float _t)
        : a(_a)
        , b(_b)
        , t(_t)
    {}

    Collision()
        : a(nullptr)
        , b(nullptr)
        , t(2.0)
    {}

    bool IsNull() { return a == nullptr || b == nullptr; }

    const Unit *a;
    const Unit *b;
    float t;
};


class Unit : public Point
{
public:
    enum class UnitType { POD, CHECKPOINT };
    UnitType type;
    int id;
    float r;
    float vx;
    float vy;

    virtual Collision collision(const Unit &u)
    {
        // square of the distance
        float dist = this->distance2(u);

        // sum of the radii squared
        float sr = (this->r + u.r) * (this->r + u.r);

        // Using everythign squared to avoid sqrt.

        if (dist < sr)
        {
            // objects are touching
            return Collision(this, &u, 0.0f);
        }

        // moving in parallel will never collide
        if (this->vx == u.vx && this->vy == u.vy)
        {
             return Collision();
        }

        // reference frame of u. U is stationary, while we move.
        float x = this->x - u.x;
        float y = this->y - u.y;
        Point myp = Point(x, y);
        float vx = this->vx - u.vx;
        float vy = this->vy - u.vy;
        Point up = Point(0.0f, 0.0f);

        // Find the closest point to u on the line that is our velocity vector.
        Point p = up.closest(myp, Point(x + vx, y + vy));

        // squre the distance between u and the closest point
        float pdist = up.distance2(p);

        // square of the distance between u and the closest point
        float mypdist = myp.distance2(p);

        if (pdist < sr)
        {
            // our speed on the line
            float length = std::sqrt(vx*vx + vy*vy);

            // move along the line to find the point of impact.
            float backdist = std::sqrt(sr - pdist);
            p.x = p.x - backdist * (vx / length);
            p.y = p.y - backdist * (vy / length);

            // if the point is now further away, we are moving away and therefore no collision will occur.
            if (myp.distance2(p) > mypdist)
            {
                return Collision();
            }

            pdist = p.distance(myp);

            // Impact further than possible travel in one turn.
            if (pdist > length)
            {
                return Collision();
            }

            float t = pdist / length;

            return Collision(this, &u, t);
        }

        return Collision();
    }

    virtual void bounce(Unit u)
    {

    }
};

#endif 
