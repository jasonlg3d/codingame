#include "Pod.h"

const float PI = 3.1415927f;

Pod::Pod()
    : angle(0)
    , next_check_pt_id(0)
    , checked(0)
    , timeout(0)
    , partner()
    , shield(false)
{
    type = UnitType::POD;
}

float Pod::getAngle(const Point &p)
{
    float d = this->distance(p);
    float dx = (p.x - this->x) / d;
    float dy = (p.y - this->y) / d;

    float a = std::acos(dx) * 180.0f / PI;

    if (dy < 0) {
        a = 360.0f - a;
    }

    return a;
}
float Pod::diffAngle(const Point &p)
{
    float a = this->getAngle(p);
    float right = this->angle <= a ? a - this->angle : 360.0f - this->angle + a;
    float left  = this->angle >= a ? this->angle - a : this->angle + 360.0f - a;

    if (right < left)
    {
        return right;
    }
    else
    {
        return - left;
    }
}
void Pod::rotate(const Point &p)
{
    float a = this->diffAngle(p);

    if      (a >  18.0f) { a =  18.0f; }
    else if (a < -18.0f) { a = -18.0f; }

    this->angle += a;

    if (this->angle >= 360.0f)
    {
         this->angle = this->angle - 360.0f;
    }
    else if (this->angle < 0.0f)
    {
        this->angle += 360.0f;
    }
}
void Pod::boost(int thrust)
{
    if (this->shield)
    {
        return;
    }

    float ra = this->angle * PI / 180.0f;

    this->vx += std::cos(ra) * thrust;
    this->vy += std::sin(ra) * thrust;
}

void Pod::move(float t)
{
    this->x += this->vx * t;
    this->y += this->vy * t;
}

void Pod::end()
{
    this->x = std::round(this->x);
    this->y = std::round(this->y);
    this->vx = static_cast<int>(this->vx * 0.85f);
    this->vx = static_cast<int>(this->vy * 0.85f);
}
void Pod::play(const Point &p, int thrust)
{
    this->rotate(p);
    this->boost(thrust);
    this->move(1.0f);
    this->end();
}

void Pod::bounce(Unit *u)
{
    if (u->type == UnitType::CHECKPOINT)
    {
        this->bounceWithCheckpoint(dynamic_cast<Checkpoint*>(u));
    }
    else
    {
        Pod * pod2 = dynamic_cast<Pod*>(u);
        // Shield mass == 1
        float m1 = this->shield ? 10.0f : 1.0f;
        float m2 = pod2->shield ? 10.0f : 1.0f;
        float mcoeff = (m1 + m2) / (m1 * m2);

        float nx = this->x - pod2->x;
        float ny = this->y - pod2->y;

        float nxnysquare = nx*nx + ny*ny;

        float dvx = this->vx - pod2->vy;
        float dvy = this->vy - pod2->vy;

        float product = nx*dvx + ny*dvy;
        float fx = (nx * product) / (nxnysquare * mcoeff);
        float fy = (ny * product) / (nxnysquare * mcoeff);

        // Apply the impact vectors
        this->vx -= fx / m1;
        this->vy -= fy / m1;
        pod2->vx += fx / m2;
        pod2->vy += fy / m2;

        // if the nomr of the impact vector is less than 120, we normalize it to 120
        float impulse = sqrt(fx*fx + fy*fy);
        if (impulse < 120.0f)
        {
            fx = fx * 120.0f / impulse;
            fy = fy * 120.0f / impulse;
        }

        // Apply the impact vector a second time
        this->vx -= fx / m1;
        this->vy -= fy / m1;

        pod2->vx += fx / m2;
        pod2->vy += fy / m2;
    }
}

void Pod::bounceWithCheckpoint(const Checkpoint *cp)
{
    next_check_pt_id = cp->next_cp_id;
    this->timeout = 100;
}

void Pod::output(Move move)
{

}