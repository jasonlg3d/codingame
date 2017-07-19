#if !defined POD_H
#define POD_H

#include "Checkpoint.h"
#include "Unit.h"

class Pod : public Unit
{
public:
    Pod();

    float angle;
    int next_check_pt_id;
    int checked;
    int timeout;
    Pod partner;
    bool shield;

    float getAngle(const Point &p);
    float diffAngle(const Point &p);
    void rotate(const Point &p);
    void boost(int thrust);
    void move(float t);
    void end();
    void play(const Point &p, int thrust);
    void bounce(Unit *u);
    void output(Move move);
    void bounceWithCheckpoint(const Checkpoint *cp)
}

#endif
