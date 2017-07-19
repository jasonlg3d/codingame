#if !defined CHECKPOINT_H
#define CHECKPOINT_H

#include "Unit.h"

class Checkpoint : public Unit
{
public:
    Checkpoint()
    {
        type = UnitType::CHECKPOINT;
    }

    void bounce(Unit u)
    {

    }

    int next_cp_id;
};

#endif
