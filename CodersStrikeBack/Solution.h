#if !defined SOLUTION_H
#define SOLUTION_H

#include "Pod.h"
#include "Checkpoint.h"
#include "Collision.h"
#include <vector>

class Solution
{
    void randomize();

    void play(std::vector<Pod*> &pods, std::vector<Checkpoint*> &cps)
    {
        float t = 0.0f;

        while (t < 1.0)
        {
            Collision first_collision;
            Collision temp_collision;

            // Find all the collisions that are going to occur
            for (int i = 0; i < pods.size(); ++i)
            {
                /// Collision with another pod?
                for (int j = 0; j < pods.size(); ++j)
                {
                    temp_collision = pods[i]->collision(*pods[j]);

                    // If the collision is earlier than the current earliest collision, or if we havent found a collision yet.
                    if (!temp_collision.IsNull() && temp_collision->t + t < 1.0f && (first_collision.IsNull() || temp_collision->t < first_collision->t))
                    {
                        first_collision = std::move(temp_collision);
                    }

                }

                temp_collision = pods[i]->collision(*cps[pods[i]->next_check_pt_id]);

                if (!temp_collision.IsNull() && temp_collision->t + t < 1.0f && (first_collision.IsNull() || temp_collision->t < first_collision->t))
                {
                    first_collision = std::move(temp_collision);
                }
            }
            
            if (first_collision.IsNull())
            {
                // No collision, update the pods
                for (int i = 0; i < pods.size(); ++i)
                {
                    pods[i]->move(1.0 - t);
                }
                t = 0.0f;
            else
            {
                // Move the pods up to the first collision.
                for (int i = 0; i < pods.size(); ++i)
                {
                    pods[i]->move(first_collision->t);
                }

                first_collision->a->bounce(first_collision->b);

                t += first_collision->t;                
            }

            for (int i = 0; i < pods.size(); ++i)
            {
                pods[i]->end();
            }
        }
    }
};

#endif
