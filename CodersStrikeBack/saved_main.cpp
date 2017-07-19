#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <chrono>
#include <thread>
#include <algorithm>
#include <sstream>

const float PI = 3.1415927f;

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

class Unit;

class Collision
{
public:
    Collision(Unit *_a, Unit *_b, float _t)
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

    bool operator==(const Collision &other) const
    {
        return other.a == this->a && other.b == this->b && other.t == this->t;
    }

    bool operator!=(const Collision &other) const
    {
        return !(other == (*this));
    }

    Unit *a;
    Unit *b;
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

    Unit(UnitType _type, int _id, float _r, float _vx, float _vy)
        : Point()
        , type(_type)
        , id(_id)
        , r(_r)
        , vx(_vx)
        , vy(_vy)
    {}

    virtual Collision collision(Unit &u)
    {
        if (u.id == this->id) { return Collision(); }

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

    virtual void bounce(Unit *u) = 0;
};

class Checkpoint : public Unit
{
public:
    Checkpoint(int _x, int _y, int _id)
        : Unit(UnitType::CHECKPOINT, id, 600, 0, 0)
    {
        x = _x;
        y = _y;
    }

    void bounce(Unit *u)
    {

    }

    int next_cp_id;
};

class Pod : public Unit
{
public:
    Pod();

    float angle;
    int next_check_pt_id;
    int num_of_cps_passed;
    int timeout;
    Pod *partner;
    bool shield;
    Checkpoint *next_checkpoint;

    float target_angle;
    float target_thrust;


    float getAngle(const Point &p);
    float diffAngle(const Point &p);
    void rotate(float _angle);
    void boost(int thrust);
    void move(float t);
    void end();
    virtual void bounce(Unit *u) override;
    void bounceWithCheckpoint(const Checkpoint *cp);
    Point getTarget(float delta_angle);
    float score();
    void output();
};

void Pod::output()
{
    Point p = getTarget(target_angle);

    if (!shield)
    {
        std::cout << static_cast<int>(p.x) << " " << static_cast<int>(p.y) << " " << static_cast<int>(target_thrust) << std::endl;
    }
    else
    {
        std::cout << static_cast<int>(p.x) << " " << static_cast<int>(p.y) << " " << "SHIELD" << std::endl;
    }
}

Point Pod::getTarget(float delta_angle)
{
    float a = angle + delta_angle;

    if (a >= 360.0f)
    {
        a = a - 360.0f;
    }
    else if (a < 0.0f)
    {
        a += 360.0f;
    }

    a = a * PI / 180.0;
    float px = x + cos(a) * 10000.0;
    float py = y + sin(a) * 10000.0;

    return Point(px, py);
}

Pod::Pod()
    : Unit(UnitType::POD, 0, 400.0f, 0.0f, 0.0f)
    , angle(0)
    , next_check_pt_id(0)
    , num_of_cps_passed(0)
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
void Pod::rotate(float _angle)
{

    float a = this->diffAngle(getTarget(_angle));

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

float Pod::score()
{
    return 50000.0f * num_of_cps_passed - distance2(*next_checkpoint);
}


static unsigned int g_seed;

// Used to seed the generator.           
inline void fast_srand(int seed) {
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline float fast_rand(float min, float max) {
    g_seed = (214013*g_seed+2531011);
    int result = (g_seed>>16)&0x7FFF;

    if (min < 0.0f)
    {
        if (result < 16384.0f)
        {
            return static_cast<float>(result / 16384.0f) * min;
        }
        else
        {
            return static_cast<float>(result - 16384.0f) / 16384.0f * max;
        }
    }
    else
    {
        return min + static_cast<float>(result / 32767.0f) * (max - min);
    }
}

class Gene
{
public:
    float thrust_gene[2];
    float angle_gene[2];
    float shield_gene[2];

    void randomize()
    {        
        thrust_gene[0] = fast_rand( 0.0f, 1.0f);
        thrust_gene[1] = fast_rand( 0.0f, 1.0f);
        angle_gene[0]  = fast_rand(-1.0f, 1.0f);
        angle_gene[1]  = fast_rand(-1.0f, 1.0f);
        shield_gene[0] = fast_rand( 0.0f, 1.0f);
        shield_gene[1] = fast_rand( 0.0f, 1.0f);
    }

    void mutate(float amplitude)
    {
        mutate_angle (amplitude, 0);
        mutate_thrust(amplitude, 0);
        mutate_shield(amplitude, 0);
        mutate_angle (amplitude, 1);
        mutate_thrust(amplitude, 1);
        mutate_shield(amplitude, 1);
    }

    void mutate_angle(float amplitude, int id)
    {        
        float a_min = angle_gene[id] - 1.0 * amplitude;
        float a_max = angle_gene[id] + 1.0 * amplitude;

        if(a_min < -1.0) { a_min = -1.0; }
        if(a_max >  1.0) { a_max =  1.0; }

        angle_gene[id] = fast_rand(a_min, a_max);
    }

    void mutate_thrust(float amplitude, int id)
    {
        float min = thrust_gene[id] - 1.0 * amplitude;
        float max = thrust_gene[id] + 1.0 * amplitude;

        if (min < 0.0f) { min = 0.0f; }
        if (max > 1.0f) { max = 1.0f; }

        thrust_gene[id] = fast_rand(min, max);
    }

    void mutate_shield(float amplitude, int id)
    {
        float min = shield_gene[id] - 1.0 * amplitude;
        float max = shield_gene[id] + 1.0 * amplitude;

        if (min < 0.0f) { min = 0.0f; }
        if (max > 1.0f) { max = 1.0f; }

        shield_gene[id] = fast_rand(min, max);
    }
};

class Genome
{
public:
    std::array<Gene, 3> genes;

    void randomize()
    {
        for (auto & gene : genes)
        {
            gene.randomize();
        }
    }

    int getThrust(int id, int turn)
    {
        int thrust = static_cast<int>(100.0f * genes[turn].thrust_gene[id]);
        if (thrust > 100) { return 100; }
        else if (thrust < 0) { return 0; }
        else { return thrust; }
    }

    float getAngle(int id, int turn)
    {
        float angle = genes[turn].angle_gene[id] * 18.0f;
        if (angle < -18.0f) { return -18.0f; }
        else if (angle > 18.0f) { return 18.0f; }
        else { return angle; }
    }

    bool getShield(int id, int turn)
    {
        return genes[turn].shield_gene[id] >= 0.95f;
    }


    void crossover(const Genome &p1, const Genome &p2)
    {
        genes[0].thrust_gene[0] = p1.genes[0].thrust_gene[0];
        genes[0].thrust_gene[1] = p2.genes[0].thrust_gene[1];

        genes[1].thrust_gene[0] = p2.genes[1].thrust_gene[0];
        genes[1].thrust_gene[1] = p1.genes[1].thrust_gene[1];

        genes[2].thrust_gene[0] = p1.genes[2].thrust_gene[0];
        genes[2].thrust_gene[1] = p2.genes[2].thrust_gene[1];
        
        genes[0].angle_gene[0] = p2.genes[0].angle_gene[0];
        genes[0].angle_gene[1] = p1.genes[0].angle_gene[1];

        genes[1].angle_gene[0] = p1.genes[1].angle_gene[0];
        genes[1].angle_gene[1] = p2.genes[1].angle_gene[1];

        genes[2].angle_gene[0] = p2.genes[2].angle_gene[0];
        genes[2].angle_gene[1] = p1.genes[2].angle_gene[1];

        genes[0].shield_gene[0] = p1.genes[0].shield_gene[0];
        genes[0].shield_gene[1] = p2.genes[0].shield_gene[1];

        genes[1].shield_gene[0] = p2.genes[1].shield_gene[0];
        genes[1].shield_gene[1] = p1.genes[1].shield_gene[1];

        genes[2].shield_gene[0] = p1.genes[2].shield_gene[0];
        genes[2].shield_gene[1] = p2.genes[2].shield_gene[1];
    }

    void mutate(float amplitude)
    {
        for (auto & gene : genes)
        {
            gene.mutate(amplitude);
        }
    }
};

class Solution
{
public:
    Genome my_genes;
    Genome opp_genes;

    std::vector<Pod> pods;
    float score;

    void randomize()
    {
        my_genes.randomize();
    }
    
    float evaluate(Pod my_pod_1, Pod my_pod_2, Pod opp_pod_1, Pod opp_pod_2, std::vector<Checkpoint*> &checkpoints)
    {
        pods.clear();
        pods.push_back(my_pod_1);
        pods.push_back(my_pod_2);
        pods.push_back(opp_pod_1);
        pods.push_back(opp_pod_2);

        for (int i = 0; i < my_genes.genes.size(); ++i)
        {
            playTurn(checkpoints, i);
        }

        score = evaluation();

        pods[0].target_angle = my_genes.getAngle(0, 0);
        pods[1].target_angle = my_genes.getAngle(1, 0);
        pods[0].target_thrust = my_genes.getThrust(0, 0);
        pods[1].target_thrust = my_genes.getThrust(1, 0);
        pods[0].shield = my_genes.getShield(0, 0);
        pods[1].shield = my_genes.getShield(1, 0);
        return score;
    }

    float evaluation()
    {
        return pods[0].score() + pods[1].score() - pods[2].score() - pods[3].score();
    }

    void output()
    {
        pods[0].output();
        pods[1].output();
    }

    void playTurn(std::vector<Checkpoint*> &checkpoints, int turn_num)
    {
        float t = 0.0f;
        
        pods[0].rotate(my_genes.getAngle(0, turn_num));
        pods[1].rotate(my_genes.getAngle(1, turn_num));
        pods[2].rotate(opp_genes.getAngle(0, turn_num));
        pods[3].rotate(opp_genes.getAngle(1, turn_num));

        pods[0].boost(my_genes.getThrust(0, turn_num));
        pods[1].boost(my_genes.getThrust(1, turn_num));
        pods[2].boost(opp_genes.getThrust(0, turn_num));
        pods[3].boost(opp_genes.getThrust(1, turn_num));

        pods[0].next_checkpoint = checkpoints[pods[0].next_check_pt_id];
        pods[1].next_checkpoint = checkpoints[pods[1].next_check_pt_id];
        pods[2].next_checkpoint = checkpoints[pods[2].next_check_pt_id];
        pods[3].next_checkpoint = checkpoints[pods[3].next_check_pt_id];

        Collision prev_found_collision;

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
                    temp_collision = pods[i].collision(pods[j]);

                    // If the collision is earlier than the current earliest collision, or if we havent found a collision yet.
                    if (!temp_collision.IsNull() && temp_collision.t + t < 1.0f && (first_collision.IsNull() || temp_collision.t < first_collision.t) && temp_collision != prev_found_collision && temp_collision.t != 0.0)
                    {
                        first_collision = std::move(temp_collision);
                    }

                }

                temp_collision = pods[i].collision(*checkpoints[pods[i].next_check_pt_id]);

                if (!temp_collision.IsNull() && temp_collision.t + t < 1.0f && (first_collision.IsNull() || temp_collision.t < first_collision.t) && temp_collision != prev_found_collision && temp_collision.t != 0.0)
                {
                    first_collision = std::move(temp_collision);
                }
            }
            
            if (first_collision.IsNull() || temp_collision == prev_found_collision || temp_collision.t == 0.0)
            {
                // No collision, update the m_pods
                for (int i = 0; i < pods.size(); ++i)
                {
                    pods[i].move(1.0 - t);
                }
                t = 1.0f;
            }
            else
            {
                // Move the pods up to the first collision.
                for (int i = 0; i < pods.size(); ++i)
                {
                    pods[i].move(first_collision.t);
                }

                first_collision.a->bounce(first_collision.b);

                t += first_collision.t;

                prev_found_collision = first_collision;
            }

            for (int i = 0; i < pods.size(); ++i)
            {
                pods[i].end();
            }
        }
    }

};

//typedef std::chrono::high_resolution_clock clock;


class Controller
{
public:
    Controller()        
        : m_my_pods()
        , m_my_pods_clone()
        , m_opp_pods()
        , m_opp_pods_clone()
        , m_pods()
        , m_checkpoints()
        , m_solutions()
        , m_winning_solution()
        , m_number_of_laps(0)
        , m_current_lap(0)
    {
        fast_srand(519881);
        m_winning_solution.randomize();

        m_solutions[0] = new Solution();
        m_solutions[1] = new Solution();
        m_solutions[2] = new Solution();

        m_my_pods[0].id = 100;
        m_my_pods[1].id = 101;
        m_opp_pods[0].id = 200;
        m_opp_pods[1].id = 201;

    }
    void setNumberOfLaps(int laps)
    {
        m_number_of_laps = laps;
    }

    void output()
    {
        m_solutions[0]->output();
    }

    void addCheckpoint(int x, int y)
    {
        m_checkpoints.push_back(new Checkpoint(x, y, m_checkpoints.size()));
        m_checkpoints.back()->next_cp_id = m_checkpoints.size();
    }

    void generatePopulation()
    {
        fast_srand(16519887);
        m_solutions[0]->randomize();
        m_solutions[1]->randomize();
        *m_solutions[2] = m_winning_solution;
    }  
    
    void findSolution()
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        generatePopulation();

        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        
        std::cerr << "Evaluating Solutions..." << std::endl;
        m_solutions[0]->evaluate(m_my_pods[0], m_my_pods[1], m_opp_pods[0], m_opp_pods[1], m_checkpoints);
        m_solutions[1]->evaluate(m_my_pods[0], m_my_pods[1], m_opp_pods[0], m_opp_pods[1], m_checkpoints);
        m_solutions[2]->evaluate(m_my_pods[0], m_my_pods[1], m_opp_pods[0], m_opp_pods[1], m_checkpoints);

        std::cerr << "Ranking Solutions..." << std::endl;
        rankSolutions();
        int num_of_evolutions = 0;
           
        std::cerr << "Evolving Solutions..." << std::endl; 
        float amplitude = 1.0;
        while (elapsed < 100.0)
        {
            fast_srand(elapsed);
            amplitude = 1.0; //std::min(std::max(1.0 - elapsed / 50.0, 0.01), 1.0);

            evolve(amplitude);
            rankSolutions();
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
            num_of_evolutions++;
        }
        std::cerr << elapsed << " elapsed." << std::endl;
        
        auto &gene = m_solutions[0]->my_genes.genes[0];
        
        std::cerr << num_of_evolutions << " evolutions." << std::endl;
        std::cerr << "Score: " << m_solutions[0]->score << std::endl;
        std::cerr << "Gene 0: " << gene.thrust_gene[0] << ", " << gene.angle_gene[0] << ", " << gene.shield_gene[0] << std::endl;        
        std::cerr << "Gene 1: " << gene.thrust_gene[1] << ", " << gene.angle_gene[1] << ", " << gene.shield_gene[1] << std::endl;
    }

    void rankSolutions()
    {
        int i = 0;
        i++;
        std::sort(m_solutions.begin(), m_solutions.end(), [](const Solution *a, const Solution *b) { return a->score > b->score; });
    }

    void evolve(float amplitude)
    {
        crossover();
        mutate(amplitude);
    }

    void mutate(float amplitude)
    {
        m_solutions[2]->my_genes.mutate(amplitude);
        m_solutions[2]->evaluate(m_my_pods[0], m_my_pods[1], m_opp_pods[0], m_opp_pods[1], m_checkpoints);
    }

    void crossover()
    {
        m_buffer = *m_solutions[2];
        m_buffer.my_genes.crossover(m_solutions[0]->my_genes, m_solutions[1]->my_genes);

        m_buffer.evaluate(m_my_pods[0], m_my_pods[1], m_opp_pods[0], m_opp_pods[1], m_checkpoints);
        checkForBetterSolution();
    }

    void checkForBetterSolution()
    {
        if (m_buffer.score > m_solutions[0]->score)
        {
            delete m_solutions[2];
            m_solutions[2] = m_solutions[1];
            m_solutions[1] = m_solutions[0];
            m_solutions[0] = new Solution(m_buffer);
        }
        else if (m_buffer.score > m_solutions[1]->score)
        {
            delete m_solutions[2];
            m_solutions[2] = m_solutions[1];
            m_solutions[1] = new Solution(m_buffer);
        }
        else if (m_buffer.score > m_solutions[2]->score)
        {
            delete m_solutions[2];
            m_solutions[2] = new Solution(m_buffer);
        }
    }

    float score()
    {
        
    }

    Pod & MyPod(int id)
    {
        return m_my_pods[id];
    }

    Pod & OppPod(int id)
    {
        return m_opp_pods[id];
    }

private:
    Pod m_my_pods[2];
    Pod m_my_pods_clone[2];
    Pod m_opp_pods[2];
    Pod m_opp_pods_clone[2];
    std::vector<Pod*> m_pods;
    std::vector<Checkpoint*> m_checkpoints;
    std::array<Solution*, 3> m_solutions;
    Solution m_buffer;
    Solution m_winning_solution;
    int m_number_of_laps;
    int m_current_lap;
    int m_ranked_sol_ids[3];
};

void test(Controller & ctrl)
{
    ctrl.setNumberOfLaps(3);
    ctrl.addCheckpoint( 3500, 5100);
    ctrl.addCheckpoint(13500, 7500);
    ctrl.addCheckpoint(12500, 1300);
    ctrl.addCheckpoint(10500, 5900);


    ctrl.MyPod(0).x = 3600;
    ctrl.MyPod(0).y = 4600;
    ctrl.MyPod(0).vx = 10;
    ctrl.MyPod(0).vy = 15;
    ctrl.MyPod(0).angle = 0;
    ctrl.MyPod(0).next_check_pt_id = 1;

    ctrl.MyPod(1).x = 3400;
    ctrl.MyPod(1).y = 5500;
    ctrl.MyPod(1).vx = 30;
    ctrl.MyPod(1).vy = 20;
    ctrl.MyPod(1).angle = 1;
    ctrl.MyPod(1).next_check_pt_id = 1;

    ctrl.OppPod(0).x = 3950;
    ctrl.OppPod(0).y = 6950;
    ctrl.OppPod(0).vx = 0;
    ctrl.OppPod(0).vy = 0;
    ctrl.OppPod(0).angle = 0;
    ctrl.OppPod(0).next_check_pt_id = 1;

    ctrl.OppPod(1).x = 3200;
    ctrl.OppPod(1).y = 6500;
    ctrl.OppPod(1).vx = 0;
    ctrl.OppPod(1).vy = 0;
    ctrl.OppPod(1).angle = 1;
    ctrl.OppPod(1).next_check_pt_id = 1;


    ctrl.findSolution();
    ctrl.findSolution();
    ctrl.findSolution();
    ctrl.findSolution();
    ctrl.findSolution();
    ctrl.output();
}

int main()
{
    Controller ctrl;
    //test(ctrl);
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();
    ctrl.setNumberOfLaps(laps);
    for (int i = 0; i < checkpointCount; i++) {
        int checkpointX;
        int checkpointY;
        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
        ctrl.addCheckpoint(checkpointX, checkpointY);
    }

    // game loop
    while (1) 
    {
        for (int i = 0; i < 2; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            std::cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; std::cin.ignore();

            ctrl.MyPod(i).x = x;
            ctrl.MyPod(i).y = y;
            ctrl.MyPod(i).vx = vx;
            ctrl.MyPod(i).vy = vy;
            ctrl.MyPod(i).angle = angle;
            ctrl.MyPod(i).next_check_pt_id = nextCheckPointId;
        }
        for (int i = 0; i < 2; i++) {
            int x2; // x position of the opponent's pod
            int y2; // y position of the opponent's pod
            int vx2; // x speed of the opponent's pod
            int vy2; // y speed of the opponent's pod
            int angle2; // angle of the opponent's pod
            int nextCheckPointId2; // next check point id of the opponent's pod
            std::cin >> x2 >> y2 >> vx2 >> vy2 >> angle2 >> nextCheckPointId2; std::cin.ignore();

            ctrl.OppPod(i).x = x2;
            ctrl.OppPod(i).y = y2;
            ctrl.OppPod(i).vx = vx2;
            ctrl.OppPod(i).vy = vy2;
            ctrl.OppPod(i).angle = angle2;
            ctrl.OppPod(i).next_check_pt_id = nextCheckPointId2;
        }

        ctrl.findSolution();

        ctrl.output();
    }
}
