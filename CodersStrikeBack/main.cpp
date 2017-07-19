#include <cmath>
#include <vector>
#include <array>
#include <iostream>
#include <chrono>
#include <thread>
#include <algorithm>
#include <sstream>

const float PI = 3.1415927f;
const int NUM_OF_TURNS_TO_SIM = 6;

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

    float distance(const Point &p) const
    {
        return std::sqrt((this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y));
    }
    float distance2(const Point &p) const
    {
        return (this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y);
    }

    float angleTo(const Point &target)
    {
        float d = distance(target);
        float dx = (target.x - this->x) / d;
        float dy = (target.y - this->y) / d;

        float a = std::acos(dx) * 180.0f / PI;

        if (dy < 0) {
            a = 360.0f - a;
        }

        return a;
    }

    Point closest(const Point &a, const Point &b) const
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

Point get_target_for_angle_change(const Point &loc, float curr_angle, float delta_angle) 
{
    float a = curr_angle + delta_angle;

    if (a >= 360.0f)
    {
        a = a - 360.0f;
    }
    else if (a < 0.0f)
    {
        a += 360.0f;
    }

    a = a * PI / 180.0;
    float px = loc.x + cos(a) * 5000.0;
    float py = loc.y + sin(a) * 5000.0;

    return Point(px, py);
}

class AI
{
public:
    AI()
        : m_in()
        , m_out()
    {}

    struct Input
    {
        int turn_num;
        Point loc;            // position of your pod
        int vx;               // x speed of your pod
        int vy;               // y speed of your pod
        int angle;            // angle of your pod
        int nextCheckPointId; // next check point id of your pod
        Point next_cp;
    };

    struct Output
    {
        Point target;
        int thrust;
        bool shield;

        Output()
            : target()
            , thrust(0)
            , shield(false)
        {}
    };

    Input & In() { return m_in; }
    const Output & Out(int pair_id) const { return m_out[pair_id]; }

    virtual void Update() = 0;
protected:
    Input m_in;
    Output  m_out[2];
};

class VectorAI : public AI
{
public:
    VectorAI()
        : AI()
    {}

    void Update()
    {
        m_out[0].target = m_in.next_cp;
        m_out[0].thrust = 80;
        m_out[0].shield = false;

        m_out[1].target = m_in.next_cp;
        m_out[1].thrust = 80;
        m_out[1].shield = false;
    }
};

class Gene
{
public:
    void SetThrustGene(float val, int id)
    {
        m_genes[id] = val;
    }

    float ThrustGene(int id) const
    {
        return m_genes[id];
    }

    float AngleGene(int id) const
    {
        return m_genes[2 + id];
    }

    void SetAngleGene(float val, int id)
    {
        m_genes[2 + id] = val;
    }

    float ShieldGene(int id) const
    {
        return m_genes[4 + id];
    }

    void SetShieldGene(float val, int id)
    {
        m_genes[4 + id] = val;
    }

    void randomize()
    {        

        m_genes[0] = fast_rand( 0.0f, 1.0f);
        m_genes[1] = fast_rand( 0.0f, 1.0f);
        m_genes[2] = fast_rand(-1.0f, 1.0f);
        m_genes[3] = fast_rand(-1.0f, 1.0f);
        m_genes[4] = fast_rand( 0.0f, 1.0f);
        m_genes[5] = fast_rand( 0.0f, 1.0f);
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

private:
    float m_genes[6]; // 0,1 - Thrust | 2,3 - Angle | 4,5 - Sheild

    void mutate_angle(float amplitude, int id)
    {        
        float a_min = m_genes[2 + id] - 1.0 * amplitude;
        float a_max = m_genes[2 + id] + 1.0 * amplitude;

        if(a_min < -1.0) { a_min = -1.0; }
        if(a_max >  1.0) { a_max =  1.0; }

        m_genes[2 + id] = fast_rand(a_min, a_max);
    }

    void mutate_thrust(float amplitude, int id)
    {
        float min = m_genes[id] - 1.0 * amplitude;
        float max = m_genes[id] + 1.0 * amplitude;

        if (min < 0.0f) { min = 0.0f; }
        if (max > 1.0f) { max = 1.0f; }

        m_genes[id] = fast_rand(min, max);
    }

    void mutate_shield(float amplitude, int id)
    {
        float min = m_genes[4 + id] - 1.0 * amplitude;
        float max = m_genes[4 + id] + 1.0 * amplitude;

        if (min < 0.0f) { min = 0.0f; }
        if (max > 1.0f) { max = 1.0f; }

        m_genes[4 + id] = fast_rand(min, max);
    }
};

class Genome
{
public:
    Gene genes;

    void randomize()
    {
        genes.randomize();
    }

    int getThrust(int id, int turn) const
    {
        int thrust = static_cast<int>(100.0f * genes.ThrustGene(id));
        if (thrust > 85) { return 100; }
        else if (thrust < 15) { return 15; }
        else { return thrust; }
    }

    float getAngle(int id, int turn) const
    {
        float angle = genes.AngleGene(id) * 18.0f;
        if (angle < -18.0f) { return -18.0f; }
        else if (angle > 18.0f) { return 18.0f; }
        else { return angle; }
    }

    bool getShield(int id, int turn) const
    {
        return genes.ShieldGene(id) >= 0.95f;
    }


    Genome crossover(const Genome &p1, const Genome &p2)
    {
        Genome result;
        
        if (fast_rand(0.0f, 10.0f) > 5.0f)
        {
            result.genes.SetThrustGene(p1.genes.ThrustGene(0), 0);
            result.genes.SetThrustGene(p1.genes.ThrustGene(1), 1);

            result.genes.SetAngleGene(p2.genes.ThrustGene(0), 0);
            result.genes.SetAngleGene(p2.genes.ThrustGene(1), 1);

            result.genes.SetShieldGene(p1.genes.ThrustGene(0), 0);
            result.genes.SetShieldGene(p1.genes.ThrustGene(1), 1);
        }
        else
        {
            result.genes.SetThrustGene(p2.genes.ThrustGene(0), 0);
            result.genes.SetThrustGene(p2.genes.ThrustGene(1), 1);

            result.genes.SetAngleGene(p1.genes.ThrustGene(0), 0);
            result.genes.SetAngleGene(p1.genes.ThrustGene(1), 1);

            result.genes.SetShieldGene(p2.genes.ThrustGene(0), 0);
            result.genes.SetShieldGene(p2.genes.ThrustGene(1), 1);
        }

        return result;
    }

    void mutate(float amplitude)
    {
        genes.mutate(amplitude);
    }
};

class GeneticAI : public AI
{
public:
    GeneticAI()
        : AI()
        , m_genome()
    {}
    
    void LoadGenome(const Genome &g)
    {
        m_genome = g;
    }

    Genome GetGenome() const { return m_genome; }

    void Update()
    {
        m_out[0].shield = getShield(0);
        m_out[0].target = getTarget(0);
        m_out[0].thrust = getThrust(0);
                
        m_out[1].shield = getShield(1);
        m_out[1].target = getTarget(1);
        m_out[1].thrust = getThrust(1);
    }

    void Randomize()
    {
        m_genome.randomize();
    }

    void Mutate(float amplitude)
    {
        m_genome.mutate(amplitude);
    }

    void SetScore(float val) { m_score = val; }
    float Score() const { return m_score; }

    Genome BreedWith(const GeneticAI *partner)
    {
        return m_genome.crossover(m_genome, partner->GetGenome());
    }

private:
    
    Point getTarget(int pair_id)
    {
        return get_target_for_angle_change(Point(m_in.loc.x, m_in.loc.y), m_in.angle, m_genome.getAngle(pair_id, m_in.turn_num));
    }

    int getThrust(int pair_id) 
    {
        return m_genome.getThrust(pair_id, m_in.turn_num);
    }

    bool getShield(int pair_id) 
    {
        return m_genome.getShield(pair_id, m_in.turn_num);
    }

    Genome m_genome;
    float m_score;
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
    Pod(int pair_id);

    struct State
    {
        int angle;
        Point loc;
        Point next_cp_loc;
        int next_cp_id;
        int vx;
        int vy;
        bool shield;
    };

    const AI* podAI() const { return m_ai; }

    void ReloadState()
    {
        m_state = m_initial_state;
    }

    void SaveState()
    {
        m_initial_state = m_state;
    }

    void Update(int turn_num);

    void SetAI(AI *ai)
    {
        m_ai = ai;
    }

    float diffAngle(const Point &p);
    void rotate(const Point &target);
    void boost(int thrust);
    void move(float t);
    void end();
    virtual void bounce(Unit *u) override;
    void bounceWithCheckpoint(const Checkpoint *cp);
    float score();
    void output();

    State & St() { return m_state; }

private:
    int m_num_of_cps_passed;
    int m_timeout;
    int m_pair_id;
    AI *m_ai;
    State m_state;
    State m_initial_state;

};

void Pod::Update(int turn_num)
{
    m_ai->In().angle = m_state.angle;
    m_ai->In().loc = m_state.loc;
    m_ai->In().nextCheckPointId = m_state.next_cp_id;
    m_ai->In().next_cp.x = m_state.next_cp_loc.x;
    m_ai->In().next_cp.y = m_state.next_cp_loc.y;
    m_ai->In().turn_num = turn_num;
    m_ai->In().vx = m_state.vx;
    m_ai->In().vy = m_state.vy;
    m_ai->Update();

    rotate(m_ai->Out(m_pair_id).target);
    boost(m_ai->Out(m_pair_id).thrust);
    m_state.shield = m_ai->Out(m_pair_id).shield;
    
}

void Pod::output()
{
    Point p = m_ai->Out(m_pair_id).target;

    if (!m_ai->Out(m_pair_id).shield)
    {
        std::cout << static_cast<int>(p.x) << " " << static_cast<int>(p.y) << " " << static_cast<int>(m_ai->Out(m_pair_id).thrust) << std::endl;
    }
    else
    {
        std::cout << static_cast<int>(p.x) << " " << static_cast<int>(p.y) << " " << "SHIELD" << std::endl;
    }
}

Pod::Pod(int pair_id)
    : Unit(UnitType::POD, 0, 400.0f, 0.0f, 0.0f)
    , m_num_of_cps_passed(0)
    , m_timeout(0)
    , m_pair_id(pair_id)
    , m_ai(nullptr)
    , m_state()
    , m_initial_state()
{
    type = UnitType::POD;
}

float Pod::diffAngle(const Point &p)
{
    float a = m_state.loc.angleTo(p);
    float right = m_state.angle <= a ? a - m_state.angle : 360.0f - m_state.angle + a;
    float left  = m_state.angle >= a ? m_state.angle - a : m_state.angle + 360.0f - a;

    if (right < left)
    {
        return right;
    }
    else
    {
        return - left;
    }
}
void Pod::rotate(const Point &target)
{

    float a = this->diffAngle(target);

    if      (a >  18.0f) { a =  18.0f; }
    else if (a < -18.0f) { a = -18.0f; }

    m_state.angle += a;

    if (m_state.angle >= 360.0f)
    {
         m_state.angle = m_state.angle - 360.0f;
    }
    else if (m_state.angle < 0.0f)
    {
        m_state.angle += 360.0f;
    }
}
void Pod::boost(int thrust)
{
    if (m_state.shield)
    {
        return;
    }

    float ra = m_state.angle * PI / 180.0f;

    m_state.vx += std::cos(ra) * thrust;
    m_state.vy += std::sin(ra) * thrust;
}

void Pod::move(float t)
{
    m_state.loc.x += m_state.vx * t;
    m_state.loc.y += m_state.vy * t;
}

void Pod::end()
{
    m_state.loc.x = std::round(m_state.loc.x);
    m_state.loc.y = std::round(m_state.loc.y);
    m_state.vx = static_cast<int>(m_state.vx * 0.85f);
    m_state.vx = static_cast<int>(m_state.vy * 0.85f);
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
        float m1 = m_state.shield ? 10.0f : 1.0f;
        float m2 = pod2->St().shield ? 10.0f : 1.0f;
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
        m_state.vx -= fx / m1;
        m_state.vy -= fy / m1;
        pod2->St().vx += fx / m2;
        pod2->St().vy += fy / m2;

        // if the nomr of the impact vector is less than 120, we normalize it to 120
        float impulse = sqrt(fx*fx + fy*fy);
        if (impulse < 120.0f)
        {
            fx = fx * 120.0f / impulse;
            fy = fy * 120.0f / impulse;
        }

        // Apply the impact vector a second time
        m_state.vx -= fx / m1;
        m_state.vy -= fy / m1;

        pod2->St().vx += fx / m2;
        pod2->St().vy += fy / m2;
    }
}

void Pod::bounceWithCheckpoint(const Checkpoint *cp)
{
    if (cp->id == m_state.next_cp_id)
    {
        m_num_of_cps_passed++;
    }
    m_state.next_cp_id = cp->next_cp_id;
    m_timeout = 100;
}

float Pod::score()
{
    float score = 50000.0f * m_num_of_cps_passed - m_state.loc.distance2(m_state.next_cp_loc) * 0.01;
    return score;
}

struct PodID
{
    enum
    {
        MY_POD_1,
        MY_POD_2,
        OPP_POD_1,
        OPP_POD_2
    };
};

class Controller
{
public:
    Controller()        
        : m_pods()
        , m_checkpoints()
        , m_number_of_laps(0)
        , m_current_lap(0)
    {
        fast_srand(519881);

        m_buffer_ai = new GeneticAI();

        m_basic_ai[0] = new VectorAI();
        m_basic_ai[1] = new VectorAI();

        m_pods.push_back(Pod(0));
        m_pods.push_back(Pod(1));
        m_pods.push_back(Pod(0));
        m_pods.push_back(Pod(1));

        m_pods[PodID::MY_POD_1].id = 100;
        m_pods[PodID::MY_POD_2].id = 101;
        m_pods[PodID::OPP_POD_1].id = 200;
        m_pods[PodID::OPP_POD_2].id = 201;

        m_pods[PodID::MY_POD_1].SetAI(m_basic_ai[0]);
        m_pods[PodID::MY_POD_2].SetAI(m_basic_ai[1]);
        m_pods[PodID::OPP_POD_1].SetAI(m_basic_ai[0]);
        m_pods[PodID::OPP_POD_2].SetAI(m_basic_ai[1]);

        for (int i = 0; i < m_ai_solutions.size(); ++i)
        {
            m_ai_solutions[i] = new GeneticAI();
            m_ai_solutions[i]->Randomize();
        }

    }
    void SetNumberOfLaps(int laps)
    {
        m_number_of_laps = laps;
    }

    void Update()
    {
        for (auto &pod : m_pods)
        {
            pod.SaveState();
        }
                
        ScoreAI(m_ai_solutions[0]);
        ScoreAI(m_ai_solutions[1]);
        ScoreAI(m_ai_solutions[2]);
        
        SortAIs();

        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();

        float amplitude = 1.0;
        double time_to_evolve_me = 140.0;
        double time_to_evolve_them = 20.0;
        
        int num_of_generations = 0;

        while (elapsed < time_to_evolve_me)
        {
            fast_srand(static_cast<int>(now));
            amplitude = std::min(std::max(1.0 - elapsed / time_to_evolve_me, 0.01), 0.3);
            Evolve(amplitude);
            ++num_of_generations;
            elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        }
        
        for (auto &pod : m_pods)
        {
            pod.ReloadState();
        }

        std::cerr << "POD 1 Angle: " << m_pods[0].St().angle << "°" << std::endl;
        std::cerr << "POD 2 Angle: " << m_pods[1].St().angle << "°" << std::endl;
        std::cerr << "POD 1 Loc: " << m_pods[0].St().loc.x << "," << m_pods[0].St().loc.y << std::endl;
        std::cerr << "POD 2 Loc: " << m_pods[1].St().loc.x << "," << m_pods[1].St().loc.y << std::endl;
        std::cerr << "Number of Generations: " << num_of_generations << std::endl;
        std::cerr << "Winning Score: " << m_ai_solutions[0]->Score() << std::endl;
        std::cerr << "Solution: " << m_ai_solutions[0]->GetGenome().genes.AngleGene(0) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().genes.ThrustGene(0) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().genes.ShieldGene(0) << std::endl
                  << "          " << m_ai_solutions[0]->GetGenome().getAngle(0, 0) << "°, "
                  << "          " << m_ai_solutions[0]->GetGenome().getThrust(0, 0) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().getShield(0, 0) << std::endl
                  << "          " << m_ai_solutions[0]->GetGenome().genes.AngleGene(1) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().genes.ThrustGene(1) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().genes.ShieldGene(1) << std::endl
                  << "          " << m_ai_solutions[0]->GetGenome().getAngle(1, 0) << "°, "
                  << "          " << m_ai_solutions[0]->GetGenome().getThrust(1, 0) << ", "
                  << "          " << m_ai_solutions[0]->GetGenome().getShield(1, 0) << std::endl;
        m_pods[PodID::MY_POD_1].SetAI(m_ai_solutions[0]);
        m_pods[PodID::MY_POD_2].SetAI(m_ai_solutions[0]);
        
        m_pods[PodID::MY_POD_1].Update(0);
        m_pods[PodID::MY_POD_2].Update(0);
    }

    void Evolve(float amplitude)
    {

        Mutate(0, amplitude * 0.1);
        Mutate(1, amplitude * 0.1);
        Mutate(2, amplitude * 0.1);
        Mutate(3, amplitude);
        Mutate(4, amplitude);
        Mutate(5, amplitude);
        
        Breed(0, 1);
        Breed(0, 2);
        Breed(0, 3);
        Breed(0, 4);
        Breed(0, 5);
        Breed(1, 2);
        Breed(1, 3);
        Breed(1, 4);
        Breed(1, 5);

        SortAIs();
    }

    void Mutate(int id, float amplitude)
    {
        m_ai_solutions[id]->Mutate(amplitude);
        ScoreAI(m_ai_solutions[id]);
    }

    void ScoreAI(GeneticAI *ai)
    {
        for (auto & pod : m_pods)
        {
            pod.ReloadState();
        }

        m_pods[PodID::MY_POD_1].SetAI(ai);
        m_pods[PodID::MY_POD_2].SetAI(ai);

        for (int i = 0; i < NUM_OF_TURNS_TO_SIM; ++i)
        {
            SimulateTurn(i);
        }

        ai->SetScore(Score());
    }

    void Breed(int id_1, int id_2)
    {
        m_buffer_ai->LoadGenome(m_ai_solutions[id_1]->BreedWith(m_ai_solutions[id_2]));

        ScoreAI(m_buffer_ai);

        if (m_buffer_ai->Score() > m_ai_solutions.back()->Score())
        {
            m_ai_solutions.back()->LoadGenome(m_buffer_ai->GetGenome());
        }

        SortAIs();
    }

    void SortAIs()
    {
        std::sort(m_ai_solutions.begin(), m_ai_solutions.end(), [](const GeneticAI *a, const GeneticAI *b) { return a->Score() > b->Score(); });
    }

    void Output()
    {
        m_pods[PodID::MY_POD_1].output();
        m_pods[PodID::MY_POD_2].output();
    }

    void AddCheckpoint(int x, int y)
    {
        m_checkpoints.push_back(new Checkpoint(x, y, m_checkpoints.size()));
        m_checkpoints.back()->next_cp_id = m_checkpoints.size();
    }
    
    Pod & MyPod(int id)
    {
        return m_pods[PodID::MY_POD_1 + id];
    }

    Pod & OppPod(int id)
    {
        return m_pods[PodID::OPP_POD_1 + id];
    }

private:   

    float Score()
    {
        return 1.5 * (m_pods[PodID::MY_POD_1].score() + m_pods[PodID::MY_POD_2].score()) -
               0.0 * (m_pods[PodID::OPP_POD_1].score() + m_pods[PodID::OPP_POD_2].score());
    }

    void SimulateTurn(int turn_num)
    {
        float t = 0.0f;

        for (int i = 0; i < m_pods.size(); ++i)
        {
            m_pods[i].St().next_cp_loc.x = m_checkpoints[m_pods[i].St().next_cp_id]->x;
            m_pods[i].St().next_cp_loc.y = m_checkpoints[m_pods[i].St().next_cp_id]->y;
        }
        
        for (auto &pod : m_pods)
        {
            pod.Update(turn_num);
        }

        Collision prev_found_collision;

        while (t < 1.0)
        {
            Collision first_collision;
            Collision temp_collision;

            // Find all the collisions that are going to occur
            for (int i = 0; i < m_pods.size(); ++i)
            {
                /// Collision with another pod?
                for (int j = 0; j < m_pods.size(); ++j)
                {
                    temp_collision = m_pods[i].collision(m_pods[j]);

                    // If the collision is earlier than the current earliest collision, or if we havent found a collision yet.
                    if (!temp_collision.IsNull() && temp_collision.t + t < 1.0f && (first_collision.IsNull() || temp_collision.t < first_collision.t) && temp_collision != prev_found_collision && temp_collision.t != 0.0)
                    {
                        first_collision = std::move(temp_collision);
                    }

                }

                temp_collision = m_pods[i].collision(*m_checkpoints[m_pods[i].St().next_cp_id]);

                if (!temp_collision.IsNull() && temp_collision.t + t < 1.0f && (first_collision.IsNull() || temp_collision.t < first_collision.t) && temp_collision != prev_found_collision && temp_collision.t != 0.0)
                {
                    first_collision = std::move(temp_collision);
                }
            }
            
            if (first_collision.IsNull() || temp_collision == prev_found_collision || temp_collision.t == 0.0)
            {
                // No collision, update the m_pods
                for (int i = 0; i < m_pods.size(); ++i)
                {
                    m_pods[i].move(1.0 - t);
                }
                t = 1.0f;
            }
            else
            {
                // Move the pods up to the first collision.
                for (int i = 0; i < m_pods.size(); ++i)
                {
                    m_pods[i].move(first_collision.t);
                }

                first_collision.a->bounce(first_collision.b);

                t += first_collision.t;

                prev_found_collision = first_collision;
            }

            for (int i = 0; i < m_pods.size(); ++i)
            {
                m_pods[i].end();
            }
        }
    }

    std::vector<Pod> m_pods;
    std::vector<Checkpoint*> m_checkpoints;
    int m_number_of_laps;
    int m_current_lap;
    int m_ranked_sol_ids[3];
    std::array<VectorAI*, 2> m_basic_ai;
    std::array<GeneticAI*, 6> m_ai_solutions;
    std::array<GeneticAI*, 6> m_opponent_ai_solutions;
    GeneticAI* m_buffer_ai;
};

void test(Controller & ctrl)
{
    ctrl.SetNumberOfLaps(3);
    ctrl.AddCheckpoint( 3500, 5100);
    ctrl.AddCheckpoint(13500, 7500);
    ctrl.AddCheckpoint(12500, 1300);
    ctrl.AddCheckpoint(10500, 5900);


    ctrl.MyPod(0).St().loc.x      = 3600;
    ctrl.MyPod(0).St().loc.y      = 4600;
    ctrl.MyPod(0).St().vx         = 10;
    ctrl.MyPod(0).St().vy         = 15;
    ctrl.MyPod(0).St().angle      = 0;
    ctrl.MyPod(0).St().next_cp_id = 1;

    ctrl.MyPod(1).St().loc.x      = 3400;
    ctrl.MyPod(1).St().loc.y      = 5500;
    ctrl.MyPod(1).St().vx         = 30;
    ctrl.MyPod(1).St().vy         = 20;
    ctrl.MyPod(1).St().angle      = 1;
    ctrl.MyPod(1).St().next_cp_id = 1;

    ctrl.OppPod(0).St().loc.x      = 3950;
    ctrl.OppPod(0).St().loc.y      = 6950;
    ctrl.OppPod(0).St().vx         = 0;
    ctrl.OppPod(0).St().vy         = 0;
    ctrl.OppPod(0).St().angle      = 0;
    ctrl.OppPod(0).St().next_cp_id = 1;

    ctrl.OppPod(1).St().loc.x      = 3200;
    ctrl.OppPod(1).St().loc.y      = 6500;
    ctrl.OppPod(1).St().vx         = 0;
    ctrl.OppPod(1).St().vy         = 0;
    ctrl.OppPod(1).St().angle      = 1;
    ctrl.OppPod(1).St().next_cp_id = 1;

    ctrl.Update();
    ctrl.Output();
}

int main()
{
    Controller ctrl;
    //test(ctrl);
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();
    ctrl.SetNumberOfLaps(laps);
    for (int i = 0; i < checkpointCount; i++) {
        int checkpointX;
        int checkpointY;
        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
        ctrl.AddCheckpoint(checkpointX, checkpointY);
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

            ctrl.MyPod(i).St().loc.x      = x;
            ctrl.MyPod(i).St().loc.y      = y;
            ctrl.MyPod(i).St().vx         = vx;
            ctrl.MyPod(i).St().vy         = vy;
            ctrl.MyPod(i).St().angle      = angle;
            ctrl.MyPod(i).St().next_cp_id = nextCheckPointId;
        }
        for (int i = 0; i < 2; i++) {
            int x2; // x position of the opponent's pod
            int y2; // y position of the opponent's pod
            int vx2; // x speed of the opponent's pod
            int vy2; // y speed of the opponent's pod
            int angle2; // angle of the opponent's pod
            int nextCheckPointId2; // next check point id of the opponent's pod
            std::cin >> x2 >> y2 >> vx2 >> vy2 >> angle2 >> nextCheckPointId2; std::cin.ignore();

            ctrl.OppPod(i).St().loc.x      = x2;
            ctrl.OppPod(i).St().loc.y      = y2;
            ctrl.OppPod(i).St().vx         = vx2;
            ctrl.OppPod(i).St().vy         = vy2;
            ctrl.OppPod(i).St().angle      = angle2;
            ctrl.OppPod(i).St().next_cp_id = nextCheckPointId2;
        }
        
        ctrl.Update();

        ctrl.Output();
    }
}
