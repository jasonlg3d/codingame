#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <chrono>
#include <numeric>

const double PI = 3.1415927;
const double DEGREES_TO_RAD = PI / 180.0;
const double RAD_TO_DEGREES = 180.0 / PI;

const double SHIP_RADIUS = 100.0;

static unsigned int g_seed;
// Used to seed the generator.           
inline void fast_srand(int seed) {
    g_seed = seed;
}

// Compute a pseudorandom integer.
// Output value in range [0, 32767]
inline double fast_rand(double min, double max) {
    g_seed = (214013*g_seed+2531011);
    int result = (g_seed>>16)&0x7FFF;

    if (min < 0.0)
    {
        if (result < 16384.0f)
        {
            return static_cast<double>(result / 16384.0) * min;
        }
        else
        {
            return static_cast<double>(result - 16384.0) / 16384.0f * max;
        }
    }
    else
    {
        return min + static_cast<double>(result / 32767.0) * (max - min);
    }
}

class Vector
{
public:
    Vector()
        : dx(0.0)
        , dy(0.0)
    {}

    Vector(double _dx, double _dy)
        : dx(_dx)
        , dy(_dy)
    {}

    double dx;
    double dy;

    Vector rotate(double angle_rad)
    {
        return Vector(
            this->dx * std::cos(angle_rad) - this->dy * std::sin(angle_rad),
            this->dx * std::sin(angle_rad) + this->dy * std::cos(angle_rad)
        );
    }

    Vector operator+(const Vector &v) const
    {
        return Vector(this->dx + v.dx, this->dy + v.dy);
    }

    Vector operator-(const Vector &v) const
    {
        return Vector(this->dx - v.dx, this->dy - v.dy);
    }

    Vector operator*(double m) const
    {
        return Vector(this->dx * m, this->dy *m);
    }
};

class Point
{
public:
    Point()
        : x(0.0)
        , y(0.0)
    {}

    Point(double _x, double _y)
        : x(_x)
        , y(_y)
    {}

    double x;
    double y;

    Point operator+(const Vector &v) const
    {
        return Point(this->x + v.dx, this->y + v.dy);
    }

    Point operator+(const Point &pt) const
    {
        return Point(this->x + pt.x, this->y + pt.y);
    }
    Point operator-(const Point &pt) const
    {
        return Point(this->x - pt.x, this->y - pt.y);
    }

    double distanceTo(const Point &pt) const
    {
        double dx = this->x - pt.x;
        double dy = this->y - pt.y;
        return std::sqrt(dx * dx + dy * dy);
    }
};
 
class Line
{
public:
    Line()
        : pts()
    {}

    bool IsHorizontalAt(double loc) const
    {
        for (int i = 0; i < pts.size() - 1; ++i)
        {
            if (pts.at(i).x + 5 < loc && loc < pts.at(i + 1).x - 5)
            {
                return pts.at(i).y == pts.at(i + 1).y;
            }
        }
        return false;
    }

    std::pair<Point, Point> GetSegmentFor(double loc) const
    {
        for (int i = 0; i < pts.size() - 1; ++i)
        {
            if (pts.at(i).x < loc && loc < pts.at(i + 1).x)
            {
                return std::make_pair(pts.at(i), pts.at(i + 1));
            }
        }
    }

    double GetYAtX(double loc) const
    {
        for (int i = 0; i < pts.size() - 1; ++i)
        {
            if (pts.at(i).x < loc && loc < pts.at(i + 1).x)
            {
                const Point &first = pts.at(i);
                const Point &second = pts.at(i + 1);

                return first.y + (loc - first.x) * (second.y - first.y) / (second.x - first.x);
            }
        }

        return 0.0;
    }

    std::vector<Point> pts;
};

class Particle
{
public:
    void Accelerate(const Vector & accel_vec, double time)
    {
        pos.x += velocity.dx * time + 0.5 * accel_vec.dx * time * time;
        pos.y += velocity.dy * time + 0.5 * accel_vec.dy * time * time;
        velocity.dx += accel_vec.dx * time;
        velocity.dy += accel_vec.dy * time;

        //pos.x = static_cast<int>(pos.x);
        //pos.y = static_cast<int>(pos.y);
        //velocity.dx = static_cast<int>(velocity.dx);
        //velocity.dy = static_cast<int>(velocity.dy);
    }

    Vector velocity;
    Point pos;
};

Vector GRAVITY(0.0, -3.711);
double MAX_X = 6999.0;
double MIN_X = 0.0;

struct ControlCommand
{
    ControlCommand(int ang, int pwr, int dur)
        : angle(ang)
        , power(pwr)
        , duration(dur)
    {}

    int angle;
    int power;
    int duration;
};


struct State
{
    State(int f, int pwr, int ang, Particle ldr)
        : fuel(f)
        , power(pwr)
        , angle(ang)
        , lander(ldr)
    {}

    State()
        : fuel(0)
        , power(0)
        , angle(0)
        , lander()
    {}

    int fuel;
    int power;
    int angle;
    Particle lander;
};

enum FlyState
{
    FLYING,
    CRASHED,
    LANDED
};

struct Lander
{
    Lander(const State& is, const std::vector<ControlCommand> &cmds, const Line &gnd)
        : initial_state(is)
        , commands(cmds)
        , ground(gnd)
        , fly_state(FLYING)
    {}

    Lander()
        : initial_state()
        , commands()
        , ground()
        , fly_state(FLYING)
    {}

    void ComputeTrajectory(double time)
    {
        State next_state = initial_state;
        for (int i = 0; i < commands.size(); ++i)
        {
            ControlCommand &cmd = commands.at(i);
            int new_ang = next_state.angle;
            int new_pwr = next_state.power;
    
            double elapsed = 0.0;
            int new_fuel = next_state.fuel;
            Particle new_particle = next_state.lander;
            while (elapsed <= cmd.duration)
            {
                new_ang = next_state.angle + std::min(std::max((cmd.angle - next_state.angle), -15), 15);
                new_pwr = next_state.power + std::min(std::max((cmd.power - next_state.power), -1), 1);
                Vector thrust_accel = (Vector(0.0, 1.0) * new_pwr).rotate(new_ang * DEGREES_TO_RAD);
                Vector accel = thrust_accel + GRAVITY;
                new_particle.Accelerate(accel, time);
                new_fuel = next_state.fuel - new_pwr;
                elapsed += time;
                next_state.angle = new_ang;
                next_state.fuel = new_fuel;
                next_state.power = new_pwr;
                next_state.lander = new_particle;                
                trajectory.push_back(next_state);
    
                if (EvaluateOutside(next_state)) { return; }
                if (EvaluateHitTheGround(next_state)) { return; }
                if (EvaluateNoFuel(next_state)) { return; }
            }
        }
    }

    State ComputeNextState(const State &st, const ControlCommand &cmd, double time)
    {
        int new_ang = st.angle;
        int new_pwr = st.power;

        double elapsed = 0.0;
        int new_fuel = st.fuel;
        Particle new_particle = st.lander;
        State next_state = st;
        while (elapsed <= cmd.duration)
        {
            new_ang = next_state.angle + std::min(std::max((cmd.angle - next_state.angle), -15), 15);
            new_pwr = next_state.power + std::min(std::max((cmd.power - next_state.power), -1), 1);
            Vector thrust_accel = (Vector(0.0, 1.0) * new_pwr).rotate(new_ang * DEGREES_TO_RAD);
            Vector accel = thrust_accel + GRAVITY;
            new_particle.Accelerate(accel, time);
            new_fuel = next_state.fuel - new_pwr;
            elapsed += time;
            next_state.angle = new_ang;
            next_state.fuel = new_fuel;
            next_state.power = new_pwr;
            next_state.lander = new_particle;
            
            trajectory.push_back(next_state);
        }
        return State(new_fuel, new_pwr, new_ang, new_particle);
    }

    bool EvaluateOutside(const State &st)
    {
        if (MIN_X > st.lander.pos.x || st.lander.pos.x > MAX_X)
        {
            fly_state = CRASHED;
            return true;
        }
        return false;
    }

    bool EvaluateHitTheGround(const State &st)
    {
        double height_at = ground.GetYAtX(st.lander.pos.x);
        bool is_flat = ground.IsHorizontalAt(st.lander.pos.x);
        double clearance = is_flat ? 0.0 : 100.0;
        if (height_at + clearance >= st.lander.pos.y)
        {
            if (st.angle == 0 &&
                st.lander.velocity.dy > -39.0 &&
                std::abs(st.lander.velocity.dx) <= 19.0 &&
                ground.IsHorizontalAt(st.lander.pos.x))
            {
                fly_state = LANDED;
            }
            else
            {
                fly_state = CRASHED;
            }
            return true;
        }
        return false;
    }

    bool EvaluateNoFuel(const State &st)
    {
        if (st.fuel <= 0)
        {
            fly_state = CRASHED;
            return true;
        }
        return false;
    }
    std::vector<ControlCommand> commands;
    State initial_state;
    Line ground;
    std::vector<State> trajectory;
    FlyState fly_state;
};

const int GENERATION_COUNT   = 400;
const int POPULATION_SIZE    = 20;
const int GENOME_SIZE        = 80;
const double UNIFORM_RATE    = 0.5;
const double MUTATION_RATE   = 0.06;
const double SELECTION_RATIO = 0.2;
struct Gene
{
    double power;
    double angle;
    double duration;
};

struct GenomeAndCommands
{
    GenomeAndCommands()
        : genes()
        , commands()
        , fitness_score(0.0)
        , lander()
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        fast_srand(static_cast<int>(now));

        for(int i = 0; i < GENOME_SIZE; ++i)
        {
            Gene gene;
            gene.angle = fast_rand(-90.0, 90.0);
            gene.power = fast_rand(0.0, 5.0);
            gene.duration = fast_rand(0.0, 20.0);

            genes.push_back(gene);
        }
    }

    void GenerateCommandsFromGenes()
    {
        commands.clear();
        for (auto gene : genes)
        {
            commands.push_back(ControlCommand(static_cast<int>(gene.angle), static_cast<int>(gene.power), static_cast<int>(gene.duration)));
            if (commands.back().power == 5) { commands.back().power = 4; }

        }
    }


    std::vector<Gene> genes;
    std::vector<ControlCommand> commands;
    double fitness_score;
    Lander lander;
};


struct Autopilot
{
    Autopilot()
        : init_state()
        , ground()
        , dt(1.0)
        , landing_elev(0.0)
    {}

    void Init(const State &is)
    {
        init_state = is;
    }

    ControlCommand CommandAtTime(int t)
    {
        
    }

    void FindSolution()
    {
        CalculateLandingZone();
        std::vector<GenomeAndCommands> population;
        for (int i = 0; i < POPULATION_SIZE; ++i)
        {
            population.push_back(GenomeAndCommands());
        }

        SortPopulationByFitness(population);

        std::vector<GenomeAndCommands> next_pop;
        double mutation_amplitude = 0.5;
        double mutation_rate = 0.8;
        for (int i = 0; i < GENERATION_COUNT; ++i)
        {
            next_pop.clear();

            /// Top 10% get saved.
            double best_score = population.at(0).fitness_score;

            for (auto & g : population)
            {
                if (next_pop.size() >= 5) { break; }
                if (best_score >= 0.0)
                {
                    if (g.fitness_score > best_score * 0.5)
                    {
                        next_pop.push_back(g);
                    }
                }
                else
                {
                    if (g.fitness_score > best_score * 1.5)
                    {
                        next_pop.push_back(g);
                    }
                }
            }

            if (next_pop.size() < 2)
            {
                next_pop.push_back(population.at(1));
            }

            while (next_pop.size() < POPULATION_SIZE)
            {
                int g1_idx = Select(population);
                int g2_idx = Select(population);
                if (g1_idx == g2_idx)
                {
                    g2_idx++;
                }
                next_pop.push_back(Crossover(population.at(g1_idx), population.at(g2_idx)));
                Mutate(next_pop.back(), mutation_rate, mutation_amplitude);
                EvaluateGenome(next_pop.back());
            }
            
            if (best_score > 10000)
            {
                mutation_amplitude = 0.05;
                mutation_rate = 0.5;
            }
            else
            {
                mutation_rate = std::max(0.6, mutation_rate - (static_cast<double>(i) / static_cast<double>(GENERATION_COUNT)) * 0.1);
                mutation_amplitude = std::max(0.08, mutation_amplitude - (static_cast<double>(i) / static_cast<double>(GENERATION_COUNT)) * 0.1);
            }
            
            population = std::move(next_pop);
            std::sort(population.begin(), population.end(), [](const GenomeAndCommands &lh, const GenomeAndCommands &rh) { return lh.fitness_score > rh.fitness_score; });
            double total_fitness = std::accumulate(population.begin(), population.end(), 0.0, [](double tot, const GenomeAndCommands &g) { return tot + g.fitness_score; });
            double best_xvel = population.at(0).lander.trajectory.back().lander.velocity.dx;
            double best_yvel = population.at(0).lander.trajectory.back().lander.velocity.dy;
            double best_xpos = population.at(0).lander.trajectory.back().lander.pos.x;
            double best_ypos = population.at(0).lander.trajectory.back().lander.pos.y;
            double best_ang  = population.at(0).lander.trajectory.back().angle;
            double best_fuel = population.at(0).lander.trajectory.back().fuel;
            total_fitness = total_fitness;
        }

        for (int i = 0; i < population.at(0).commands.size(); ++i)
        {
            ControlCommand &cmd = population.at(0).commands.at(i);
            double elapsed = 0.0;
            while (elapsed <= cmd.duration)
            {
                commands.push_back(cmd);
                elapsed += 1.0;
            }
        }

        std::cerr << population.at(0).fitness_score << std::endl;
        for(auto & t : population.at(0).lander.trajectory)
        {
            std::cerr << "h @ " << t.lander.pos.x << " = " << ground.GetYAtX(t.lander.pos.x) << std::endl;
        }
        //std::cin.ignore();
    }

    void PrintCommand(int i)
    {
        std::cout << commands.at(i).angle << " " << commands.at(i).power << std::endl;
    }

    int Select(std::vector<GenomeAndCommands> &pop)
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        fast_srand(static_cast<int>(now));
        for (int i = 0; i < pop.size() - 1; ++i)
        {
            double rand = fast_rand(0.0, 1.0);

            if (rand <= SELECTION_RATIO * (pop.size() - i) / pop.size())
            {
                return i;
            }
        }
        return 0;
    }

    GenomeAndCommands Crossover(const GenomeAndCommands &g1, const GenomeAndCommands &g2)
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        fast_srand(static_cast<int>(now));
        GenomeAndCommands offspring;
        offspring.genes.clear();
        for (int i = 0; i < g1.genes.size(); ++i)
        {
            double rand = fast_rand(0.0, 1.0);
            if (rand <= UNIFORM_RATE)
            {
                offspring.genes.push_back(g1.genes.at(i));
            }
            else
            {
                offspring.genes.push_back(g2.genes.at(i));
            }
        }

        return offspring;
    }

    void Mutate(GenomeAndCommands &genome, double rate, double amplitude)
    {
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        fast_srand(static_cast<int>(now));
        for (auto & gene : genome.genes)
        {
            double rand = fast_rand(0.0, 1.0);

            if (rand <= rate)
            {
                double ang_min = std::max(-90.0, gene.angle - 90.0 * amplitude);
                double ang_max = std::min(90.0, gene.angle + 90.0 * amplitude);
                gene.angle = fast_rand(ang_min, ang_max);

                double pwr_min = std::max(0.0, gene.power - 5.0 * amplitude);
                double pwr_max = std::min(5.0, gene.power + 5.0 * amplitude);
                gene.power = fast_rand(pwr_min, pwr_max);

                double dur_min = std::max(0.0, gene.duration - 20.0 * amplitude);
                double dur_max = std::min(20.0, gene.duration + 20.0 * amplitude);
                gene.duration = fast_rand(dur_min, dur_max);
            }
        }
    }

    void AddSurfacePoint(const Point &p)
    {
        ground.pts.push_back(p);
    }

    void CalculateLandingZone()
    {
        for (int i = 0; i < ground.pts.size() - 1; ++i)
        {
            if (ground.pts.at(i).y == ground.pts.at(i + 1).y)
            {
                landing_elev = ground.pts.at(i).y;
                lz_min = ground.pts.at(i).x;
                lz_max = ground.pts.at(i+1).x;
                return;
            }
        }
    }

    void SortPopulationByFitness(std::vector<GenomeAndCommands> &pop)
    {
        for (auto &p : pop)
        {
            EvaluateGenome(p);
        }

        std::sort(pop.begin(), pop.end(), [](const GenomeAndCommands &lh, const GenomeAndCommands &rh) { return lh.fitness_score > rh.fitness_score; });
    }

    void EvaluateGenome(GenomeAndCommands &ind)
    {
        ind.GenerateCommandsFromGenes();
        Lander lander(init_state, ind.commands, ground);
        lander.ComputeTrajectory(dt);

        ind.lander = lander;
        ind.fitness_score = 0.0;
        if (lander.fly_state == LANDED)
        {
            ind.fitness_score = lander.trajectory.back().fuel * 10000.0;
        }
        else if (lander.fly_state == FLYING)
        {
            double xpos = lander.trajectory.back().lander.pos.x;
            double ypos = lander.trajectory.back().lander.pos.y;
            double xvel = lander.trajectory.back().lander.velocity.dx;
            double yvel = lander.trajectory.back().lander.velocity.dy;
            
            ind.fitness_score += yvel;
            ind.fitness_score += (40.0 - std::abs(xvel));
            
            if (lz_min < xpos && xpos < lz_max) // within the lz gets more points
            {
                double dist_from_td = ypos - landing_elev;
                ind.fitness_score += 50.0 / dist_from_td;
            }

        }
        else if (lander.fly_state == CRASHED)
        {
            double xpos = lander.trajectory.back().lander.pos.x;
            double ypos = lander.trajectory.back().lander.pos.y;
            double xvel = lander.trajectory.back().lander.velocity.dx;
            double yvel = lander.trajectory.back().lander.velocity.dy;
            int angle = lander.trajectory.back().angle;

            
            if (-39.0 <= yvel && yvel < 0.0)
            {
                ind.fitness_score += 200.0;
            }
            else if(yvel < -39.0)
            {
                double t = -39 - yvel;
                ind.fitness_score -= t * t * 1.1;
            }
            else if (yvel >= 0.0)
            {
                ind.fitness_score -= yvel;
            }

            if (std::abs(xvel) > 19.0)
            {
                double t = std::abs(19.0 - std::abs(xvel));
                ind.fitness_score -= t*t;
            }
            else if (std::abs(xvel) < 19.0)
            {
                ind.fitness_score += 100.0;
            }
                        
            if (lz_min < xpos && xpos < lz_max) // within the lz gets more points
            {
                ind.fitness_score += 500.0;

                if (angle == 0)
                {
                    ind.fitness_score += 500.0;
                }
                if (std::abs(angle) < 5)
                {
                    ind.fitness_score += 200.0;
                }
                else
                {
                    ind.fitness_score -= std::abs(angle)*3.0;
                }
            }
            else if (xpos < lz_min)
            {
                ind.fitness_score -= (lz_min - xpos);
            }
            else if (xpos > lz_max)
            {
                ind.fitness_score -= (xpos - lz_max);
            }
        }
    }

    State init_state;
    Line ground;
    double dt;
    double landing_elev;
    double lz_min;
    double lz_max;
    std::vector<ControlCommand> commands;
};

void test()
{
    Autopilot ap;
    State is;
    is.angle = -90.0;
    is.fuel = 8001;
    is.power = 0;
    is.lander.pos.x = 500;
    is.lander.pos.y = 2700;
    is.lander.velocity.dx = 100;
    is.lander.velocity.dy = 0;
    ap.AddSurfacePoint(Point(   0, 1000));
    ap.AddSurfacePoint(Point( 300, 1500));
    ap.AddSurfacePoint(Point( 350, 1400));
    ap.AddSurfacePoint(Point( 500, 2000));
    ap.AddSurfacePoint(Point( 800, 1800));
    ap.AddSurfacePoint(Point(1000, 2500));
    ap.AddSurfacePoint(Point(1200, 2100));
    ap.AddSurfacePoint(Point(1500, 2400));
    ap.AddSurfacePoint(Point(2000, 1000));
    ap.AddSurfacePoint(Point(2200,  500));
    ap.AddSurfacePoint(Point(2500,  100));
    ap.AddSurfacePoint(Point(2900,  800));
    ap.AddSurfacePoint(Point(3000,  500));
    ap.AddSurfacePoint(Point(3200, 1000));
    ap.AddSurfacePoint(Point(3500, 2000));
    ap.AddSurfacePoint(Point(3800,  800));
    ap.AddSurfacePoint(Point(4000,  200));
    ap.AddSurfacePoint(Point(5000,  200));
    ap.AddSurfacePoint(Point(5500, 1500));
    ap.AddSurfacePoint(Point(6999, 2800));

    ap.Init(is);
    ap.FindSolution();
    

}

void execute()
{
    Autopilot ap;
    int surfaceN; // the number of points used to draw the surface of Mars.
    std::cin >> surfaceN; std::cin.ignore();
    for (int i = 0; i < surfaceN; i++) {
        int landX; // X coordinate of a surface point. (0 to 6999)
        int landY; // Y coordinate of a surface point. By linking all the points together in a sequential fashion, you form the surface of Mars.
        std::cin >> landX >> landY; std::cin.ignore();
        
        ap.AddSurfacePoint(Point(landX, landY));
    }
    
    bool first_pass = true;
    int step = 0;
    while (1) {
    
        int x;
        int y;
        int hspeed;
        int vspeed;
        int fuel;
        int rotate;
        int power;
        std::cin >> x >> y >> hspeed >> vspeed >> fuel >> rotate >> power; std::cin.ignore();
        if (first_pass)
        {
            State is;
            is.angle = rotate;
            is.fuel = fuel;
            is.power = power;
            is.lander.pos.x = x;
            is.lander.pos.y = y;
            is.lander.velocity.dx = hspeed;
            is.lander.velocity.dy = vspeed;
    
            ap.Init(is);
            ap.FindSolution();
            first_pass = false;
        }
        
        ap.PrintCommand(step);
        ++step;
    }
}
int main()
{
    //test();
    
    execute();

    return 0;
}



