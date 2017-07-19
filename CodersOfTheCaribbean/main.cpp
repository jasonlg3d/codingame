#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <array>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <functional>

// The game is played on a hexagonal grid 23 cells wide and 21 high.
// 
// Both players have one ship (and up to 3 in later leagues). Every ship is 3 cells long and 1 cell wide.
// 
// On the grid you will find barrels of rum (BARREL). The barrels contain between 10 and 20 units of rum. 
// Each ship can carry a maximum of 100 units of rum (surplus will be lost).
// 
// Each turn, the players can decide to move their ship towards the cell of their choice using the MOVE command. 
// Ships can gain speed, slow down, turn left (port) or right (starboard). The MOVE action uses a very simplified 
// algorithm to reach the destination.
// 
// Ships can place mines on the grid with the MINE command. This action spawns a mine in the cell directly behind 
// the ship. After placing a mine, a ship cannot place another for the next 4 turns.
// 
// If a mine is touched by a ship or hit by a cannon ball, it explodes. A ship touching a mine loses 25 units of 
// rum. A ship in a cell adjacent to an exploding mine loses 10 units of rum. Players can only see mines if they 
// are at a distance less than 5 cells away from the center of one of their ships.
// 
// A ship can fire a cannon ball with the FIRE x y command where x and y are the grid coordinates of the target 
// cell. The target must be within 10 cells of the front of the ship. The cannon ball is launched from the front 
// of the ship and will take 1 + (distance to target) / 3 turns to reach the target (the result is rounded). 
// If the cannon ball lands on the front or back of a ship, that ship loses 25 units. The ship loses 50 units if 
// the cannon ball lands on its center. After shooting a cannon ball, a ship cannot fire during the next turn.
// 
// 
// Note: when using a command such as MINE, FIRE and WAIT, the ship will still be moving with the same direction
// and speed as last turn.
// 
// Game turns:
// 
// One game turn is computed as follows:
// The amount of rum each ship is carrying is decreased by 1 unit.
// The players' commands are applied (spawning of cannon balls, mines and ship movement).
// Ships move forward by the same number of cells as their speed.
// Ships turn.
// Damage from cannon balls is computed.
// Elimination of ships with no more rum.
// If at any point during its movement a ship shares a cell with a barrel of rum, the ship collects that rum. In 
// the same way, if a ship touches a mine, that mine explodes and the loss of rum is applied immediately.

static unsigned int g_seed;
const size_t TURNS_TO_CREATE = 176;

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
        if (result < 16384.0)
        {
            return static_cast<double>(result / 16384.0) * min;
        }
        else
        {
            return static_cast<double>(result - 16384.0) / 16384.0 * max;
        }
    }
    else
    {
        return min + static_cast<double>(result / 32767.0) * (max - min);
    }
}


enum class EntityType : int8_t
{
    SHIP, BARREL, MINE, NONE
};

EntityType StringToEntityType(const std::string &str)
{
    if (str.compare("SHIP") == 0)
    {
        return EntityType::SHIP;
    }
    else if (str.compare("BARREL") == 0)
    {
        return EntityType::BARREL;
    }
    else if (str.compare("MINE") == 0)
    {
        return EntityType::MINE;
    }
    else
    {
        return EntityType::NONE;
    }
}

struct Cube
{
    int8_t x;
    int8_t y;
    int8_t z;

    Cube()
        : x(-1)
        , y(-1)
        , z(-1)
    {}

    Cube(int8_t _x, int8_t _y, int8_t _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {}

    bool operator==(const Cube& c)
    {
        return c.x == this->x && c.y == this->y && c.z == this->z;
    }

    bool operator!=(const Cube &c)
    {
        return !(*this == c);
    }

    bool operator<(const Cube &c)
    {
        return false;
    }

    bool operator>(const Cube &c)
    {
        return false;
    }
};
bool operator<(const Cube &lhs, const Cube &rhs)
{
    return false;
}
bool operator==(const Cube &lhs, const Cube &rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool operator!=(const Cube &lhs, const Cube &rhs)
{
    return !(lhs == rhs);
}
std::ostream & operator<<(std::ostream &os, const Cube &c) {
    return os << "(" << static_cast<int>(c.x) << "," << static_cast<int>(c.y) << "," << static_cast<int>(c.z) << ")";
}

struct CubeHash
{
    std::size_t operator()(Cube const &c) const
    {
        std::size_t h1 = std::hash<int>{}(c.x);
        std::size_t h2 = std::hash<int>{}(c.z);
        return h1 ^ (h2 << 1);
    }
};

namespace std
{
template<> struct hash<Cube>
{
    typedef Cube argument_type;
    typedef std::size_t result_type;
    result_type operator() (argument_type const &c) const
    {
        result_type const h1(std::hash<int>{}(c.x));
        result_type const h2(std::hash<int>{}(c.z));
        return h1 ^ (h2 << 1);
    }
};
}

Cube cube_add(const Cube &a, const Cube &b)
{
    return Cube(a.x + b.x, a.y + b.y, a.z + b.z);
}

Cube cube_sub(const Cube &a, const Cube &b)
{
    return Cube(a.x - b.x, a.y - b.y, a.z - b.z);
}

struct Coord
{
    int8_t col;
    int8_t row;

    Coord()
        : row(-1)
        , col(-1)
    {}

    Coord(int8_t _col, int8_t _row)
        : col(_col)
        , row(_row)
    {}

    bool IsValid() const
    {
        return col >= 0 && row >= 0 && col < 23 && row < 21;
    }
};
std::ostream & operator<<(std::ostream &os, const Coord &c) {
    return os << "(" << static_cast<int>(c.col) << "," << static_cast<int>(c.row) << ")";
}

Coord cube_to_oddr(const Cube &cube)
{
    int8_t col = cube.x + (cube.z - (cube.z&1)) / 2;
    int8_t row = cube.z;
    return Coord(col, row);
}

Cube oddr_to_cube(const Coord &coord)
{
    int8_t x = coord.col - (coord.row - (coord.row & 1)) / 2;
    int8_t z = coord.row;
    int8_t y = -x-z;

    return Cube(x,y,z);
}

bool cube_is_valid(const Cube &cube)
{
    return cube_to_oddr(cube).IsValid();
}

enum class Direction : int8_t
{
    W, NW, NE, E, SE, SW, NONE
};
std::ostream & operator<<(std::ostream &os, Direction o) {
    switch (o)
    {
    case Direction::W:
        return os << "W";
    case Direction::NW:
        return os << "NW";
    case Direction::NE:
        return os << "NE";
    case Direction::E:
        return os << "E";
    case Direction::SE:
        return os << "SE";
    case Direction::SW:
        return os << "SW";
    default:
        return os << "NONE";
    }
}

int cube_distance(const Cube &start, const Cube &end)
{
    return std::max(std::max(abs(start.x - end.x), abs(start.y - end.y)), abs(start.z - end.z));
}

std::vector<Cube> in_range(const Cube &center, int range)
{
    std::vector<Cube> result;
    for (int dx = -range; dx <= range; ++dx)
    {
        for (int dy = -range; dy <= range; ++dy)
        {
            int dz = -dx - dy;
            result.push_back(cube_add(center, Cube(dx, dy, dz)));
        }
    }

    return result;
}

struct Barrel
{
    Barrel(int x, int y, int z, int qty_)
        : loc(x, y, z)
        , qty(qty_)
    {}

    Cube loc;
    int qty;
};

bool operator == (const Barrel &lhs, const Barrel &rhs)
{
    return lhs.loc.x == rhs.loc.x && lhs.loc.y == rhs.loc.y && lhs.loc.z == rhs.loc.z;
}
struct BarrelHash
{
    std::size_t operator()(Barrel const &c) const
    {
        std::size_t h1 = std::hash<int>{}(c.loc.x);
        std::size_t h2 = std::hash<int>{}(c.loc.z);
        return h1 ^ (h2 << 1);
    }
};

namespace std
{
template<> struct hash<Barrel>
{
    typedef Barrel argument_type;
    typedef std::size_t result_type;
    result_type operator() (argument_type const &c) const
    {
        result_type const h1(std::hash<int>{}(c.loc.x));
        result_type const h2(std::hash<int>{}(c.loc.z));
        return h1 ^ (h2 << 1);
    }
};
}


struct PriorityQueue
{
    typedef std::pair<double, Cube> PQElement;
    std::priority_queue<PQElement, std::vector<PQElement>, std::greater<PQElement>> elements;

    inline bool empty() const { return elements.empty(); }
    inline void put(const Cube &c, double priority)
    {
        elements.emplace(priority, c);
    }

    Cube get()
    {
        Cube best_item = elements.top().second;
        elements.pop();
        return best_item;
    }
};


struct Ship
{
    Cube loc;
    uint8_t rum_stock;
    int8_t speed;
    Direction orientation;
    int8_t id;
};
std::ostream & operator<<(std::ostream &os, Ship ship) {
    return os << "SHIP " << ship.id << ": " << " @ (" << ship.loc << std::endl
              << "Rum: " << ship.rum_stock  << std::endl
              << "Traveling: " << ship.orientation << "@" << ship.speed << std::endl;
}

class MapGraph
{
public:
    MapGraph(int width_, int height_)
        : m_width(width_)
        , m_height(height_)
    {
        m_cube_directions[static_cast<size_t>(Direction::W)] = Cube(-1, 1, 0);
        m_cube_directions[static_cast<size_t>(Direction::NW)] = Cube(0, 1, -1);
        m_cube_directions[static_cast<size_t>(Direction::NE)] = Cube(1, 0, -1);
        m_cube_directions[static_cast<size_t>(Direction::E)] = Cube(1, -1, 0);
        m_cube_directions[static_cast<size_t>(Direction::SE)] = Cube(0, -1, 1);
        m_cube_directions[static_cast<size_t>(Direction::SW)] = Cube(-1, 0, 1);
    }

    void Clear()
    {
        m_mines.clear();
        m_barrels.clear();
        m_enemy_ships.at(0) = oddr_to_cube(Coord(-1, -1));
        m_enemy_ships.at(1) = oddr_to_cube(Coord(-1, -1));
        m_enemy_ships.at(2) = oddr_to_cube(Coord(-1, -1));
    }

    void AddMyShipLocation(int id, int col, int row)
    {
        m_my_ships[id] = oddr_to_cube(Coord(col, row));
    }

    void AddEnemyShipLocation(int id, int col, int row)
    {
        m_enemy_ships[id] = oddr_to_cube(Coord(col, row));
    }

    void AddBarrelLocation(int col, int row, int qty)
    {
        Cube loc = oddr_to_cube(Coord(col, row));
        Barrel barrel(loc.x, loc.y, loc.z, qty);
        m_barrels.emplace(barrel);
    }

    void AddMineLocation(int col, int row)
    {
        m_mines.emplace(oddr_to_cube(Coord(col, row)));
    }

    Cube CubeInDirection(Direction dir) const
    {
        return m_cube_directions[static_cast<size_t>(dir)];
    }

    Cube Neighbor(const Cube &c, Direction dir) const
    {
        return cube_add(c, CubeInDirection(dir));
    }

    std::vector<Cube> Neighbors(const Cube &loc) const
    {
        std::vector<Cube> results;

        for (int dir = 0; dir < static_cast<int>(Direction::NONE); ++dir)
        {
            auto neighbor = Neighbor(loc, static_cast<Direction>(dir));
            results.push_back(neighbor);
        }

        return results;
    }

    Cube EstimateEnemyShipLocation(const Ship &ship, int turn_count)
    {
        Cube result = ship.loc;
        Coord resultc = cube_to_oddr(result);

        // Have to go two turns to get future state.
        for (int turn = 0; turn < turn_count; ++turn)
        {
            for (int i = 1; i <= ship.speed; ++i)
            {
                result = cube_add(result, CubeInDirection(ship.orientation));
                resultc = cube_to_oddr(result);
            }
        }
        return result;
    }

    void GetPath(const Cube &start, const Cube &end, std::vector<Cube> &path, int ship_id) const
    {
        std::unordered_map<Cube, Cube> came_from;
        std::unordered_map<Cube, double> cost_so_far;
        a_star_search(start, end, came_from, cost_so_far, ship_id);
        auto full_path = reconstruct_path(start, end, came_from);

        for (auto &c : full_path)
        {
            if(c == start) { continue; }
            std::vector<Cube> neighbors = Neighbors(c);
            for (auto &n : neighbors)
            {
                if (m_mines.count(n))
                {
                    path.push_back(c);
                }
            }
        }

        auto found_goal = std::find(path.begin(), path.end(), end);
        if (found_goal == path.end())
        {
            path.push_back(end);
        }
    }

    Direction GetDirection(const Cube &start, const Cube &end) const
    {
        int shortest_dist = 50;
        int dir_count = 0;
        Direction dir = Direction::NONE;
        for (auto &pt : Neighbors(start))
        {
            if (cube_distance(pt, end) < shortest_dist)
            {
                shortest_dist = cube_distance(pt, end);
                dir = static_cast<Direction>(dir_count);
            }
            ++dir_count;
        }
        return dir;
    }

    Cube GetBestShipTarget(const Cube &start) const
    {
        typedef std::pair<double, Cube> PQElement;
        std::priority_queue<PQElement, std::vector<PQElement>, std::less<PQElement>> pq;
        for (auto &s : m_enemy_ships)
        {
            int distance = cube_distance(start, s);
            pq.emplace(distance, s);
        }

        if (!pq.empty())
        {
            std::cerr << "Target ship at " << cube_to_oddr(pq.top().second) << std::endl;
            return pq.top().second;
        }
        else
        {
            std::cerr << "No ships found." << std::endl;
            return oddr_to_cube(Coord(-1, -1));
        }
    }

    Cube GetBestBarrelTarget(const Cube &start, Direction curr_dir, std::vector<Cube> &do_not_target_list) const
    {
        typedef std::pair<double, Cube> PQElement;
        std::priority_queue<PQElement, std::vector<PQElement>, std::less<PQElement>> pq;
        for (auto &b : m_barrels)
        {
            double distance_fact = 170.0;
            double qty_fact      = 1.0;
            int distance = std::max(1, cube_distance(start, b.loc));

            Direction dir_to_barrel = GetDirection(start, b.loc);

            double direction_bonus = 0.0;
            if (dir_to_barrel == curr_dir)
            {
                direction_bonus = 150.0;
            }
            else if ((curr_dir == Direction::W  && (dir_to_barrel == Direction::SW || dir_to_barrel == Direction::NW)) ||
                     (curr_dir == Direction::NW && (dir_to_barrel == Direction::NE || dir_to_barrel == Direction::W )) ||
                     (curr_dir == Direction::NE && (dir_to_barrel == Direction::NW || dir_to_barrel == Direction::E )) ||
                     (curr_dir == Direction::E  && (dir_to_barrel == Direction::NE || dir_to_barrel == Direction::SE)) ||
                     (curr_dir == Direction::SE && (dir_to_barrel == Direction::E  || dir_to_barrel == Direction::SW)) ||
                     (curr_dir == Direction::SW && (dir_to_barrel == Direction::SE || dir_to_barrel == Direction::W ))
                    )
            {
                direction_bonus = 50.0;
            }
            
            auto dnt_found = std::find(do_not_target_list.begin(), do_not_target_list.end(), b.loc);

            auto inrange = in_range(b.loc, 3);

            double qty_total = 0.0;
            for (auto & c : inrange)
            {
                auto found_b = std::find_if(m_barrels.begin(), m_barrels.end(), [&](const Barrel &b) { return b.loc == c; });
                if (found_b != m_barrels.end())
                {
                    qty_total += found_b->qty * 0.5;
                }
            }

            if (dnt_found == do_not_target_list.end())
            {
                qty_total += b.qty;
                pq.emplace(1.0 / distance * distance_fact + qty_total * qty_fact + direction_bonus, b.loc);
            }
        }

        if (!pq.empty())
        {
            std::cerr << "Targeting barrel at " << cube_to_oddr(pq.top().second) << std::endl;
            return pq.top().second;
        }
        else
        {
            std::cerr << "No barrel target found." << std::endl;
            return oddr_to_cube(Coord(-1, -1));
        }
    }

    bool PathObstructed(std::vector<Cube> path, int offset, int ship_id) const
    {
        for (auto it = path.begin() + offset; it < path.end(); ++it)
        {
            if (!Passable(*it, ship_id))
            {
                return true;
            }
        }

        return false;
    }

private:

    inline bool InBounds(const Cube &id) const
    {
        return 0 <= id.x && id.x < m_width && 0 <= id.y && id.y < m_height;
    }

    inline bool Passable(const Cube &id, int ship_id) const
    {
        if (!m_mines.count(id))
        {
            if(id == m_my_ships.at(ship_id))
            { 
                return true; 
            }
            else
            {
                bool ship_blocking = false;

                for (int i = 0; i < 3; ++i)
                {
                    for (auto & c : Neighbors(m_my_ships.at(i)))
                    {
                        if (id == c && ship_id != i)
                        {
                            return false;
                        }
                    }
                    for (auto & c : Neighbors(m_enemy_ships.at(i)))
                    {
                        if (id == c)
                        {
                            return false;
                        }
                    }
                }
            }
        }
        else
        {
            return false;
        }
    }

    double Cost(const Cube &from, const Cube &to) const
    {
        return m_mines.count(to) ? 20.0 : 1.0;
    }

    inline double heuristic(const Cube &a, const Cube &b) const
    {
        return cube_distance(a, b);
    }

    void a_star_search(const Cube &start, const Cube &goal,
        std::unordered_map<Cube, Cube> &came_from,
        std::unordered_map<Cube, double> &cost_so_far, int ship_id) const
    {
        PriorityQueue frontier;
        frontier.put(start, 0.0);

        came_from[start] = start;
        cost_so_far[start] = 0.0;

        while (!frontier.empty())
        {
            auto current = frontier.get();

            if (current == goal)
            {
                break;
            }

            for (auto &next : Neighbors(current))
            {
                double new_cost = cost_so_far[current] + Cost(current, next);

                if (!cost_so_far.count(next) || new_cost < cost_so_far[next] && Passable(next, ship_id))
                {
                    cost_so_far[next] = new_cost;
                    double priority = new_cost + heuristic(next, goal);
                    came_from[next] = current;
                    frontier.put(next, priority);
                }
            }
        }
    }

    void dijkstra_search(const Cube &start, const Cube &goal,
        std::unordered_map<Cube, Cube> &came_from,
        std::unordered_map<Cube, double> &cost_so_far) const
    {
        PriorityQueue frontier;
        frontier.put(start, 0.0);

        came_from[start] = start;
        cost_so_far[start] = 0.0;

        while (!frontier.empty())
        {
            auto current = frontier.get();

            if (current == goal)
            {
                break;
            }

            for (auto &next : Neighbors(current))
            {
                double new_cost = cost_so_far[current] + Cost(current, next);

                if (!cost_so_far.count(next) || new_cost < cost_so_far[next])
                {
                    cost_so_far[next] = new_cost;
                    came_from[next] = current;
                    frontier.put(next, new_cost);
                }
            }
        }
    }

    std::unordered_map<Cube, Cube> breadth_first_search(const Cube &start, const Cube &goal) const
    {
        std::queue<Cube> frontier;
        frontier.push(start);

        std::unordered_map<Cube, Cube> came_from;
        came_from[start] = start;

        while (!frontier.empty())
        {
            auto current = frontier.front();
            frontier.pop();

            if (current == goal)
            {
                break;
            }

            for (auto &next : Neighbors(current))
            {
                if (!came_from.count(next))
                {
                    frontier.push(next);
                    came_from[next] = current;
                }
            }
        }

        return came_from;
    }

    std::vector<Cube> reconstruct_path(const Cube &start, const Cube &goal, std::unordered_map<Cube, Cube> &came_from) const
    {
        std::vector<Cube> path;
        Cube current = goal;
        path.push_back(current);
        while (current != start)
        {
            current = came_from[current];
            path.push_back(current);
        }

        std::reverse(path.begin(), path.end());
        return path;
    }
    int m_width, m_height;
    std::unordered_set<Cube> m_mines;    
    std::unordered_set<Barrel> m_barrels;
    std::array<Cube, 3> m_enemy_ships;
    std::array<Cube, 3> m_my_ships;
    std::array<Cube, 6> m_cube_directions;
    
};

enum MoveType
{
    MOVE, SLOWER, WAIT, FIRE, NONE
};

std::ostream & operator<<(std::ostream &os, MoveType mt) {
    switch (mt)
    {
    case MOVE:
        return os << "MOVE";
    case SLOWER:
        return os << "SLOWER";
    case WAIT:
        return os << "WAIT";
    case FIRE:
        return os << "FIRE";
    default:
        return os << "NONE";
    }
}

std::ostream & operator<< (std::ostream &os, const std::vector<Cube> &cube_list)
{
    for (auto &cube : cube_list)
    {
        os << cube_to_oddr(cube) << std::endl;
    }

    return os;
}
class AI
{
public:
    AI()
        : m_map(23, 21)
        , m_curr_path()
        , m_path_offset()
        , m_my_ships()
        , m_enemy_ships()
        , m_barrels_gone(false)
    {
        for (int idx = 0; idx < 3; ++idx) {
            m_path_offset[idx] = 0;
        }
    }

    void Process()
    {
        Cube fire_at_barrel_tgt[2];
        std::vector<Cube> dnt;
        fire_at_barrel_tgt[0] = m_map.GetBestBarrelTarget(m_enemy_ships.at(0).loc, m_enemy_ships.at(0).orientation, dnt);//oddr_to_cube(Coord(-1, -1));
        fire_at_barrel_tgt[1] = m_map.GetBestBarrelTarget(m_enemy_ships.at(1).loc, m_enemy_ships.at(1).orientation, dnt);//oddr_to_cube(Coord(-1, -1));

        for (int idx = 0; idx < m_my_ship_count; ++idx)
        {
            m_curr_path.at(idx).clear();
            std::vector<Cube> do_not_target_list;

            for (int jdx = 0; jdx < m_my_ship_count; ++jdx)
            {
                if (idx != jdx && !m_curr_path.at(jdx).empty())
                {
                    do_not_target_list.push_back(m_curr_path.at(jdx).back());
                }
            }
            Cube bt = m_map.GetBestBarrelTarget(m_my_ships.at(idx).loc, m_my_ships.at(idx).orientation, do_not_target_list);

            
            m_curr_path.at(idx).clear();
            if (cube_to_oddr(bt).IsValid())
            {
                m_curr_path.at(idx).clear();
                m_map.GetPath(m_my_ships.at(idx).loc, bt, m_curr_path.at(idx), idx);
                m_path_offset.at(idx) = 0;
                std::cerr << "Ship[" << idx << "] Found Path:" << std::endl << m_curr_path.at(idx);
                m_barrels_gone = false;
            }
            else
            {
                m_barrels_gone = true;
                Cube st = m_map.GetBestShipTarget(m_my_ships.at(idx).loc);

                if (cube_is_valid(st))
                {
                    m_map.GetPath(m_my_ships.at(idx).loc, st, m_curr_path.at(idx), idx);
                    m_path_offset.at(idx) = 0;
                    std::cerr << "Ship[" << idx << "] Found Path:" << std::endl << m_curr_path.at(idx);
                }
                else
                {
                    m_map.GetPath(m_my_ships.at(idx).loc, oddr_to_cube(Coord(fast_rand(0, 22), fast_rand(0, 20))), m_curr_path.at(idx), idx);
                    m_path_offset.at(idx) = 0;
                    std::cerr << "Ship[" << idx << "] Found Path:" << std::endl << m_curr_path.at(idx);
                }
            }
        }

        for (int idx = 0; idx < m_my_ship_count; ++idx)
        {
            if (!m_curr_path.at(idx).empty() && m_curr_path.at(idx).at(m_path_offset.at(idx)) == m_my_ships.at(idx).loc)
            {
                m_path_offset.at(idx)++;
            }

            Cube fire_tgt;            

            auto inrg_list = in_range(m_my_ships.at(idx).loc, 10);
            auto found_tgt = std::find(inrg_list.begin(), inrg_list.end(), fire_at_barrel_tgt[0]);
            if (found_tgt == inrg_list.end())
            {
                found_tgt = std::find(inrg_list.begin(), inrg_list.end(), fire_at_barrel_tgt[1]);
            }

            
            if (CalculateFiringSolution(idx, fire_tgt) && (OnCourse(idx) || m_barrels_gone) && m_fire_tgt_prev != fire_tgt)
            {
                Coord fire_coord = cube_to_oddr(fire_tgt);
                std::cout << "FIRE " << static_cast<int>(fire_coord.col) << " " << static_cast<int>(fire_coord.row) << std::endl;
            }
            else if (found_tgt != inrg_list.end() && m_fire_tgt_prev != *found_tgt && m_my_ships.at(idx).rum_stock > 50)
            {
                fire_tgt = *found_tgt;
                Coord fire_coord = cube_to_oddr(fire_tgt);
                std::cout << "FIRE " << static_cast<int>(fire_coord.col) << " " << static_cast<int>(fire_coord.row) << std::endl;
            }
            else if (m_path_offset.at(idx) < m_curr_path.at(idx).size())
            {
                Coord c = cube_to_oddr(m_curr_path.at(idx).at(m_path_offset.at(idx)));
                std::cout << "MOVE " << static_cast<int>(c.col) << " " << static_cast<int>(c.row) << std::endl;                
            }
            else
            {
                std::cout << "MINE" << std::endl;
            }

            m_fire_tgt_prev = fire_tgt;
        }
    }

    void Reset()
    {
        m_my_ship_array_offset = 0;
        m_enemy_ship_array_offset = 0;
        m_map.Clear();
    }

    void SetShipCount(int sc)
    {
        m_my_ship_count = sc;
    }

    void AddEntityToCurrentState(int id, const std::string &type, int x, int y, int arg1, int arg2, int arg3, int arg4)
    {
        EntityType etype = StringToEntityType(type);
        switch (etype)
        {
        case EntityType::BARREL:
            m_map.AddBarrelLocation(x, y, arg1);
            break;
        case EntityType::MINE:
            m_map.AddMineLocation(x, y);
            break;
        case EntityType::SHIP:
            if (arg4 == 1)
            {
                m_map.AddMyShipLocation(m_my_ship_array_offset, x, y);
                m_my_ships.at(m_my_ship_array_offset).loc = oddr_to_cube(Coord(x, y));
                m_my_ships.at(m_my_ship_array_offset).rum_stock = arg3;
                m_my_ships.at(m_my_ship_array_offset).orientation = static_cast<Direction>(arg1);
                m_my_ships.at(m_my_ship_array_offset).speed = arg2;
                m_my_ships.at(m_my_ship_array_offset).id = id;
                ++m_my_ship_array_offset;
            }
            else
            {
                m_map.AddEnemyShipLocation(m_enemy_ship_array_offset, x, y);
                m_enemy_ships.at(m_enemy_ship_array_offset).loc = oddr_to_cube(Coord(x, y));
                m_enemy_ships.at(m_enemy_ship_array_offset).rum_stock = arg3;
                m_enemy_ships.at(m_enemy_ship_array_offset).orientation = static_cast<Direction>(arg1);
                m_enemy_ships.at(m_enemy_ship_array_offset).speed = arg2;
                m_enemy_ships.at(m_enemy_ship_array_offset).id = id;
                ++m_enemy_ship_array_offset;
            }
            break;
        case EntityType::NONE:
        default:
            break;
        }
    }
private:
    bool OnCourse(int ship_id)
    {
        Ship &ship = m_my_ships.at(ship_id);
        Cube dest = m_curr_path.at(ship_id).at(0);

        Cube next_turn_loc = ship.loc;
        for (int spd = 1; spd <= ship.speed; ++spd)
        {
            next_turn_loc = cube_add(next_turn_loc, m_map.CubeInDirection(ship.orientation));
        }

        return cube_distance(next_turn_loc, dest) < cube_distance(ship.loc, dest);
    }
    
    bool CalculateFiringSolution(int my_ship_id, Cube &fire_tgt)
    {
        // A ship can fire a cannon ball with the FIRE x y command where x and y are the grid coordinates of the target cell. 
        // The target must be within 10 cells of the front of the ship. The cannon ball is launched from the front of the 
        // ship and will take 1 + (distance to target) / 3 turns to reach the target (the result is rounded). If the cannon 
        // ball lands on the front or back of a ship, that ship loses 25 units. The ship loses 50 units if the cannon ball 
        // lands on its center. After shooting a cannon ball, a ship cannot fire during the next turn.
        
        // Ships move forward by the same number of cells as their speed.
        Ship &my_ship = m_my_ships.at(my_ship_id);

        Cube target;
        for (auto &eship : m_enemy_ships)
        {
            for (int turn = 1; turn < 5; ++turn)
            {
                target = m_map.EstimateEnemyShipLocation(eship, turn);

                int dist = cube_distance(my_ship.loc, target);

                if (static_cast<int>(std::floor((1.0 + static_cast<double>(dist) / 3.0) + 0.5)) == turn)
                {
                    fire_tgt = target;
                    Coord fire_c = cube_to_oddr(fire_tgt);
                    std::cerr << "Fire solution: " << fire_c << " Dist: " << dist << std::endl;
                    return true;
                }
            }
        }    

        return false;
    }

    MapGraph m_map;
    std::array<std::vector<Cube>, 3> m_curr_path;
    std::array<int, 3> m_path_offset;
    std::array<Ship, 3> m_my_ships;
    std::array<Ship, 3> m_enemy_ships;
    int m_my_ship_count;
    int m_my_ship_array_offset;
    int m_enemy_ship_array_offset;
    bool m_barrels_gone;
    Cube m_fire_tgt_prev;
};


void test()
{
    AI ai;
    ai.SetShipCount(1);;
    ai.AddEntityToCurrentState(0, "SHIP", 10, 2, static_cast<int>(Direction::E), 0, 100, 1);

    ai.AddEntityToCurrentState(1, "SHIP", 9, 5, static_cast<int>(Direction::W), 2, 100, 0);

    ai.AddEntityToCurrentState(0, "MINE", 15, 2, 0, 0, 0, 0);

    ai.AddEntityToCurrentState(0, "BARREL", 20, 2, 20, 0, 0, 0);

    ai.Process();
    ai.Process();
}
int main()
{
    //test();
    AI ai;
    // game loop
    while (1) {

        ai.Reset();
        int myShipCount; // the number of remaining ships
        std::cin >> myShipCount; std::cin.ignore();
        ai.SetShipCount(myShipCount);

        int entityCount; // the number of entities (e.g. ships, mines or cannonballs)
        std::cin >> entityCount; std::cin.ignore();
        
        for (int i = 0; i < entityCount; i++) {
            int entityId;
            std::string entityType;
            int x;
            int y;
            int arg1;
            int arg2;
            int arg3;
            int arg4;
            std::cin >> entityId >> entityType >> x >> y >> arg1 >> arg2 >> arg3 >> arg4; std::cin.ignore();
            //std::cerr << entityId << " " << entityType << " " << x << " " << y << " " << arg1 << " " << arg2 << " " << arg3 << " " << arg4 << std::endl;
            ai.AddEntityToCurrentState(entityId, entityType, x, y, arg1, arg2, arg3, arg4);
        }

        ai.Process();
    }
}