#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <array>
#include <cstring>
#include <memory>
#include <chrono>

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

enum Direction : int8_t {
    N,   // 0
    NE,  // 1
    E,   // 2
    SE,  // 3
    S,   // 4
    SW,  // 5
    W,   // 6
    NW,  // 7
    NONE // 8
};

bool PushIsValid(Direction from, Direction to) {
    return std::abs(to-from) <= 1 ||
        from == Direction::N && to == Direction::NW ||
        from == Direction::NW && to == Direction::N;
}

enum Type : int8_t {
    MOVE_BUILD,
    PUSH_BUILD
};

Type StringToType(const std::string& type_str)
{
    if (type_str.compare("MOVE&BUILD") == 0)
    {
        return Type::MOVE_BUILD;
    }
    else if(type_str.compare("PUSH&BUILD") == 0)
    {
        return Type::PUSH_BUILD;
    }
    else
    {
        return Type::MOVE_BUILD;
    }
}
struct Point
{
    int8_t x, y;

    Point()
        : x(-1)
        , y(-1)
    {}

    Point(int8_t _x, int8_t _y)
        : x(_x)
        , y(_y)
    {}

    Point operator+(int dir) const {
        return (*this) + static_cast<Direction>(dir);
    }

    Point operator+=(int dir) const
    {
        return (*this) + static_cast<Direction>(dir);
    }

    Point operator+(Direction dir) const {
        switch(dir) {
            case Direction::N:
                return Point(x,y-1);
            case Direction::NE:
                return Point(x+1,y-1);
            case Direction::E:
                return Point(x+1,y);
            case Direction::SE:
                return Point(x+1, y+1);
            case Direction::S:
                return Point(x, y+1);
            case Direction::SW:
                return Point(x-1, y+1);
            case Direction::W:
                return Point(x-1, y);
            case Direction::NW:
                return Point(x-1, y-1);
            case Direction::NONE:
            default:
                return Point(x,y);
        }
    }

    bool operator==(const Point&b) const
    {
        return this->x == b.x && this->y == b.y;
    }

    bool operator!=(const Point&b) const
    {
        return !(*this == b);
    }


};

bool AreAdjacent(const Point &pta, const Point &ptb)
{
    return std::abs(pta.x - ptb.x) < 2 && std::abs(pta.y - ptb.y) < 2;
}

std::ostream & operator<<(std::ostream &os, Point pt) {
    return os << "(" << static_cast<int>(pt.x) << "," << static_cast<int>(pt.y) << ")";
}

struct GameState
{
    GameState()
        : player_loc()
        , opp_loc()
        , grid_size(0)
        , blocks()
        , is_valid(true)
    {}

    int Block(const Point & pt) const {
        return blocks[pt.x][pt.y];
    }

    bool Visited(const Point & pt) const {
        return visited[pt.x][pt.y] != 0;
    }

    void Clear()
    {
        for (int i = 0; i < 2; ++i)
        {
            player_loc[i].x = -1;
            player_loc[i].y = -1;
            opp_loc[i].x = -1;
            opp_loc[i].y = -1;
        }

        for (int i = 0; i < 7; ++i)
        {
            for (int j = 0; j < 7; ++j)
            {
                blocks[i][j] = -1;
            }
        }
    }

    bool IsOccupied(const Point & pt, int ignore_player = -1) const {
        return 
            ignore_player != 0 && player_loc[0] == pt ||
            ignore_player != 1 && player_loc[1] == pt ||
            opp_loc[0] == pt ||
            opp_loc[1] == pt;
    }

    bool CanMove(const Point & pt_from, Direction dir, int unit_idx) const {
        int opp_idx = unit_idx == 0 ? 1 : 0;
        Point new_pt = pt_from + dir;
        if(new_pt.x < 0 || new_pt.x >= grid_size) return false;
        if(new_pt.y < 0 || new_pt.y >= grid_size) return false;
        int8_t frm_ht = blocks[pt_from.x][pt_from.y];
        int8_t dest_ht = blocks[new_pt.x][new_pt.y];
        return  dest_ht < 4 && dest_ht >= 0 && (dest_ht - frm_ht) <= 1 && !IsOccupied(new_pt);
    }

    bool CanMove(const Point & pt_from, int dir, int unit_idx) const {
        return CanMove(pt_from, static_cast<Direction>(dir), unit_idx);
    }

    bool CanBuild(const Point & pt_from, Direction dir, int unit_idx) const {
        int opp_idx = unit_idx == 0 ? 1 : 0;
        Point new_pt = pt_from + dir;
        if(new_pt.x < 0 || new_pt.x >= grid_size) return false;
        if(new_pt.y < 0 || new_pt.y >= grid_size) return false;
        int8_t dest_ht = blocks[new_pt.x][new_pt.y];
        return dest_ht >= 0 && dest_ht < 4 && !IsOccupied(new_pt, unit_idx);
    }
    bool CanBuild(const Point & pt_from, int dir, int unit_idx) const {
        return CanBuild(pt_from, static_cast<Direction>(dir), unit_idx);
    }

    
    bool CanPush(const Point & pt_from, int push_from_dir, int push_to_dir) const {
        return CanPush(pt_from, static_cast<Direction>(push_from_dir), static_cast<Direction>(push_to_dir));
    }
    

    bool CanPush(const Point &pt_from, Direction push_from_dir, Direction push_to_dir) const 
    {
        if(!PushIsValid(push_from_dir, push_to_dir)) return false;
        Point push_loc = pt_from + push_from_dir;
        for(int i = 0; i < 2; ++i) {
            if (push_loc == opp_loc[i])
            {
                Point new_pt = opp_loc[i] + push_to_dir;
                if (new_pt.x < 0 || new_pt.x >= grid_size) return false;
                if (new_pt.y < 0 || new_pt.y >= grid_size) return false;
                int8_t frm_ht = blocks[opp_loc[i].x][opp_loc[i].y];
                int8_t dest_ht = blocks[new_pt.x][new_pt.y];
                return  dest_ht < 4 && dest_ht >= 0 && (dest_ht - frm_ht) <= 1 && !IsOccupied(new_pt);
            }
        }
        return false;
    }
    
    double Score(bool use_visitable, int curr_turn_num)
    {
        if (!is_valid)
        {
            return 0;
        }
        else
        {
            double defensive_weight = std::min(2.5, std::max(1.0, 2.0 / curr_turn_num));
            score = ScorePlayer(player_loc[0], use_visitable, curr_turn_num) + 
                    ScorePlayer(player_loc[1], use_visitable, curr_turn_num) +
                    defensive_weight * 
                    (ScoreOpponent(opp_loc[0], use_visitable, curr_turn_num) + 
                     ScoreOpponent(opp_loc[1], use_visitable, curr_turn_num));
            return score;
        }
    }

    Point HighestPointPotential()
    {
        Point best_point;
        double best_score = -10000000.0;
        for (int x = 0; x < 7; ++x)
        {
            for (int y = 0; y < 7; ++y)
            {
                double score = ScorePoint(x, y);
                if (score > best_score)
                {
                    best_score = score;
                    best_point.x = x;
                    best_point.y = y;
                }
            }
        }
        return best_point;
    }

    std::vector<Point> DeathTraps()
    {
        std::vector<Point> traps;
        for (int x = 0; x < 7; ++x)
        {
            for (int y = 0; y < 7; ++y)
            {
                std::vector<Point> adj_pts;
                Point pt_cons(x,y);
                if (Block(pt_cons) < 0) { continue; }
                for (int dir = Direction::N; dir < Direction::NONE; ++dir) {
                    Point adj_pt = pt_cons + dir;
                    if (adj_pt.x < 0 || adj_pt.x >= grid_size) continue;
                    if (adj_pt.y < 0 || adj_pt.y >= grid_size) continue;
                    adj_pts.push_back(adj_pt);                    
                }
                int pt_ht = blocks[x][y];
                bool is_trap = true;
                for (auto &pt : adj_pts)
                {
                    if ((pt_ht >= Block(pt) || pt_ht + 1 == Block(pt)) && Block(pt) > 0)
                    {
                        is_trap = false;
                    }
                }
                if (is_trap)
                {
                    traps.push_back(pt_cons);
                }
            }
        }
        return traps;
    }

    double ScorePoint(int x, int y)
    {
        ScorePlayer(Point(x, y), false, 0);
    }


    double ScoreOpponent(const Point &unit_loc, bool use_visitable, int curr_turn_num)
    {
        double sc = 0.0;

        if(unit_loc == Point(-1,-1)) return sc;
        int values[3][3] = { {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1} };
        for (int dx = -1; dx < 2; ++dx) {
            for (int dy = -1; dy < 2; ++dy) {
                int x = unit_loc.x + dx;
                int y = unit_loc.y + dy;
                if (x < 0 || x >= grid_size) continue;
                if (y < 0 || y >= grid_size) continue;
                if (blocks[x][y] == 4) {
                    values[dx + 1][dy + 1] = -2;
                    continue;
                }
                else
                {
                    values[dx + 1][dy + 1] = blocks[x][y];
                }
            }
        }

        sc -= values[1][1] * 10;

        for (int x = 0; x < 3; ++x)
        {
            for (int y = 0; y < 3; ++y)
            {
                if (values[x][y] - values[1][1] > 1)
                {
                    sc += 10.0;
                }
                else
                {
                    sc -= values[x][y] * values[x][y] * values[x][y];
                }
            }
        }

        if (use_visitable)
        {
            int visitable = CountVisitablePoints(unit_loc);
            if (visitable > 0)
            {
                sc += 50.0 / visitable;

            }

            auto dtvec = DeathTraps();
            for (auto & dt : dtvec)
            {
                if (AreAdjacent(unit_loc, dt))
                {
                    sc += curr_turn_num * 10.0;
                }

                if (unit_loc == dt)
                {
                    sc += 50.0;
                }
            }
        }
        return sc;
    }

    double ScorePlayer(const Point &unit_loc, bool use_visitable, int curr_turn_num)
    {
        if(unit_loc == Point(-1,-1)) return 0;
        int values[3][3] = { {-1, -1, -1}, {-1, -1, -1}, {-1, -1, -1} };
        for (int dx = -1; dx < 2; ++dx) {
            for (int dy = -1; dy < 2; ++dy) {
                int x = unit_loc.x + dx;
                int y = unit_loc.y + dy;
                if (x < 0 || x >= grid_size) continue;
                if (y < 0 || y >= grid_size) continue;
                if (blocks[x][y] == 4) {
                    values[dx + 1][dy + 1] = -1;
                    continue;
                }
                else
                {
                    values[dx + 1][dy + 1] = blocks[x][y];
                }
            }
        }

        double sc = 0;
        double point_fact = 25;
        double potential_point_fact = 10;
        sc += values[1][1] * 10;

        if (values[1][1] == 3)
        {
            sc += point_fact;
        }

        for (int x = 0; x < 3; ++x)
        {
            for (int y = 0; y < 3; ++y)
            {
                if (values[x][y] - values[1][1] > 1)
                {
                    sc += -10;
                }
                else if (values[x][y] == 3 && values[1][1] > 1)
                {
                    sc += potential_point_fact;
                }
                else
                {
                    sc += values[x][y] * values[x][y] * values[x][y];
                }
            }
        }

        if (use_visitable)
        {
            sc += CountVisitablePoints(unit_loc) * curr_turn_num;

            auto dtvec = DeathTraps();
            for (auto & dt : dtvec)
            {
                if (AreAdjacent(unit_loc, dt))
                {
                    sc -= curr_turn_num * 0.0;
                }

                if (unit_loc == dt)
                {
                    sc -= 1000.0;
                }
            }
        }
        return sc;
    }
    
    void Visit(const Point & pt) {
        visited[pt.x][pt.y] = 1;
    }

    void AdjacentPoints(const Point & pt, std::vector<Point> & adj_pts) {
        for(int dir = Direction::N; dir < Direction::NONE; ++dir) {
            Point adj_pt = pt + dir;
            if (Visited(adj_pt)) continue;
            if (adj_pt.x < 0 || adj_pt.x >= grid_size) continue;
            if (adj_pt.y < 0 || adj_pt.y >= grid_size) continue;
            Visit(adj_pt);
            int value = Block(adj_pt);
            if(value >= 0 && value < 4 && !IsOccupied(adj_pt)) {
                adj_pts.push_back(adj_pt);
                AdjacentPoints(adj_pt, adj_pts);
            }
        }
    }

    int CountVisitablePoints(const Point & pt) {
        memset(visited, 0, sizeof(visited));
        Visit(pt);
        std::vector<Point> visitable;
        AdjacentPoints(pt, visitable);
        return static_cast<int>(visitable.size() + 1);
    }

    Point player_loc[2];
    Point opp_loc[2];
    double score;
    int8_t grid_size;
    int8_t blocks[7][7];
    int8_t visited[7][7];
    bool is_valid;
};

std::ostream & operator<<(std::ostream &os, const GameState &gs) {
    return os << "Player Unit 0:   " << gs.player_loc[0] << std::endl
              << "Player Unit 1:   " << gs.player_loc[1] << std::endl
              << "Opponent Unit 0: " << gs.opp_loc[0] << std::endl
              << "Opponent Unit 1: " << gs.opp_loc[1] << std::endl;
              //<< "Grid:" << std::endl
              //<< static_cast<int>(gs.blocks[0][0]) << static_cast<int>(gs.blocks[1][0]) << static_cast<int>(gs.blocks[2][0]) << static_cast<int>(gs.blocks[3][0]) << static_cast<int>(gs.blocks[4][0]) << static_cast<int>(gs.blocks[5][0]) << static_cast<int>(gs.blocks[6][0]) << std::endl
              //<< static_cast<int>(gs.blocks[0][1]) << static_cast<int>(gs.blocks[1][1]) << static_cast<int>(gs.blocks[2][1]) << static_cast<int>(gs.blocks[3][1]) << static_cast<int>(gs.blocks[4][1]) << static_cast<int>(gs.blocks[5][1]) << static_cast<int>(gs.blocks[6][1]) << std::endl
              //<< static_cast<int>(gs.blocks[0][2]) << static_cast<int>(gs.blocks[1][2]) << static_cast<int>(gs.blocks[2][2]) << static_cast<int>(gs.blocks[3][2]) << static_cast<int>(gs.blocks[4][2]) << static_cast<int>(gs.blocks[5][2]) << static_cast<int>(gs.blocks[6][2]) << std::endl
              //<< static_cast<int>(gs.blocks[0][3]) << static_cast<int>(gs.blocks[1][3]) << static_cast<int>(gs.blocks[2][3]) << static_cast<int>(gs.blocks[3][3]) << static_cast<int>(gs.blocks[4][3]) << static_cast<int>(gs.blocks[5][3]) << static_cast<int>(gs.blocks[6][3]) << std::endl
              //<< static_cast<int>(gs.blocks[0][4]) << static_cast<int>(gs.blocks[1][4]) << static_cast<int>(gs.blocks[2][4]) << static_cast<int>(gs.blocks[3][4]) << static_cast<int>(gs.blocks[4][4]) << static_cast<int>(gs.blocks[5][4]) << static_cast<int>(gs.blocks[6][4]) << std::endl
              //<< static_cast<int>(gs.blocks[0][5]) << static_cast<int>(gs.blocks[1][5]) << static_cast<int>(gs.blocks[2][5]) << static_cast<int>(gs.blocks[3][5]) << static_cast<int>(gs.blocks[4][5]) << static_cast<int>(gs.blocks[5][5]) << static_cast<int>(gs.blocks[6][5]) << std::endl
              //<< static_cast<int>(gs.blocks[0][6]) << static_cast<int>(gs.blocks[1][6]) << static_cast<int>(gs.blocks[2][6]) << static_cast<int>(gs.blocks[3][6]) << static_cast<int>(gs.blocks[4][6]) << static_cast<int>(gs.blocks[5][6]) << static_cast<int>(gs.blocks[6][6]) << std::endl;
}

Direction StringToDirection(const std::string &str)
{
    if (str.compare("N") == 0)
    {
        return Direction::N;
    }
    else if (str.compare("NE") == 0)
    {
        return Direction::NE;
    }
    else if (str.compare("E") == 0)
    {
        return Direction::E;
    }
    else if (str.compare("SE") == 0)
    {
        return Direction::SE;
    }
    else if (str.compare("S") == 0)
    {
        return Direction::S;
    }
    else if (str.compare("SW") == 0)
    {
        return Direction::SW;
    }
    else if (str.compare("W") == 0)
    {
        return Direction::W;
    }
    else if (str.compare("NW") == 0)
    {
        return Direction::NW;
    }
    else
    {
        return Direction::NONE;
    }
}

std::string DirectionToString(Direction dir)
{
    switch (dir)
    {
    case Direction::N:
        return "N";
    case Direction::NE:
        return "NE";
    case Direction::E:
        return "E";
    case Direction::SE:
        return "SE";
    case Direction::S:
        return "S";
    case Direction::SW:
        return "SW";
    case Direction::W:
        return "W";
    case Direction::NW:
        return "NW";
    case Direction::NONE:
    default:
        return "";
    }
}

std::ostream & operator<<(std::ostream & os, Direction dir) {
    return os << DirectionToString(dir);
}

std::string TypeToString(int8_t type)
{
    switch (type)
    {
    case Type::MOVE_BUILD:
        return "MOVE&BUILD";
    case Type::PUSH_BUILD:
        return "PUSH&BUILD";
    default:
        return "MOVE&BUILD";
    }
}

std::string TypeToString(Type type)
{
    switch (type)
    {
    case Type::MOVE_BUILD:
        return "MOVE&BUILD";
    case Type::PUSH_BUILD:
        return "PUSH&BUILD";
    default:
        return "MOVE&BUILD";
    }
}
std::ostream & operator<<(std::ostream & os, Type type) {
    return os << TypeToString(type);
}

struct Move
{
    void Clear()
    {
        move_dir  = Direction::NONE;
        build_dir = Direction::NONE;
        type      = Type::MOVE_BUILD;
        bot_idx   = 0;
    }
    int8_t type;
    int8_t move_dir;
    int8_t build_dir;
    int8_t bot_idx;
};

std::ostream & operator<<(std::ostream & os, Move move) {
    return os << TypeToString(move.type) << " " << static_cast<int>(move.bot_idx) << " " << static_cast<Direction>(move.move_dir) << " " << static_cast<Direction>(move.build_dir);
}

class AI
{
public:
    AI(int8_t size)
        : m_grid_size(size)
    {}

    virtual void Clear() = 0;
    virtual void Process() = 0;
    virtual void SetGameStateRow(int row_idx, const std::string &row_data) = 0;
    virtual void SetGameStatePlayer(int player_idx, int x, int y) = 0;
    virtual void SetGameStateOpponent(int opp_idx, int x, int y) = 0;
    virtual void AddGameStateLegalAction(int legal_action_count, int idx, const std::string &type, int unit_idx, const std::string &move_dir, const std::string &build_dir) = 0;
private:
    int8_t m_grid_size;
};

std::vector<Move> GetPossiblePlayerMoves(const GameState &state, int unit_idx)
{
    std::vector<Move> result;
    for (int m_dir = Direction::N; m_dir < Direction::NONE; ++m_dir) {
        if (state.CanMove(state.player_loc[unit_idx], m_dir, unit_idx)) {
            Point new_pt = state.player_loc[unit_idx] + m_dir;
            for (int b_dir = Direction::N; b_dir < Direction::NONE; ++b_dir) {
                if (state.CanBuild(new_pt, b_dir, unit_idx)) {
                    Move move;
                    move.build_dir = b_dir;
                    move.move_dir = m_dir;
                    move.type = Type::MOVE_BUILD;
                    move.bot_idx = unit_idx;
                    result.push_back(move);
                }
            }
        }

        for (int b_dir = Direction::N; b_dir < Direction::NONE; ++b_dir) {
            if (state.CanPush(state.player_loc[unit_idx], m_dir, b_dir))
            {
                Move move;
                move.build_dir = b_dir;
                move.move_dir = m_dir;
                move.type = Type::PUSH_BUILD;
                move.bot_idx = unit_idx;
                result.push_back(move);
            }
        }
    }

    return result;
}

struct Turn
{
    Move player_move;
    Move opponent_move;
    GameState result;
    double score;
    int depth;
    bool is_valid;
    std::vector<Turn*> next_turns;

    Turn(int _depth, int max_depth);

    ~Turn();

    void Clear();

    void GenerateSubTurns(int max_depth, int current_turn_num);

    double AggregateScore() const;

    void Simulate(const GameState &curr_state, bool use_visitable, int curr_turn_num);

    GameState Simulate(const Turn &turn, const GameState &state);

};

Turn::~Turn()
{
    for (auto ptr : next_turns)
    {
        delete ptr;
        ptr = nullptr;
    }

    next_turns.clear();
}

Turn::Turn(int _depth, int max_depth)
    : player_move()
    , opponent_move()
    , result()
    , score(0.0)
    , depth(_depth)
    , is_valid(false)
    , next_turns()
{
    if (depth < max_depth)
    {
        next_turns.reserve(TURNS_TO_CREATE);
        for (int i = 0; i < TURNS_TO_CREATE; ++i)
        {
            next_turns.push_back(new Turn(depth + 1, max_depth));
        }
    }
}

void Turn::Clear()
{
    is_valid = false;
    score = std::numeric_limits<double>::lowest();
    result = GameState();
    player_move = Move();
    opponent_move = Move();
}

void Turn::GenerateSubTurns(int max_depth, int current_turn_num)
{
    if (depth == max_depth) { return; }

    for (auto &t : next_turns)
    {
        t->Clear();
    }
    std::vector<Move> next_poss_moves[2];
    next_poss_moves[0] = GetPossiblePlayerMoves(result, 0);
    next_poss_moves[1] = GetPossiblePlayerMoves(result, 1);

    size_t total_possible_turns = next_poss_moves[0].size() + next_poss_moves[1].size();
    size_t size = next_turns.size();

    //if (size > 0)
    //{
    //    int i = 0;
    //    ++i;
    //}
    //while (next_turns.size() < total_possible_turns)
    //{
    //    next_turns.push_back(new Turn());
    //}

    size_t count = 0;
    for (int i = 0; i < 2; ++i)
    {
        for (auto & m : next_poss_moves[i])
        {
            next_turns.at(count)->player_move = m;
            ++count;
        }
    }

    for (auto &t : this->next_turns)
    {
        if (t->is_valid)
        {
            t->Simulate(this->result, t->depth < 1, current_turn_num);
        }
        else
        {
            break;
        }
    }

    std::sort(this->next_turns.begin(), this->next_turns.end(), [](const Turn *a, const Turn* b)
    {
        return a->result.score > b->result.score;
    });

    size_t number_to_keep = std::max(50 - 2 * depth, 1);
    for (int i = 0; i < std::min(number_to_keep, next_turns.size()); ++i)
    {
        next_turns.at(i)->GenerateSubTurns(max_depth, current_turn_num);
    }
}

double Turn::AggregateScore() const
{
    double s = this->score;        
    if (!next_turns.empty())
    {
        if (depth > 0)
        {
            if (next_turns.at(0)->is_valid)
            {
                s += next_turns.at(0)->AggregateScore() / (2 * depth);
            }
        }
    }
    return s;
}

void Turn::Simulate(const GameState &curr_state, bool use_visitable, int curr_turn_num)
{
    this->result = Simulate(*this, curr_state);
    this->score = this->result.Score(use_visitable, curr_turn_num);
}

GameState Turn::Simulate(const Turn &turn, const GameState &state)
{
    GameState new_state = state;
    for (int i = 0; i < 2; ++i)
    {
        // This check is to confirm the move is for the correct bot and that bot is actually active (-1,-1) indicates no bot.
        if (new_state.player_loc[i].x >= 0 && new_state.player_loc[i].y >= 0 && turn.player_move.bot_idx == i)
        {
            if (turn.player_move.type == Type::MOVE_BUILD)
            {
                if (!new_state.CanMove(new_state.player_loc[i], turn.player_move.move_dir, i))
                {
                    new_state.is_valid = false;
                    return new_state;
                }
                new_state.player_loc[i] = new_state.player_loc[i] + turn.player_move.move_dir;

            }
            else if (turn.player_move.type == Type::PUSH_BUILD)
            {
                if (!new_state.CanPush(new_state.player_loc[i], turn.player_move.move_dir, turn.player_move.build_dir))
                {
                    new_state.is_valid = false;
                    return new_state;
                }
                Point push_loc = new_state.player_loc[i] + turn.player_move.move_dir;
                    
                if (push_loc == new_state.opp_loc[0])
                {
                    new_state.opp_loc[0] = new_state.opp_loc[0] + turn.player_move.build_dir;
                }
                else if (push_loc == new_state.opp_loc[1])
                {
                    new_state.opp_loc[1] = new_state.opp_loc[1] + turn.player_move.build_dir;
                }
            }
        }

        //if (new_state.opp_loc[i].x >= 0 && new_state.opp_loc[i].y >= 0)
        //{
        //    if (!new_state.CanMove(new_state.opp_loc[i], turn.opponent_move[i].move_dir))
        //    {
        //        new_state.is_valid = false;
        //        return new_state;
        //    }
        //    new_state.opp_loc[i] = new_state.opp_loc[i] + turn.opponent_move[i].move_dir;
        //}
    }

    for (int i = 0; i < 2; ++i)
    {
        if (new_state.player_loc[i].x >= 0 && new_state.player_loc[i].y >= 0 && turn.player_move.bot_idx == i)
        {
            if (!new_state.CanBuild(new_state.player_loc[i], turn.player_move.build_dir, i))
            {
                new_state.is_valid = false;
                return new_state;
            }
            Point build_loc = new_state.player_loc[i] += turn.player_move.build_dir;        
            new_state.blocks[build_loc.x][build_loc.y] += 1;
        }

        //if (new_state.opp_loc[i].x >= 0 && new_state.opp_loc[i].y >= 0)
        //{
        //    if (!new_state.CanBuild(new_state.opp_loc[i], turn.opponent_move[i].build_dir))
        //    {
        //        new_state.is_valid = false;
        //        return new_state;
        //    }
        //    Point build_loc = new_state.opp_loc[i] += turn.opponent_move[i].build_dir;
        //    new_state.blocks[build_loc.x][build_loc.y] += 1;
        //}
    }

    return new_state;
}

class GeneticAI : public AI
{
public:
    GeneticAI(int8_t size, int max_turn_depth)
        : AI(size)
        , m_max_turn_depth(max_turn_depth)
        , m_turn_count(0)
    {
        m_curr_state.grid_size = size;
        for (int i = 0; i < TURNS_TO_CREATE; ++i)
        {
            m_possible_turns.push_back(new Turn(1, max_turn_depth));
        }
    }

    void Process() override
    {
        std::cerr << m_curr_state << std::endl;
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();

        //while (elapsed < 20)
        //{
        Move m = CalculateBestMove(m_max_turn_depth);

        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();
        //}
               
        std::cout << m << " " << m_turn_score << " points! " << elapsed << "ms" << std::endl;
        ++m_turn_count;
    }

    void Clear() override
    {
        m_curr_state.Clear();
    }

    void SetGameStateRow(int row_idx, const std::string &row_data) override
    {
        for(int i = 0; i < row_data.size(); ++i)
        { 
            m_curr_state.blocks[i][row_idx] = row_data[i] - 48;
        }
    }

    void SetGameStatePlayer(int player_idx, int x, int y) override
    {
        m_curr_state.player_loc[player_idx].x = x;
        m_curr_state.player_loc[player_idx].y = y;
    }

    void SetGameStateOpponent(int opp_idx, int x, int y) override
    {
        m_curr_state.opp_loc[opp_idx].x = x;
        m_curr_state.opp_loc[opp_idx].y = y;
    }

    void AddGameStateLegalAction(int legal_act_count, int idx, const std::string &type, int unit_idx, const std::string &move_dir, const std::string &build_dir) override
    {
        //m_population.at(idx).genome[0].genes[unit_idx].type = StringToType(type);
        //m_population.at(idx).genome[0].genes[unit_idx].move_push_dir = StringToDirection(move_dir);
        //m_population.at(idx).genome[0].genes[unit_idx].build_dir = StringToDirection(build_dir);
    }

    GameState & CurrentState() {
        return m_curr_state;
    }

    GameState & ResultState() {
        return m_result_state;
    }

private:
    void ClearPossibleTurns()
    {
        for (auto &t : m_possible_turns)
        {
            t->Clear();
        }
    }
    Move CalculateBestMove(int max_depth)
    {
        std::vector<Move> my_poss_moves[2];
        my_poss_moves[0] = GetPossiblePlayerMoves(m_curr_state, 0);
        my_poss_moves[1] = GetPossiblePlayerMoves(m_curr_state, 1);
        
        ClearPossibleTurns();
        int count = 0;
        for (int i = 0; i < 2; ++i)
        {
            for (auto & m : my_poss_moves[i])
            {
                m_possible_turns.at(count)->player_move = m;
                m_possible_turns.at(count)->is_valid = true;
                ++count;
            }
        }
        
        for (auto &t : m_possible_turns)
        {
            if (t->is_valid)
            {
                t->Simulate(m_curr_state, true, m_turn_count);
            }
            else
            {
                break;
            }
        }

        std::sort(m_possible_turns.begin(), m_possible_turns.end(), [](const Turn *a, const Turn *b)
        {
            return a->result.score > b->result.score;
        });

        for (int i = 0; i < std::min(static_cast<size_t>(75), m_possible_turns.size()); ++i)
        {
            if (m_possible_turns.at(i)->is_valid)
            {
                m_possible_turns.at(i)->GenerateSubTurns(max_depth, m_turn_count);
            }
        }

        std::sort(m_possible_turns.begin(), m_possible_turns.end(), [](const Turn *a, const Turn *b)
        {
            return a->AggregateScore() > b->AggregateScore();
        });

        if (!m_possible_turns.empty())
        {
            m_result_state = m_possible_turns[0]->result;
            m_turn_score = m_possible_turns[0]->score;
            return m_possible_turns[0]->player_move;
        }
        else
        {
            return Move();
        }
    }

    GameState m_curr_state;
    GameState m_result_state;
    double m_turn_score;
    int m_max_turn_depth;
    std::vector<Turn*> m_possible_turns;
    int m_turn_count;
};


void test()
{

    GameState st;
    
    GeneticAI ai(6, 2);
    
    
    ai.Clear();
                        // 01234567
    ai.SetGameStateRow(0, "101231");
    ai.SetGameStateRow(1, "110232");
    ai.SetGameStateRow(2, "12..22");
    ai.SetGameStateRow(3, "1....2");
    ai.SetGameStateRow(4, "0.22.0");
    ai.SetGameStateRow(5, "041220");

    auto test_vec = ai.CurrentState().DeathTraps();


    ai.SetGameStatePlayer(0, 5, 1);
    ai.SetGameStatePlayer(1, 4, 2);
    ai.SetGameStateOpponent(0, -1, -1);
    ai.SetGameStateOpponent(1, -1, -1);

    ai.Process();
}

void InferOppPositions(GameState & curr_state, GameState & prev_state) {

    bool got_pushed = 
        curr_state.player_loc[0] != prev_state.player_loc[0] ||
        curr_state.player_loc[1] != prev_state.player_loc[1];
    if(got_pushed) {
        for(int i = 0; i < 2; ++i) {
            if(curr_state.opp_loc[i] == Point(-1,-1)) {
                curr_state.opp_loc[i] = prev_state.opp_loc[i];
            }
        }
        return;
    }
}


GameState prev_state;

void run()
{
    int size;
    std::cin >> size; std::cin.ignore();
    int unitsPerPlayer;
    std::cin >> unitsPerPlayer; std::cin.ignore();
    
    GeneticAI ai(size, 2);
    
    bool first_turn = true;

    // game loop
    while (1) {
        ai.Clear();
        for (int i = 0; i < size; i++) {
            std::string row;
            std::cin >> row; std::cin.ignore();
            //std::cerr << row << std::endl;
            ai.SetGameStateRow(i, row);
        }
    
        for (int i = 0; i < unitsPerPlayer; i++) {
            int unitX;
            int unitY;
            std::cin >> unitX >> unitY; std::cin.ignore();
            ai.SetGameStatePlayer(i, unitX, unitY);
        }
        for (int i = 0; i < unitsPerPlayer; i++) {
            int otherX;
            int otherY;
            std::cin >> otherX >> otherY; std::cin.ignore();
            //std::cerr << "Opp[" << i << "] (" << otherX << "," << otherY << ")" << std::endl;
            ai.SetGameStateOpponent(i, otherX, otherY);
        }
        if(!first_turn) {
            InferOppPositions(ai.CurrentState(), prev_state);
        }
        int legalActions;
        std::cin >> legalActions; std::cin.ignore();
        std::cerr << "Legal Actions count: " << legalActions << std::endl;
        for (int i = 0; i < legalActions; i++) {
            std::string atype;            
            int index;
            std::string dir1;
            std::string dir2;
            std::cin >> atype >> index >> dir1 >> dir2; std::cin.ignore();
            ai.AddGameStateLegalAction(legalActions, i, atype, index, dir1, dir2);
        }
        
        ai.Process();
        prev_state = ai.ResultState();
        first_turn = false;
    }
}

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    //test();
    run();
    
}

