#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <array>
#include <map>

const int NUM_OF_MOL = 5;
const int NUM_OF_SAMPLES_HELD = 3;
const int NUM_OF_MOL_CAN_HOLD = 10;
const int MAX_MOLS_AVAIL = 5;

enum class Modules
{
    DIAGNOSIS,
    MOLECULES,
    LABORATORY,
    SAMPLES,
    NONE,
    NUM_OF_MODULES
};

enum CarriedBy
{
    CLOUD = -1,
    ME,
    OPPONENT
};


class Project
{
public:
    Project(int a, int b, int c, int d, int e)
        : exp_rqd()
    {
        exp_rqd[0] = a;
        exp_rqd[1] = b;
        exp_rqd[2] = c;
        exp_rqd[3] = d;
        exp_rqd[4] = e;

        score = std::accumulate(exp_rqd.begin(), exp_rqd.end(), 0);
    }
    std::array<int, NUM_OF_MOL> exp_rqd;
    int score;
};
std::vector<Project> project_library;

class Sample
{
public:
    Sample(int _id, int _carried_by, int _rank, const std::string &_gain, int _health, int ca, int cb, int cc, int cd, int ce)
        : id(_id)
        , carried_by(_carried_by)
        , rank(_rank)
        , gain(_gain)
        , mol_exp_gain(-1)
        , health(_health)
        , cost()
        , is_diagnosed(false)
        , score(-50.0)
        , times_ive_skipped(0)
    {
        cost[0] = ca;
        cost[1] = cb;
        cost[2] = cc;
        cost[3] = cd;
        cost[4] = ce;

        total_cost = std::accumulate(cost.begin(), cost.end(), 0);

        if (health >= 0) { is_diagnosed = true; }

        score = static_cast<double>(health) / static_cast<double>(total_cost);
        score *= (total_cost > NUM_OF_MOL_CAN_HOLD) ? -1.0 : 1.0;

        mol_exp_gain = static_cast<int>(gain[0]) - 65;
    }

    Sample()
        : id(-1)
        , carried_by(-1)
        , rank(-1)
        , gain("")
        , health(-1)
        , cost()
    {
    }

    int id;
    int carried_by;
    int rank;
    std::string gain;
    int mol_exp_gain;
    int health;
    int total_cost;
    std::array<int, NUM_OF_MOL> cost;
    bool is_diagnosed;
    double score;
    int times_ive_skipped;
};

std::map<int, Sample> sample_library;
std::vector<Sample*> samples_in_cloud;
std::vector<Sample*> samples_opponent_holds;
std::vector<Sample*> samples_i_hold;

std::array<int, 5> molecules_available;


void print_sample(Sample *s)
{
    std::cerr << "SAMPLE " << s->id << std::endl;
    std::cerr << "SCORE = " << s->score << std::endl;
    std::cerr << "GAIN = " << s->gain.c_str() << std::endl;
    std::cerr << "A = " << s->cost[0] << ", B = " << s->cost[1] << ", C = " << s->cost[2] << ", D = " << s->cost[3] << ", E = " << s->cost[4] << std::endl;
}

bool sample_vector_contains(std::vector<Sample*> &v, int sample_id)
{
    auto it = std::find_if(v.begin(), v.end(), [&](const Sample *s)
    {
        return s->id == sample_id;
    });

    return it != v.end();
}

class Robot
{
public:
    struct Input
    {
        Modules target_mod;
        int eta;
        int health;
        int turn;
        std::array<int, NUM_OF_MOL> mols_held;
        std::array<int, NUM_OF_MOL> expertise;

        Input()
            : target_mod(Modules::DIAGNOSIS)
            , eta(0)
            , health(0)
            , turn(0)
            , mols_held()
            , expertise()
        {
        }
    };

    Robot()
        : m_exe_fns()
        , m_goto_fns()
        , m_in()
        , m_state(Modules::NONE)
        , m_opponent(nullptr)
        , m_num_of_mols_to_acquire()
        , m_target_exp()
        , m_cloud_sample_to_grab(-1)
        , m_molecule_opp_is_hoarding()
        , m_sample_to_ditch(-1)
        , m_total_exp(0)
        , m_lowest_exp(0)
        , m_highest_exp(0)
        , m_lowest_exp_idx(0)
        , m_highest_exp_idx(0)
        , m_num_of_turns_waiting_for_mol(0)
        , m_is_intialized(false)
        , m_exp_targets_set(false)
        , m_dump_everything(false)
    {

        m_exe_fns[(size_t)Modules::NONE] = &Robot::EXE_NONE;
        m_exe_fns[(size_t)Modules::SAMPLES] = &Robot::EXE_SAMPLES;
        m_exe_fns[(size_t)Modules::DIAGNOSIS] = &Robot::EXE_DIAGNOSIS;
        m_exe_fns[(size_t)Modules::MOLECULES] = &Robot::EXE_MOLECULES;
        m_exe_fns[(size_t)Modules::LABORATORY] = &Robot::EXE_LABORATORY;

        m_goto_fns[(size_t)Modules::NONE] = &Robot::GOTO_NONE;
        m_goto_fns[(size_t)Modules::SAMPLES] = &Robot::GOTO_SAMPLES;
        m_goto_fns[(size_t)Modules::DIAGNOSIS] = &Robot::GOTO_DIAGNOSIS;
        m_goto_fns[(size_t)Modules::MOLECULES] = &Robot::GOTO_MOLECULES;
        m_goto_fns[(size_t)Modules::LABORATORY] = &Robot::GOTO_LABORATORY;

        m_insults.push_back("Bite me.");
        m_insults.push_back("Hey laser lips, your momma was a snowblower.");
        m_insults.push_back("Eat my shorts.");
        m_insults.push_back("Don't think I don't know what you're doing.");
        m_insults.push_back("Beep Blooop. I'm totally a robot.");
        m_insults.push_back("You suck, ya jackass!");
        m_insults.push_back("I'm going to hoard all of the D molecules.");
    }

    void SetOpponent(const Robot & opp)
    {
        m_opponent = &opp;
    }
    void Run()
    {
        std::cerr << "Turn " << m_in.turn << std::endl;
        m_total_exp = std::accumulate(m_in.expertise.begin(), m_in.expertise.end(), 0);
        m_total_mols_i_hold = std::accumulate(m_in.mols_held.begin(), m_in.mols_held.end(), 0);
        //m_lowest_exp_idx = std::min_element(m_in.expertise.begin(), m_in.expertise.end()) - m_in.expertise.begin();
        //m_highest_exp_idx = std::max_element(m_in.expertise.begin(), m_in.expertise.end()) - m_in.expertise.begin();
        //m_lowest_exp  = m_in.expertise[m_lowest_exp_idx];
        //m_highest_exp = m_in.expertise[m_highest_exp];

        if (m_in.turn > 200)
        {
            FindProjectGoals();
        }

        CheckForHoarding();

        if (m_in.eta == 0)
        {
            (this->*m_exe_fns[static_cast<size_t>(m_in.target_mod)])();
        }
        else
        {
            (this->*m_goto_fns[static_cast<size_t>(m_in.target_mod)])();
        }
    }

    Input & In() { return m_in; }
private:

    void FindProjectGoals()
    {
        int closest_proj_id = -1;
        int closest_proj_dev = 999;
        for (size_t proj = 0; proj < project_library.size(); ++proj)
        {
            Project &p = project_library[proj];
            int dev = 0;
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                dev += std::max(0, p.exp_rqd[mol] - m_in.expertise[mol]);
            }

            if (dev < closest_proj_dev)
            {
                closest_proj_dev = dev;
                closest_proj_id = proj;
            }
        }


        //if (closest_proj_dev < 1)
        //{
        //    m_target_exp[0] = project_library[closest_proj_id].exp_rqd[0];
        //    m_target_exp[1] = project_library[closest_proj_id].exp_rqd[1];
        //    m_target_exp[2] = project_library[closest_proj_id].exp_rqd[2];
        //    m_target_exp[3] = project_library[closest_proj_id].exp_rqd[3];
        //    m_target_exp[4] = project_library[closest_proj_id].exp_rqd[4];
        //
        //    m_exp_targets_set = true;
        //}
    }

    void CheckForHoarding()
    {
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            m_opponents_mol_history[mol].push_back(m_opponent->m_in.mols_held[mol]);
        }

        m_molecule_opp_is_hoarding.clear();
        if (m_in.turn > 50)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                double avg = static_cast<double>(std::accumulate(m_opponents_mol_history[mol].end() - std::min(std::max(m_in.turn / 2.0, 25.0), 50.0), m_opponents_mol_history[mol].end(), 0.0)) / static_cast<double>(std::min(std::max(m_in.turn / 2.0, 25.0), 50.0));
                //std::cerr << "Hoarding Average: " << avg << std::endl;
                if (avg > 3.0)
                {
                    m_molecule_opp_is_hoarding.push_back(mol);
                }
            }
        }

        if (!m_molecule_opp_is_hoarding.empty())
        {
            for (auto &m : m_molecule_opp_is_hoarding)
            {
                std::cerr << "The son of a bitch is hoarding " << static_cast<char>(m + 65) << "s." << std::endl;
            }
        }
    }

    void EXE_NONE()
    {
        std::cout << "GOTO SAMPLES" << std::endl;
    }

    void EXE_SAMPLES()
    {
        int s = ASampleExistsThatICanProcess();
        if (s < 0)
        {
            if ((samples_i_hold.size() < 3 && m_in.turn < 360) ||
                (samples_i_hold.size() < 2 && (m_in.turn < 370 || samples_i_hold[0]->times_ive_skipped > 0)) ||
                (samples_i_hold.size() < 1))

            {
                GetASample();
                return;
            }
            else
            {
                std::cerr << "I've been holding sample 0 for " << samples_i_hold[0]->times_ive_skipped << " times" << std::endl;
                std::cout << "GOTO DIAGNOSIS" << std::endl;
            }
        }
        else
        {
            if (samples_i_hold.size() < 2)
            {
                GetASample();
                return;
            }
            else
            { 
                m_cloud_sample_to_grab = s;
                std::cerr << "I'm going to get a pre-diagnosed sample." << std::endl;
                std::cout << "GOTO DIAGNOSIS" << std::endl;
            }
        }
    }

    void GetASample()
    {
        if (m_total_exp > 10 || m_in.turn > 310)
        {
            std::cout << "CONNECT 3" << std::endl;
        }
        else if (m_total_exp > 6 || m_in.turn > 240)
        {
            std::cout << "CONNECT 2" << std::endl;
        }
        else
        {
            std::cout << "CONNECT 1" << std::endl;
        }
    }

    int NumberOfRankIHold(int rank)
    {
        int count = 0;
        for (auto &s : samples_i_hold)
        {
            if (s->rank == rank)
            {
                count++;
            }
        }

        return count;
    }

    void EXE_DIAGNOSIS()
    {
        if (m_dump_everything)
        {
            if (!samples_i_hold.empty())
            {
                std::cout << "CONNECT " << samples_i_hold.back()->id << std::endl;
                return;
            }
            else
            {
                m_dump_everything = false;
                std::cout << "GOTO SAMPLES" << std::endl;
                return;
            }
        }

        if (m_cloud_sample_to_grab >= 0)
        {
            std::cerr << "Getting Cloud Sample..." << std::endl;
            if (samples_i_hold.size() > 2)
            {
                std::cout << "CONNECT " << samples_i_hold[2]->id << std::endl;
                return;
            }
            if (sample_vector_contains(samples_in_cloud, m_cloud_sample_to_grab))
            {
                std::cout << "CONNECT " << m_cloud_sample_to_grab << std::endl;
                m_cloud_sample_to_grab = -1;
                return;
            }
            else
            {
                for (auto & s : samples_i_hold)
                {
                    if (!s->is_diagnosed)
                    {
                        std::cout << "CONNECT " << s->id << std::endl;
                        return;
                    }
                }

                if (m_total_mols_i_hold == 0)
                {
                    std::cout << "GOTO MOLECULES" << std::endl;
                    return;
                }
                else
                {
                    std::cout << "GOTO LABORATORY" << std::endl;
                    return;
                }
            }
        }
        else if (m_sample_to_ditch >= 0)
        {
            std::cout << "CONNECT " << m_sample_to_ditch << std::endl;
            m_sample_to_ditch = -1;
            return;
        }
        else
        {
            for (auto & s : samples_i_hold)
            {
                if (!s->is_diagnosed)
                {
                    std::cout << "CONNECT " << s->id << std::endl;
                    return;
                }
            }

            for (auto & s : samples_i_hold)
            {
                if (!IsSampleIShouldProcess(s))
                {
                    std::cout << "CONNECT " << s->id << std::endl;
                    return;
                }
            }

            if (samples_i_hold.size() < 2 && m_in.turn < 350)
            {
                int sample_from_cloud = GetSampleFromCloud();
                if (sample_from_cloud < 0)
                {
                    std::cout << "GOTO SAMPLES" << std::endl;
                    return;
                }
                else
                {
                    std::cout << "CONNECT " << sample_from_cloud << std::endl;
                    return;
                }
                return;
            }
            else
            {
                if (CanProcessAnyOfMySamples() > 1)
                {
                    std::cout << "GOTO LABORATORY" << std::endl;
                }
                else
                {
                    std::cout << "GOTO MOLECULES" << std::endl;
                    return;
                }
            }
        }
    }

    int GetSampleFromCloud()
    {
        for (auto &s : samples_in_cloud)
        {
            if (IsSampleIShouldProcess(s) && s->rank > 1)
            {
                return s->id;
            }
        }

        return -1;
    }

    bool IsSampleIShouldProcess(Sample *s)
    {
        return IsPossibleToProcess(s) && !ContainsMoleculesBeingHoarded(s) && ContainsExpertiseINeed(s) && s->times_ive_skipped < 3;
    }

    bool ContainsExpertiseINeed(Sample *s)
    {
        bool result = (m_in.expertise[s->mol_exp_gain] < m_target_exp[s->mol_exp_gain]) || !m_exp_targets_set;
        if (!result)
        {
            std::cerr << "Rejecting due to not containing the EXP I want." << std::endl;
        }

        return result;
    }

    bool ContainsMoleculesBeingHoarded(Sample *s)
    {
        if (m_molecule_opp_is_hoarding.empty()) { return false; }

        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            if (s->cost[mol] - m_in.expertise[mol] > 0)
            {
                auto it = std::find_if(m_molecule_opp_is_hoarding.begin(), m_molecule_opp_is_hoarding.end(), [&](int m) { return m == mol; });
                if (it != m_molecule_opp_is_hoarding.end())
                {
                    std::cerr << "Rejecting due to containing a hoarded molecule." << std::endl;
                }
                return it != m_molecule_opp_is_hoarding.end();
            }
        }

        return false;
    }

    bool IsPossibleToProcess(Sample *s)
    {
        bool can_hold_enough = CalculateMoleculeCost(s) < NUM_OF_MOL_CAN_HOLD;

        bool enough_available = s->cost[0] - m_in.expertise[0] <= MAX_MOLS_AVAIL &&
            s->cost[1] - m_in.expertise[1] <= MAX_MOLS_AVAIL &&
            s->cost[2] - m_in.expertise[2] <= MAX_MOLS_AVAIL &&
            s->cost[3] - m_in.expertise[3] <= MAX_MOLS_AVAIL &&
            s->cost[4] - m_in.expertise[4] <= MAX_MOLS_AVAIL;

        if (!can_hold_enough || !enough_available)
        {
            if (!can_hold_enough && !enough_available)
            {
                std::cerr << "Rejecting due to can't hold enough and not enough available." << std::endl;
            }
            else if (!can_hold_enough)
            {
                std::cerr << "Rejecting t due to cant hold enough." << std::endl;
            }
            else if (!enough_available)
            {
                std::cerr << "Rejecting t due to not enough available." << std::endl;
            }
        }
        return can_hold_enough && enough_available;
    }

    int CalculateMoleculeCost(Sample *s)
    {
        int cost = 0;
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            cost += s->cost[mol] - m_in.expertise[mol];
        }

        return cost;
    }

    void BlockOppSample()
    {
        for (auto &s : samples_opponent_holds)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                m_num_of_mols_to_acquire[mol] += std::max(s->cost[mol], 0);
            }
        }

        int overage = std::accumulate(m_num_of_mols_to_acquire.begin(), m_num_of_mols_to_acquire.end(), 0) - NUM_OF_MOL_CAN_HOLD;
        if (overage > 0)
        {
            while (overage > 0)
            {
                for (int mol = 0; mol < NUM_OF_MOL; ++mol)
                {
                    if (m_num_of_mols_to_acquire[mol] > 0)
                    {
                        m_num_of_mols_to_acquire[mol]--;
                        --overage;
                        break;
                    }
                }
            }
        }
    }

    void EXE_MOLECULES()
    {
        std::cerr << "EXE MOLECULES" << std::endl;

        if (m_in.turn > 380 && samples_i_hold.empty())
        {
            BlockOppSample();
        }

        DecideWhichMoleculesToGet();

        if (m_total_mols_i_hold < NUM_OF_MOL_CAN_HOLD)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                if (m_in.mols_held[mol] < m_num_of_mols_to_acquire[mol] && molecules_available[mol] > 0)
                {
                    m_num_of_turns_waiting_for_mol = 0;
                    std::cout << "CONNECT " << static_cast<char>(mol + 65) << std::endl;
                    return;
                }
                else
                {
                    std::cerr << "I hold " << m_in.mols_held[mol] << " " << static_cast<char>(mol + 65) << " molecules. I need " << m_num_of_mols_to_acquire[mol] << std::endl;
                }
            }

            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                if (m_in.mols_held[mol] < m_num_of_mols_to_acquire[mol] + m_num_of_mols_to_hoard[mol] && molecules_available[mol] > 0)
                {
                    std::cout << "CONNECT " << static_cast<char>(mol + 65) << std::endl;
                    return;
                }
            }
        }
        else
        {
            std::cerr << "I cant hold any more molecules. I'm holding: " << m_total_mols_i_hold << std::endl;
        }

        if (IHaveEnoughToProcess() || CanProcessAnyOfMySamples() > 0)
        {
            m_num_of_turns_waiting_for_mol = 0;
            std::cout << "GOTO LABORATORY" << std::endl;
            return;
        }
        else
        {
            if (ASampleExistsThatICanProcess() >= 0 && m_opponent->m_in.target_mod != Modules::LABORATORY)
            {
                std::cerr << "I dont have enough to process, but a sample exists that I can process." << std::endl;
                m_num_of_turns_waiting_for_mol = 0;
                std::cout << "GOTO DIAGNOSIS" << std::endl;
                return;
            }
            else
            {
                for (auto &s : samples_i_hold)
                {
                    if (ContainsMoleculesBeingHoarded(s))
                    {
                        std::cerr << "This Contains Molecules that are beign hoarded." << std::endl;
                        m_num_of_turns_waiting_for_mol = 0;
                        m_sample_to_ditch = s->id;
                        std::cout << "GOTO DIAGNOSIS" << std::endl;
                        return;
                    }
                }

                if (m_num_of_turns_waiting_for_mol >= 6 && m_opponent->m_in.target_mod != Modules::LABORATORY)
                {
                    if (CanProcessAnyOfMySamples() > 0)
                    {
                        std::cout << "GOTO LABORATORY" << std::endl;
                        return;
                    }
                    std::cerr << "I've waited long enough, I'm getting rid of these." << std::endl;
                    m_dump_everything = true;
                    std::cout << "GOTO DIAGNOSIS" << std::endl;
                    return;
                }

                std::cout << "WAIT" << std::endl;
                m_num_of_turns_waiting_for_mol += 2;
            }
        }
    }

    int ASampleExistsThatICanProcess()
    {
        for (auto &s : samples_in_cloud)
        {
            if (CanProcessSample(s))
            {
                m_cloud_sample_to_grab = s->id;
                return s->id;
            }
        }

        return -1;
    }

    bool IHaveEnoughToProcess()
    {
        bool can_process = false;

        for (auto &s : samples_i_hold)
        {
            can_process = s->cost[0] - m_in.expertise[0] <= m_in.mols_held[0] &&
                s->cost[1] - m_in.expertise[1] <= m_in.mols_held[1] &&
                s->cost[2] - m_in.expertise[2] <= m_in.mols_held[2] &&
                s->cost[3] - m_in.expertise[3] <= m_in.mols_held[3] &&
                s->cost[4] - m_in.expertise[4] <= m_in.mols_held[4];

            if (can_process) { return true; }
        }

        return false;
    }

    void PrintNumOfMoleculesToAquire()
    {
        std::cerr << "Mols to get: A = " << m_num_of_mols_to_acquire[0] << ", B= " << m_num_of_mols_to_acquire[1] << ", C= " << m_num_of_mols_to_acquire[2] << ", D= " << m_num_of_mols_to_acquire[3] << ", E= " << m_num_of_mols_to_acquire[4] << std::endl;
    }


    void HoardMolecules()
    {
        m_num_of_mols_to_hoard.fill(0);

        if (m_total_exp > 6) { return; }

        int remain_space = NUM_OF_MOL_CAN_HOLD - std::accumulate(m_num_of_mols_to_acquire.begin(), m_num_of_mols_to_acquire.end(), 0);
        struct exp
        {
            int idx;
            int val;
        };

        std::array<exp, NUM_OF_MOL> sorted_exp;
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            sorted_exp[mol].idx = mol;
            sorted_exp[mol].val = m_in.expertise[mol];
        }

        std::sort(sorted_exp.begin(), sorted_exp.end(), [](const exp& e1, const exp& e2) { return e1.val < e2.val; });
        

        while (remain_space > 0)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                m_num_of_mols_to_hoard[mol] += 1;
                remain_space--;
                if (remain_space == 0)
                {
                    break;
                }
            }
        }
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            std::cerr << " Hoarding " << static_cast<char>(mol + 65) << "s = " << m_num_of_mols_to_hoard[mol] << std::endl;
        }
    }

    void DecideWhichMoleculesToGet()
    {
        OptimizeMySampleOrder();
        m_num_of_mols_to_acquire.fill(0);
        for (auto &s : samples_i_hold)
        {
            //print_sample(s);
        }

        if (CanIHoldEnoughMolsForAll())
        {
            std::cerr << "I can get enough for all of these." << std::endl;
            m_simulated_expertise = m_in.expertise;
            for (auto &s : samples_i_hold)
            {
                if (ContainsMoleculesBeingHoarded(s)) { continue; }
                for (int mol = 0; mol < NUM_OF_MOL; ++mol)
                {
                    m_num_of_mols_to_acquire[mol] += std::max(s->cost[mol] - m_simulated_expertise[mol], 0);
                }

                if (s->mol_exp_gain >= 0)
                {
                    m_simulated_expertise[s->mol_exp_gain]++;
                }
            }

            PrintNumOfMoleculesToAquire();
        }
        else
        {
            if (samples_i_hold.size() > 1 && CanIHoldEnoughMolForTwo(samples_i_hold[0], samples_i_hold[1]))
            {
                std::cerr << "I can get enough for 0 and 1." << std::endl;
                GetNeededMolsFor(samples_i_hold[0], samples_i_hold[1]);
                PrintNumOfMoleculesToAquire();
            }
            else if (samples_i_hold.size() > 2 && CanIHoldEnoughMolForTwo(samples_i_hold[1], samples_i_hold[2]))
            {
                std::cerr << "I can get enough for 1 and 2." << std::endl;
                GetNeededMolsFor(samples_i_hold[1], samples_i_hold[2]);
                PrintNumOfMoleculesToAquire();
            }
            else if (samples_i_hold.size() > 2 && CanIHoldEnoughMolForTwo(samples_i_hold[0], samples_i_hold[2]))
            {
                std::cerr << "I can get enough for 0 and 2." << std::endl;
                GetNeededMolsFor(samples_i_hold[0], samples_i_hold[2]);
                PrintNumOfMoleculesToAquire();
            }
            else
            {
                std::cerr << "I can get enough for 0." << std::endl;
                GetNeededMolsFor(samples_i_hold[0]);
                PrintNumOfMoleculesToAquire();
            }
        }

        int total_mols_to_aq = std::accumulate(m_num_of_mols_to_acquire.begin(), m_num_of_mols_to_acquire.end(), 0);

        HoardMolecules();
    }

    void GetNeededMolsFor(Sample *s)
    {
        m_num_of_mols_to_acquire.fill(0);

        if (ContainsMoleculesBeingHoarded(s)) { return; }

        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            m_num_of_mols_to_acquire[mol] += std::max(s->cost[mol] - m_in.expertise[mol], 0);
        }
    }

    void GetNeededMolsFor(Sample *s1, Sample *s2)
    {
        m_num_of_mols_to_acquire.fill(0);
        m_simulated_expertise = m_in.expertise;

        if (!ContainsMoleculesBeingHoarded(s1))
        {

            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                m_num_of_mols_to_acquire[mol] += std::max(s1->cost[mol] - m_simulated_expertise[mol], 0);
            }

            if (s1->mol_exp_gain >= 0)
            {
                m_simulated_expertise[s1->mol_exp_gain]++;
            }
        }

        if (!ContainsMoleculesBeingHoarded(s2))
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                m_num_of_mols_to_acquire[mol] += std::max(s2->cost[mol] - m_simulated_expertise[mol], 0);

                m_num_of_mols_to_acquire[mol] = std::min(m_num_of_mols_to_acquire[mol], MAX_MOLS_AVAIL);
            }
        }
    }

    bool CanIHoldEnoughMolForTwo(Sample *s1, Sample *s2)
    {
        int total_cost = 0;
        std::array<int, NUM_OF_MOL> cost = { 0, 0, 0, 0, 0 };
        m_simulated_expertise = m_in.expertise;
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            cost[mol] += std::max(s1->cost[mol] - m_simulated_expertise[mol], 0);
            total_cost += std::max(s1->cost[mol] - m_simulated_expertise[mol], 0);
        }

        if (s1->mol_exp_gain >= 0)
        {
            m_simulated_expertise[s1->mol_exp_gain]++;
        }

        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            cost[mol] += std::max(s2->cost[mol] - m_simulated_expertise[mol], 0);
            total_cost += std::max(s2->cost[mol] - m_simulated_expertise[mol], 0);
        }


        for (auto &c : cost)
        {
            if (c > MAX_MOLS_AVAIL) { return false; }
        }

        return total_cost <= NUM_OF_MOL_CAN_HOLD;
    }

    bool CanIHoldEnoughMolsForAll()
    {
        std::cerr << "CanIHoldEnoughMolsForAll" << std::endl;
        int total_cost = 0;
        m_simulated_expertise = m_in.expertise;
        std::array<int, NUM_OF_MOL> cost = { 0, 0, 0, 0, 0 };
        for (auto &s : samples_i_hold)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                cost[mol] += std::max(s->cost[mol] - m_simulated_expertise[mol], 0);
                total_cost += std::max(s->cost[mol] - m_simulated_expertise[mol], 0);
            }
            if (s->mol_exp_gain >= 0)
            {
                m_simulated_expertise[s->mol_exp_gain]++;
            }
        }

        for (auto &c : cost)
        {
            if (c > MAX_MOLS_AVAIL) { return false; }
        }

        return total_cost <= NUM_OF_MOL_CAN_HOLD;
    }

    void EXE_LABORATORY()
    {
        int idx = 0;
        for (auto &s : samples_i_hold)
        {
            if (CanProcessSample(s))
            {
                std::cout << "CONNECT " << s->id << std::endl;
                return;
            }
            else
            {
                std::cerr << "Skipped sample " << s->id << std::endl;
                s->times_ive_skipped++;
                std::cerr << "I've been holding sample " << s->id << " for " << s->times_ive_skipped << "trips." << std::endl;
            }
            ++idx;
        }

        if (samples_i_hold.size() > 1 || (samples_i_hold.size() == 1 && m_in.turn > 360))
        {
            std::cout << "GOTO MOLECULES" << std::endl;
            return;
        }
        else if (m_in.turn > 370)
        {
            int s = ASampleExistsThatICanProcess();
            if (s >= 0)
            {
                m_cloud_sample_to_grab = s;
                std::cout << "GOTO DIAGNOSIS" << std::endl;
                return;
            }
            else
            {
                std::cout << "GOTO MOLECULES" << std::endl;
            }
        }
        else
        {
            std::cout << "GOTO SAMPLES" << std::endl; return;
        }
    }

    int CanProcessAnyOfMySamples()
    {
        int count = 0;
        for (auto &s : samples_i_hold)
        {
            if (CanProcessSample(s))
            {
                count++;
            }
        }

        return count;
    }

    bool CanProcessSample(Sample *sample)
    {
        for (int mol = 0; mol < NUM_OF_MOL; ++mol)
        {
            if ((sample->cost[mol] - m_in.expertise[mol]) > m_in.mols_held[mol]) { return false; }
        }

        return true;
    }

    void GOTO_NONE()
    {
        std::cout << RandomInsult() << std::endl;
    }

    void GOTO_SAMPLES()
    {
        std::cout << RandomInsult() << std::endl;
    }

    void GOTO_DIAGNOSIS()
    {
        std::cout << RandomInsult() << std::endl;
    }

    void GOTO_MOLECULES()
    {
        std::cout << RandomInsult() << std::endl;
    }

    void GOTO_LABORATORY()
    {
        std::cout << RandomInsult() << std::endl;
    }

    std::string RandomInsult()
    {
        srand(m_in.turn + m_total_exp);
        int idx = rand() % m_insults.size();
        return m_insults[idx];
    }

    void OptimizeMySampleOrder()
    {
        Sample *buff[3];

        std::array<int, 3> sample_order = { 0, 1, 2 };
        std::array<int, 3> best_sample_order = sample_order;

        int lowest_cost = 9999;

        int cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        sample_order[0] = 0;
        sample_order[0] = 2;
        sample_order[0] = 1;
        cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        sample_order[0] = 1;
        sample_order[0] = 0;
        sample_order[0] = 2;
        cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        sample_order[0] = 1;
        sample_order[0] = 2;
        sample_order[0] = 0;
        cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        sample_order[0] = 2;
        sample_order[0] = 0;
        sample_order[0] = 1;
        cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        sample_order[0] = 2;
        sample_order[0] = 1;
        sample_order[0] = 0;
        cost = CalculateMoleculeCost();
        if (cost < lowest_cost)
        {
            lowest_cost = cost;
            best_sample_order = sample_order;
        }

        buff[0] = samples_i_hold[0];
        buff[1] = samples_i_hold[1];
        buff[2] = samples_i_hold[2];

        samples_i_hold[0] = buff[best_sample_order[0]];
        samples_i_hold[1] = buff[best_sample_order[1]];
        samples_i_hold[2] = buff[best_sample_order[2]];
    }

    int CalculateMoleculeCost()
    {
        int total_cost = 0;
        m_simulated_expertise = m_in.expertise;
        for (auto &s : samples_i_hold)
        {
            for (int mol = 0; mol < NUM_OF_MOL; ++mol)
            {
                total_cost += std::max(s->cost[mol] - m_simulated_expertise[mol], 0);
            }
            if (s->mol_exp_gain >= 0)
            {
                m_simulated_expertise[s->mol_exp_gain]++;
            }
        }

        return total_cost;
    }

    typedef void (Robot::*FUNC)(void);
    std::array<FUNC, static_cast<size_t>(Modules::NUM_OF_MODULES)> m_exe_fns;
    std::array<FUNC, static_cast<size_t>(Modules::NUM_OF_MODULES)> m_goto_fns;
    Input m_in;
    Modules m_state;
    const Robot *m_opponent;
    std::array<int, NUM_OF_MOL> m_num_of_mols_to_acquire;
    std::array<int, NUM_OF_MOL> m_num_of_mols_to_hoard;
    std::array<int, NUM_OF_MOL> m_simulated_expertise;
    std::array<std::vector<int>, NUM_OF_MOL> m_opponents_mol_history;
    std::array<int, NUM_OF_MOL> m_target_exp;
    std::vector<int> m_molecule_opp_is_hoarding;
    int m_cloud_sample_to_grab;
    int m_sample_to_ditch;
    int m_total_exp;
    int m_lowest_exp;
    int m_highest_exp;
    int m_lowest_exp_idx;
    int m_highest_exp_idx;
    int m_num_of_turns_waiting_for_mol;
    int m_total_mols_i_hold;
    bool m_is_intialized;
    bool m_exp_targets_set;
    bool m_dump_everything;
    std::vector<std::string> m_insults;
};

Modules StringToModule(const std::string &str)
{
    if (str.compare("DIAGNOSIS") == 0) { return Modules::DIAGNOSIS; }
    else if (str.compare("LABORATORY") == 0) { return Modules::LABORATORY; }
    else if (str.compare("MOLECULES") == 0) { return Modules::MOLECULES; }
    else if (str.compare("SAMPLES") == 0) { return Modules::SAMPLES; }
    else { return Modules::NONE; }
}

/**
* Bring data on patient samples from the diagnosis machine to the laboratory with enough molecules to produce medicine!
**/
int main()
{
    Robot robots[2];
    robots[0].SetOpponent(robots[1]);

    int projectCount;
    std::cin >> projectCount; std::cin.ignore();
    project_library.clear();
    std::cerr << "Project Count: " << projectCount << std::endl;
    for (int i = 0; i < projectCount; i++) {
        int a;
        int b;
        int c;
        int d;
        int e;
        std::cin >> a >> b >> c >> d >> e; std::cin.ignore();
        project_library.push_back(Project(a, b, c, d, e));
        std::cerr << "Project " << i << ": " << a << ", " << b << ", " << c << ", " << d << ", " << e << std::endl;
    }

    // game loop
    while (1) {
        for (int i = 0; i < 2; i++) {
            std::string target;
            int eta;
            int score;
            int storageA;
            int storageB;
            int storageC;
            int storageD;
            int storageE;
            int expertiseA;
            int expertiseB;
            int expertiseC;
            int expertiseD;
            int expertiseE;
            std::cin >> target >> eta >> score >> storageA >> storageB >> storageC >> storageD >> storageE >> expertiseA >> expertiseB >> expertiseC >> expertiseD >> expertiseE; std::cin.ignore();
            robots[i].In().target_mod = StringToModule(target);
            robots[i].In().eta = eta;
            robots[i].In().health = score;
            robots[i].In().mols_held[0] = storageA;
            robots[i].In().mols_held[1] = storageB;
            robots[i].In().mols_held[2] = storageC;
            robots[i].In().mols_held[3] = storageD;
            robots[i].In().mols_held[4] = storageE;
            robots[i].In().expertise[0] = expertiseA;
            robots[i].In().expertise[1] = expertiseB;
            robots[i].In().expertise[2] = expertiseC;
            robots[i].In().expertise[3] = expertiseD;
            robots[i].In().expertise[4] = expertiseE;

        }

        int availableA;
        int availableB;
        int availableC;
        int availableD;
        int availableE;
        std::cin >> availableA >> availableB >> availableC >> availableD >> availableE; std::cin.ignore();

        molecules_available[0] = availableA;
        molecules_available[1] = availableB;
        molecules_available[2] = availableC;
        molecules_available[3] = availableD;
        molecules_available[4] = availableE;

        int sampleCount;
        std::cin >> sampleCount; std::cin.ignore();
        sample_library.clear();
        samples_in_cloud.clear();
        samples_i_hold.clear();
        samples_opponent_holds.clear();
        for (int i = 0; i < sampleCount; i++) {
            int sampleId;
            int carriedBy;
            int rank;
            std::string expertiseGain;
            int health;
            int costA;
            int costB;
            int costC;
            int costD;
            int costE;
            std::cin >> sampleId >> carriedBy >> rank >> expertiseGain >> health >> costA >> costB >> costC >> costD >> costE; std::cin.ignore();
            sample_library.emplace(sampleId, Sample(sampleId, carriedBy, rank, expertiseGain, health, costA, costB, costC, costD, costE));

            if (carriedBy == CLOUD) { samples_in_cloud.push_back(&sample_library[sampleId]); }
            if (carriedBy == OPPONENT) { samples_opponent_holds.push_back(&sample_library[sampleId]); }
            if (carriedBy == ME) { samples_i_hold.push_back(&sample_library[sampleId]); }
        }
        std::sort(samples_in_cloud.begin(), samples_in_cloud.end(), [](const Sample *a, const Sample * b)
        {
            return a->score < b->score;
        });

        std::sort(samples_opponent_holds.begin(), samples_opponent_holds.end(), [](const Sample *a, const Sample * b)
        {
            return a->score < b->score;
        });

        std::sort(samples_i_hold.begin(), samples_i_hold.end(), [](const Sample *a, const Sample * b)
        {
            return a->score < b->score;
        });

        std::cerr << "Number of samples I have: " << samples_i_hold.size() << std::endl;
        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;

        robots[0].In().turn += 2;
        robots[0].Run();
    }
}