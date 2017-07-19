#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>


struct Graph
{
    std::map<int, std::list<int>> nodes;

    void Initialize(int number_of_relations)
    {
        
    }

    void AddRelationship(int a, int b)
    {
        nodes[a].push_back(b);
        nodes[b].push_back(a);
    }

    int Analyze()
    {
        std::cerr << "Number of Nodes: " << nodes.size() << std::endl;
        
        int depth = 0;
        while (nodes.size() > 1)
        {
            std::vector<int> list_to_del;
            for (auto &node : nodes)
            {
                if (node.second.size() == 1)
                {
                    list_to_del.push_back(node.first);
                }
            }

            for (auto & n : list_to_del)
            {
                if (!nodes.at(n).empty())
                {
                    int parent = *nodes.at(n).begin();
                    nodes.at(parent).remove(n);
                }
                nodes.erase(n);
            }

            ++depth;
            //std::cout << "Depth = " << depth << std::endl;
        }

        return depth;
    }
};

void TestCase1()
{
    Graph g;
    
    g.AddRelationship(0, 1);
    g.AddRelationship(1, 2);
    g.AddRelationship(2, 3);
    g.AddRelationship(2, 4);

    std::cout << ((g.Analyze() == 2) ? "Passed" : "Failed") << std::endl;
}

void TestCase2()
{
    Graph g;
    g.AddRelationship(0, 1);
    g.AddRelationship(1, 2);
    g.AddRelationship(1, 4);
    g.AddRelationship(2, 3);
    g.AddRelationship(4, 5);
    g.AddRelationship(4, 6);
    
    std::cout << ((g.Analyze() == 3) ? "Passed" : "Failed") << std::endl;
}

void TestCase3()
{
    Graph g;
    g.AddRelationship(0, 1);
    g.AddRelationship(0, 8);
    g.AddRelationship(0, 15);
    g.AddRelationship(1, 2);
    g.AddRelationship(1, 5);
    g.AddRelationship(2, 3);
    g.AddRelationship(2, 4);
    g.AddRelationship(5, 6);
    g.AddRelationship(5, 7);
    g.AddRelationship(8, 9);
    g.AddRelationship(8, 12);
    g.AddRelationship(9, 10);
    g.AddRelationship(9, 11);
    g.AddRelationship(12, 13);
    g.AddRelationship(12, 14);
    g.AddRelationship(15, 16);
    g.AddRelationship(15, 19);
    g.AddRelationship(16, 17);
    g.AddRelationship(16, 18);
    g.AddRelationship(19, 20);
    g.AddRelationship(19, 21);

    std::cout << ((g.Analyze() == 3) ? "Passed" : "Failed") << std::endl;
}

int main()
{
    //TestCase1();
    //TestCase2();
    //TestCase3();
    Graph g;
    int n; // the number of adjacency relations
    std::cin >> n; std::cin.ignore();
    std::cerr << n << " relationships." << std::endl;
    for (int i = 0; i < n; i++) {
        int xi; // the ID of a person which is adjacent to yi
        int yi; // the ID of a person which is adjacent to xi
        std::cin >> xi >> yi; std::cin.ignore();
        //std::cerr << xi << " " << yi << std::endl;
        g.AddRelationship(xi, yi);
    }

    // Write an action using cout. DON'T FORGET THE "<< endl"
    // To debug: cerr << "Debug messages..." << endl;


    // The minimal amount of steps required to completely propagate the advertisement
    std::cout << g.Analyze() << std::endl;
}