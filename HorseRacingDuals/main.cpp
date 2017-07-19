// Casablanca’s hippodrome is organizing a new type of horse racing: duals. During a dual, only two horses will participate in the race. In order for the race to be interesting, it is necessary to try to select two horses with similar strength.
// 
// Write a program which, using a given number of strengths, identifies the two closest strengths and shows their difference with an integer (≥ 0).
//  	Game Input
// 
// Input
// Line 1: Number N of horses
// 
// The N following lines: the strength Pi of each horse. Pi is an integer.
// 
// Output
// The difference D between the two closest strengths. D is an integer greater than or equal to 0.
// Constraints
// 1 < N  < 100000
// 0 < Pi ≤ 10000000
// Example
// Input
// 3
// 5
// 8
// 9
// Output
// 1

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    std::vector<int> list;
    int N;
    cin >> N; cin.ignore();
    for (int i = 0; i < N; i++) {
        int Pi;
        cin >> Pi; cin.ignore();
        list.push_back(Pi);
    }


    std::sort(list.begin(), list.end());

    int smallest_diff = 10000000;
    for (size_t i = 0; i < list.size() - 1; ++i)
    {
        int diff = std::abs(list.at(i) - list.at(i + 1));
        if (diff < smallest_diff)
        {
            smallest_diff = diff;
        }
    }

    // Write an action using cout. DON'T FORGET THE "<< endl"
    // To debug: cerr << "Debug messages..." << endl;

    cout << smallest_diff << endl;
}