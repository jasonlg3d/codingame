// In stations and airports you often see this type of screen:
// 
// Have you ever asked yourself how it might be possible to simulate this display on a good old terminal? We have: with ASCII art!
// Rules
// 
// ASCII art allows you to represent forms by using characters. To be precise, in our case, these forms are words. For example, the word "MANHATTAN" could be displayed as follows in ASCII art:
// 
//  
// # #  #  ### # #  #  ### ###  #  ###
// ### # # # # # # # #  #   #  # # # #
// ### ### # # ### ###  #   #  ### # #
// # # # # # # # # # #  #   #  # # # #
// # # # # # # # # # #  #   #  # # # #
//  
//​ Your mission is to write a program that can display a line of text in ASCII art in a style you are given as input.
// 
//  	Game Input
// 
// Input
// Line 1: the width L of a letter represented in ASCII art. All letters are the same width.
// 
// Line 2: the height H of a letter represented in ASCII art. All letters are the same height.
// 
// Line 3: The line of text T, composed of N ASCII characters.
// 
// Following lines: the string of characters ABCDEFGHIJKLMNOPQRSTUVWXYZ? Represented in ASCII art.
// 
// Output
// The text T in ASCII art.
// The characters a to z are shown in ASCII art by their equivalent in upper case.
// The characters that are not in the intervals [a-z] or [A-Z] will be shown as a question mark in ASCII art.
// Constraints
// 0 < L < 30
// 0 < H < 30
// 0 < N < 200
// Example 1
// Input
// 4
// 5 
// E
//  #  ##   ## ##  ### ###  ## # # ###  ## # # #   # # ###  #  ##   #  ##   ## ### # # # # # # # # # # ### ### 
// # # # # #   # # #   #   #   # #  #    # # # #   ### # # # # # # # # # # #    #  # # # # # # # # # #   #   # 
// ### ##  #   # # ##  ##  # # ###  #    # ##  #   ### # # # # ##  # # ##   #   #  # # # # ###  #   #   #   ## 
// # # # # #   # # #   #   # # # #  #  # # # # #   # # # # # # #    ## # #   #  #  # # # # ### # #  #  #       
// # # ##   ## ##  ### #    ## # # ###  #  # # ### # # # #  #  #     # # # ##   #  ###  #  # # # #  #  ###  #  
// Output
// ### 
// #   
// ##  
// #   
// ### 
// Example 2
// Input
// 4
// 5
// MANHATTAN
//  #  ##   ## ##  ### ###  ## # # ###  ## # # #   # # ###  #  ##   #  ##   ## ### # # # # # # # # # # ### ### 
// # # # # #   # # #   #   #   # #  #    # # # #   ### # # # # # # # # # # #    #  # # # # # # # # # #   #   # 
// ### ##  #   # # ##  ##  # # ###  #    # ##  #   ### # # # # ##  # # ##   #   #  # # # # ###  #   #   #   ## 
// # # # # #   # # #   #   # # # #  #  # # # # #   # # # # # # #    ## # #   #  #  # # # # ### # #  #  #       
// # # ##   ## ##  ### #    ## # # ###  #  # # ### # # # #  #  #     # # # ##   #  ###  #  # # # #  #  ###  #  
// Output
// # #  #  ### # #  #  ### ###  #  ###  
// ### # # # # # # # #  #   #  # # # #  
// ### ### # # ### ###  #   #  ### # #  
// # # # # # # # # # #  #   #  # # # #  
// # # # # # # # # # #  #   #  # # # #
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <locale>
#include <algorithm>

using namespace std;

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

struct ASCIIChar
{
    std::vector<std::string> lines;
};

void test()
{
    std::array<ASCIIChar, 27> chars;
    int L = 4;
    int H = 5;
    std::string T = "TEST2MANHATTAN";
    string input_row[5];
    input_row[0] = " #  ##   ## ##  ### ###  ## # # ###  ## # # #   # # ###  #  ##   #  ##   ## ### # # # # # # # # # # ### ### ";
    input_row[1] = "# # # # #   # # #   #   #   # #  #    # # # #   ### # # # # # # # # # # #    #  # # # # # # # # # #   #   # ";
    input_row[2] = "### ##  #   # # ##  ##  # # ###  #    # ##  #   ### # # # # ##  # # ##   #   #  # # # # ###  #   #   #   ## ";
    input_row[3] = "# # # # #   # # #   #   # # # #  #  # # # # #   # # # # # # #    ## # #   #  #  # # # # ### # #  #  #       ";
    input_row[4] = "# # ##   ## ##  ### #    ## # # ###  #  # # ### # # # #  #  #     # # # ##   #  ###  #  # # # #  #  ###  #  ";

    for (int i = 0; i < H; ++i)
    {
        std::string ROW;
        ROW = input_row[i];
        int next_char = 0;
        for (int j = 0; j < L * 27; j+=L)
        {
            chars.at(next_char).lines.push_back(ROW.substr(j, L));
            ++next_char;
        }
    }


    for (int j = 0; j < H; ++j)
    {
        for (int i = 0; i < T.size(); ++i)
        {
            int char_index = T.at(i) - 65;
            if (char_index < 0 || char_index > 26) { char_index = 26; }
            std::cout << chars.at(char_index).lines.at(j);
        }
        std::cout << std::endl;
    }
}

int main()
{
    //test();
    std::array<ASCIIChar, 27> chars;
    int L;
    cin >> L; cin.ignore();
    int H;
    cin >> H; cin.ignore();
    string T;
    getline(cin, T);
    for (int i = 0; i < H; i++) {
        string ROW;
        getline(cin, ROW);

        int next_char = 0;
        for (int j = 0; j < L * 27; j+=L)
        {
            chars.at(next_char).lines.push_back(ROW.substr(j, L));
            ++next_char;
        }
    }
    std::locale loc;
    for (std::string::size_type i = 0; i < T.length(); ++i)
    {
        T.at(i) = std::toupper(T.at(i), loc);
    }
    for (int j = 0; j < H; ++j)
    {
        for (int i = 0; i < T.size(); ++i)
        {
            int char_index = T.at(i) - 65;
            if (char_index < 0 || char_index > 26) { char_index = 26; }
            std::cout << chars.at(char_index).lines.at(j);
        }
        std::cout << std::endl;
    }
}