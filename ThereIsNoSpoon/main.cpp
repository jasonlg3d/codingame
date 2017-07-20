// The game is played on a rectangular grid with a given size. Some cells contain power nodes. The rest of the cells are empty.
// 
// The goal is to find, when they exist, the horizontal and vertical neighbors of each node.
//  	Rules
// 
// To do this, you must find each (x1,y1) coordinates containing a node, and display the (x2,y2) coordinates of the next node to the right, and the (x3,y3) coordinates of the next node to the bottom within the grid.
// 
// If a neighbor does not exist, you must output the coordinates -1 -1 instead of (x2,y2) and/or (x3,y3).
// 
// You lose if:
// You give an incorrect neighbor for a node.
// You give the neighbors for an empty cell.
// You compute the same node twice.
// You forget to compute the neighbors of a node.
// Don’t forget to run the tests by launching them from the “Test cases” window.
// 
// Warning: the tests provided are similar to the validation tests used to compute the final score but remain different. This is a "hardcoding" prevention mechanism. Harcoded solutions will not get any points.
// 
// Regarding the viewer, note that:
// A debug mode is available from the settings panel (the dented wheel)
// You can zoom/unzoom with the mouse wheel and move using drag'n drop (useful for large grids)
//  	Game Input
// 
// The program must first read the initialization data from standard input. Then, provide to the standard output one line per instruction.
// Initialization input
// Line 1: one integer width for the number of cells along the x axis.
// 
// Line 2: one integer height for the number of cells along the y axis.
// 
// Next height lines: A string  line  containing  width  characters. A dot . represents an empty cell. A zero 0 represents a cell containing a node.
// 
// Output for one game turn
// One line per node. Six integers on each line:   x1  y1  x2  y2  x3  y3
// 
// Where:
// (x1,y1) the coordinates of a node
// (x2,y2) the coordinates of the closest neighbor on the right of the node
// (x3,y3) the coordinates of the closest bottom neighbor
// If there is no neighbor, the coordinates should be -1 -1.
// Constraints
// 0 < width ≤ 30
// 0 < height ≤ 30
// 0 ≤ x1 < width
// 0 ≤ y1 < height
// -1 ≤ x2, x3 < width
// -1 ≤ y2, y3 < height
// Alloted response time to first output line ≤ 1s.
// Response time between two output lines ≤ 100ms

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

struct Point
{
    Point(int _x, int _y)
        : x(_x)
        , y(_y)
    {}

    int x, y;
};


/**
 * Don't let the machines win. You are humanity's last hope...
 **/
int main()
{
    std::vector<Point> points;
    int width; // the number of cells on the X axis
    std::cin >> width; std::cin.ignore();
    int height; // the number of cells on the Y axis
    std::cin >> height; std::cin.ignore();
    for (int y = 0; y < height; y++) {
        std::string line; // width characters, each either 0 or .
        std::getline(std::cin, line);
        std::cerr << line << std::endl;
        for (size_t x = 0; x < line.size(); ++x)
        {
            if (line.at(x) == '0')
            {
                points.push_back(Point(x, y));
            }
        }
    }

    //points.push_back(Point(0,0));
    //points.push_back(Point(2,0));
    //points.push_back(Point(4,0));
    //
    //int width = 5;
    //int height = 1;
    for (auto & pt : points)
    {
        int x_offset = 1;
        auto found_x_neighbor = std::find_if(points.begin(), points.end(), [&](const Point &p) {
                return p.x == pt.x + x_offset && p.y == pt.y;
            });
        while (x_offset < width)
        {
            found_x_neighbor = std::find_if(points.begin(), points.end(), [&](const Point &p) {
                return p.x == pt.x + x_offset && p.y == pt.y;
            });
            ++x_offset;

            if (found_x_neighbor != points.end()) { break; }
        }

        int y_offset = 1;
        auto found_y_neighbor = std::find_if(points.begin(), points.end(), [&](const Point &p) {
                return p.y == pt.y + y_offset && p.x == pt.x;
            });
        while (y_offset < height)
        {
            found_y_neighbor = std::find_if(points.begin(), points.end(), [&](const Point &p) {
                return p.y == pt.y + y_offset && p.x == pt.x;
            });
            ++y_offset;

            if (found_y_neighbor != points.end()) { break; }
        }
        
        int x2 = -1;
        int y2 = -1;
        int x3 = -1;
        int y3 = -1;

        if (found_x_neighbor != points.end())
        {
            x2 = found_x_neighbor->x;
            y2 = found_x_neighbor->y;
        }

        if (found_y_neighbor != points.end())
        {
            x3 = found_y_neighbor->x;
            y3 = found_y_neighbor->y;
        }

        std::cout << pt.x << " " << pt.y << " " << x2 << " " << y2 << " " << x3 << " " << y3 << std::endl;
    }
}