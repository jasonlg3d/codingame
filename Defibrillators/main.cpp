// The city of Montpellier has equipped its streets with defibrillators to help save victims of cardiac arrests. The data corresponding to the position of all defibrillators is available online.
// 
// Based on the data we provide in the tests, write a program that will allow users to find the defibrillator nearest to their location using their mobile phone.
//  	Rules
// 
// The input data you require for your program is provided in text format.
// This data is comprised of lines, each of which represents a defibrillator. Each defibrillator is represented by the following fields:
// A number identifying the defibrillator
// Name
// Address
// Contact Phone number
// Longitude (degrees)
// Latitude (degrees)
// These fields are separated by a semicolon (;).
// 
// Beware: the decimal numbers use the comma (,) as decimal separator. Remember to turn the comma (,) into dot (.) if necessary in order to use the data in your program.
//  
// DISTANCE
// The distance d between two points A and B will be calculated using the following formula:
// 
// 
// x = (longitudeB - longitudeA) * cos((latitudeA + latitudeB) / 2)
// y = (latitudeB - latitudeA)
// d = sqrt(x^2 + y^2) * 6371
//
// Note: In this formula, the latitudes and longitudes are expressed in radians. 6371 corresponds to the radius of the earth in km.
// 
// The program will display the name of the defibrillator located the closest to the user’s position. This position is given as input to the program.
//  	Game Input
// 
// Input
// Line 1: User's longitude (in degrees)
// 
// Line 2: User's latitude (in degrees)
// 
// Line 3: The number N of defibrillators located in the streets of Montpellier
// 
// N next lines: a description of each defibrillator
// 
// Output
// The name of the defibrillator located the closest to the user’s position.
// Constraints
// 0 < N < 10000
// Example
// Input
// 3,879483
// 43,608177
// 3
// 1;Maison de la Prevention Sante;6 rue Maguelone 340000 Montpellier;;3,87952263361082;43,6071285339217
// 2;Hotel de Ville;1 place Georges Freche 34267 Montpellier;;3,89652239197876;43,5987299452849
// 3;Zoo de Lunaret;50 avenue Agropolis 34090 Mtp;;3,87388031141133;43,6395872778854
// Output
// Maison de la Prevention Sante

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

std::vector<std::string> string_split(const std::string &source, const char *delimiter, bool keepEmpty) {
        std::vector<std::string> results;
        size_t prev = 0;
        size_t next = 0;

        while((next = source.find_first_of(delimiter, prev)) != std::string::npos)
        {
            if(keepEmpty || (next - prev != 0)) {
                    results.push_back(source.substr(prev,next-prev));
            }
            prev = next + 1;
        }
        if(prev <= source.size())
        {
            results.push_back(source.substr(prev));
        }
        return results;
}

double euro_string_2_double(std::string str)
{
    auto found_comma = str.find_first_of(',');
    str.replace(found_comma, 1, ".");
    return std::stod(str);
}

struct Defibrillator
{
    Defibrillator(std::string data)
    {
        auto tokens = string_split(data, ";", true);

        Name = tokens[1];
        Address = tokens[2];
        lon = euro_string_2_double(tokens[4]);
        lat = euro_string_2_double(tokens[5]);
    }
    std::string Name;
    std::string Address;
    double lon;
    double lat;
};

double dist(double lat_a, double lon_a, double lat_b, double lon_b)
{
    double deg_to_rad = 3.1415926535897932384626433 / 180.0;
    lat_a *= deg_to_rad;
    lon_a *= deg_to_rad;
    lat_b *= deg_to_rad;
    lon_b *= deg_to_rad;
    double x = (lon_b - lon_a) * std::cos((lat_a + lat_b) / 2);
    double y = lat_b - lat_a;
    return std::sqrt(x*x + y*y) * 6371.0;
}

int main()
{
    std::cerr.precision(20);
    std::vector<Defibrillator> list;
    std::string LON;
    std::cin >> LON; std::cin.ignore();
    std::string LAT;
    std::cin >> LAT; std::cin.ignore();
    double my_lat = euro_string_2_double(LAT);
    double my_lon = euro_string_2_double(LON);
    std::cerr << "MY Pos: " << LAT << ", " << LON << std::endl;
    std::cerr << "MY Pos: " << my_lat << ", " << my_lon << std::endl;
    int N;
    std::cin >> N; std::cin.ignore();
    for (int i = 0; i < N; i++) {
        std::string DEFIB;
        getline(std::cin, DEFIB);
        std::cerr << DEFIB << std::endl;
        list.push_back(Defibrillator(DEFIB));
    }

    std::sort(list.begin(), list.end(), [&](const Defibrillator &lh, const Defibrillator &rh)
    {
        double dist_lh = dist(lh.lat, lh.lon, my_lat, my_lon);
        double dist_rh = dist(rh.lat, rh.lon, my_lat, my_lon);
        return dist_lh < dist_rh;
    });

    for (auto &d : list)
    {
        std::cerr << d.Name << " = " << dist(d.lat, d.lon, my_lat, my_lon) << std::endl;
    }

    // Write an action using cout. DON'T FORGET THE "<< endl"
    // To debug: cerr << "Debug messages..." << endl;

    std::cout << list.at(0).Name << std::endl;
}