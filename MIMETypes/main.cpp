// MIME types are used in numerous internet protocols to associate a media type (html, image, video ...) with the content sent. The MIME type is generally inferred from the extension of the file to be sent.
// 
// You have to write a program that makes it possible to detect the MIME type of a file based on its name.
//  	Rules
// You are provided with a table which associates MIME types to file extensions. You are also given a list of names of files to be transferred and for each one of these files, you must find the MIME type to be used.
// 
// The extension of a file is defined as the substring which follows the last occurrence, if any, of the dot character within the file name.
// If the extension for a given file can be found in the association table (case insensitive, e.g. TXT is treated the same way as txt), then print the corresponding MIME type. If it is not possible to find the MIME type corresponding to a file, or if the file doesn’t have an extension, print UNKNOWN.
//  	Game Input
// Input
// Line 1: Number N of elements which make up the association table.
// 
// Line 2: Number Q of file names to be analyzed.
// 
// N following lines: One file extension per line and the corresponding MIME type (separated by a blank space).
// 
// Q following lines: One file name per line.
// 
// Output
// For each of the Q filenames, display on a line the corresponding MIME type. If there is no corresponding type, then display UNKNOWN.
// Constraints
// 0 < N < 10000 
// 0 < Q < 10000
// File extensions are composed of a maximum of 10 alphanumerical ASCII characters.
// MIME types are composed of a maximum 50 alphanumerical and punctuation ASCII characters.
// File names are composed of a maximum of 256 alphanumerical ASCII characters and dots (full stops).
// There are no spaces in the file names, extensions or MIME types.
// Example
// Input
// 2
// 4
// html text/html
// png image/png
// test.html
// noextension
// portrait.png
// doc.TXT
// Output
// text/html
// UNKNOWN
// image/png
// UNKNOWN

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <locale>

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/
int main()
{
    std::map<std::string, std::string> mime_tbl;

    int N; // Number of elements which make up the association table.
    std::cin >> N; std::cin.ignore();
    int Q; // Number Q of file names to be analyzed.
    std::cin >> Q; std::cin.ignore();
    for (int i = 0; i < N; i++) {
        std::string EXT; // file extension
        std::string MT; // MIME type.
        std::cin >> EXT >> MT; std::cin.ignore();
        
        std::locale loc;
        for(auto & c : EXT)
        {
            c = std::tolower(c, loc);
        }
        
        std::cerr << EXT << " " << MT << std::endl;        
        mime_tbl.emplace(std::make_pair(EXT, MT));
    }
    for (int i = 0; i < Q; i++) {
        std::string FNAME; // One file name per line.
        getline(std::cin, FNAME);

        std::string extension("");
        auto pos = FNAME.find_last_of('.');        
        if(pos != std::string::npos)
        {
            extension = FNAME.substr(pos + 1);
            std::locale loc;
            for(auto & c : extension)
            {
                c = std::tolower(c, loc);
            }
        }

        std::cerr << FNAME << " ______ " << extension << "  " << pos << std::endl;
        auto found_ext = mime_tbl.find(extension);
        if (found_ext != mime_tbl.end())
        {
            std::cout << found_ext->second << std::endl;;
        }
        else
        {
            std::cout << "UNKNOWN" << std::endl;
        }
    }
}