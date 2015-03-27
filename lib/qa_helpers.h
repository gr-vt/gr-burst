#ifndef INCLUDED_QA_HELPERS_H
#define INCLUDED_QA_HELPERS_H

#include <gnuradio/block.h>
#include <vector>
#include <string>

// TODO: need to remove this hardcoded stuff
//static std::string testFilesDirectory = "/home/kiran/awst/pybombs/src/gr-ieee-80211b/testfiles/";
static std::string testFilesDirectory = "/tmp/";

class qa_helpers {
public:

	// helper functions to run tests
    static std::vector<gr_complex> readComplexFile(std::string filename);
    static std::vector<float> readFloatFile(std::string filename);
    static std::vector<int> readIntFile(std::string filename);
    static std::vector<unsigned char> readUCharFile(std::string filename);
    static std::vector< std::vector<gr_complex> > readComplexFileAsMat(std::string filename, int M, int N);
    static std::vector< std::vector<float> > readFloatFileAsMat(std::string filename, int M, int N);

    static void writeComplexFile(std::string filename, std::vector<gr_complex> vec);
    static void writeFloatFile(std::string filename, std::vector<float> vec);
    static void writeIntFile(std::string filename, std::vector<int> vec);
    static void writeUCharFile(std::string filename, std::vector<unsigned char> vec);

    static long getFileSize(FILE *file);
    static std::vector<gr_complex> readComplexBinFile(std::string filename);

    static std::vector<std::string> getFilesInDir(const std::string& pat);
    static bool fileExists(const char *filename);

    static void areComplexVectorsEqual(std::vector<gr_complex> expectedVec, std::vector<gr_complex> actualVec, double delta);
    static void areFloatVectorsEqual(std::vector<float> expectedVec, std::vector<float> actualVec, double delta);
    static void areIntVectorsEqual(std::vector<int> expectedVec, std::vector<int> actualVec);
    static void areUCharVectorsEqual(std::vector<unsigned char> expectedVec, std::vector<unsigned char> actualVec);
};

#endif
