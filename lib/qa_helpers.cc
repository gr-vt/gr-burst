#include "qa_helpers.h"
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include <glob.h>
#include <string>

#include <iostream>
#include <stdio.h>

std::vector<gr_complex> qa_helpers::readComplexFile(std::string filename) {
	std::vector<gr_complex> complexData;

	std::ifstream infile(filename.c_str());
	float realPart, imagPart;
	while(infile>>realPart) {
		infile >> imagPart;
		gr_complex num;
		num.real(realPart);
		num.imag(imagPart);
		complexData.push_back(num);
	}

	return complexData;
}

std::vector< std::vector<gr_complex> > qa_helpers::readComplexFileAsMat(std::string filename, int M, int N) {
	std::vector< std::vector<gr_complex> > cmplxMat;

	std::ifstream infile(filename.c_str());
	float realPart, imagPart;;
	for(int ii=0; ii<M; ii++) {
		std::vector<gr_complex> row(N, gr_complex(0,0));
		for(int jj=0; jj<N; jj++) {
			infile >> realPart;
			infile >> imagPart;
			gr_complex num;
			num.real(realPart);
			num.imag(imagPart);
			row[jj] = num;
		}
		cmplxMat.push_back(row);
	}

	return cmplxMat;
}

std::vector<float> qa_helpers::readFloatFile(std::string filename) {
	std::vector<float> floatData;

	std::ifstream infile(filename.c_str());
	float num;
	while(infile>>num) {
		floatData.push_back(num);
	}

	return floatData;
}

std::vector< std::vector<float> > qa_helpers::readFloatFileAsMat(std::string filename, int M, int N) {
	std::vector< std::vector<float> > floatMat;

	std::ifstream infile(filename.c_str());
	float num;
	for(int ii=0; ii<M; ii++) {
		std::vector<float> row(N,0.0);
		for(int jj=0; jj<N; jj++) {
			infile>>num;
			row[jj]=num;
		}
		floatMat.push_back(row);
	}

	return floatMat;
}

std::vector<int> qa_helpers::readIntFile(std::string filename) {
	std::vector<int> intData;

	std::ifstream infile(filename.c_str());
	int num;
	while(infile>>num) {
		intData.push_back(num);
	}

	return intData;
}

std::vector<unsigned char> qa_helpers::readUCharFile(std::string filename) {
	std::vector<unsigned char> ucharData;

	std::ifstream infile(filename.c_str());
	unsigned char num;
	while(infile>>num) {
		ucharData.push_back(num-'0');
	}

	return ucharData;
}

void qa_helpers::areComplexVectorsEqual(std::vector<gr_complex> expectedVec, std::vector<gr_complex> actualVec, double delta) {
//	std::cout << "expectedVec.size() = " << expectedVec.size() << "\n"; fflush(stdout);
//	std::cout << "actualVec.size() = " << actualVec.size() << "\n"; fflush(stdout);
	if(expectedVec.size()!=actualVec.size()) {
		CPPUNIT_FAIL("Expected Vector Size and Actual Vector Size are different!");
	}
	for(int ii=0; ii<expectedVec.size(); ii++) {
		float real_expected = expectedVec[ii].real();
		float imag_expected = expectedVec[ii].imag();
		float real_actual = actualVec[ii].real();
		float imag_actual = actualVec[ii].imag();

		if( (real_actual > (real_expected+delta)) || (real_actual < (real_expected-delta)) ||
			(imag_actual > (imag_expected+delta)) || (imag_actual < (imag_expected-delta)) ) {
			std::cout << "real_expected=" << real_expected << " real_actual=" << real_actual
					 << " imag_expected=" << imag_expected << " imag_actual=" << imag_actual << " " << ii << "\n";
			fflush(stdout);
		}

		CPPUNIT_ASSERT_DOUBLES_EQUAL(real_expected, real_actual, delta);
		CPPUNIT_ASSERT_DOUBLES_EQUAL(imag_expected, imag_actual, delta);

	}
}

void qa_helpers::areFloatVectorsEqual(std::vector<float> expectedVec, std::vector<float> actualVec, double delta) {
//	std::cout << "expectedVec.size() = " << expectedVec.size() << "\n"; fflush(stdout);
//	std::cout << "actualVec.size() = " << actualVec.size() << "\n"; fflush(stdout);
	if(expectedVec.size()!=actualVec.size()) {
		CPPUNIT_FAIL("Expected Vector Size and Actual Vector Size are different!");
	}
	for(int ii=0; ii<expectedVec.size(); ii++) {
		float expectedVal = expectedVec[ii];
		float actualVal = actualVec[ii];

		if( actualVal > (expectedVal+delta) || actualVal < (expectedVal-delta) ) {
			std::cout << "expected=" << expectedVal << " actual=" << actualVal << " " << ii << "\n"; fflush(stdout);
		}

		CPPUNIT_ASSERT_DOUBLES_EQUAL(expectedVal, actualVal, delta);
	}
}

void qa_helpers::areIntVectorsEqual(std::vector<int> expectedVec, std::vector<int> actualVec) {
	if(expectedVec.size()!=actualVec.size()) {
		CPPUNIT_FAIL("Expected Vector Size and Actual Vector Size are different!");
	}
	for(int ii=0; ii<expectedVec.size(); ii++) {
		int expectedVal = expectedVec[ii];
		int actualVal = actualVec[ii];

//    		std::cout << "expectedVal=" << expectedVal << " actualVal=" << actualVal << "\n";

		CPPUNIT_ASSERT_EQUAL(expectedVal, actualVal);
	}
}

void qa_helpers::areUCharVectorsEqual(std::vector<unsigned char> expectedVec, std::vector<unsigned char> actualVec) {
	if(expectedVec.size()!=actualVec.size()) {
		CPPUNIT_FAIL("Expected Vector Size and Actual Vector Size are different!");
	}
	for(int ii=0; ii<expectedVec.size(); ii++) {
		int expectedVal = expectedVec[ii];
		int actualVal = actualVec[ii];

//    	std::cout << "expectedVal=" << expectedVal << " actualVal=" << actualVal << "\n";

		CPPUNIT_ASSERT_EQUAL(expectedVal, actualVal);
	}
}

std::vector<std::string> qa_helpers::getFilesInDir(const std::string& pat) {
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    std::vector<std::string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(std::string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}

bool qa_helpers::fileExists(const char *filename) {
  std::ifstream ifile(filename);
  return ifile.good();
}

// Get the size of a file
long qa_helpers::getFileSize(FILE *file) {
    long lCurPos, lEndPos;
    lCurPos = ftell(file);
    fseek(file, 0, 2);
    lEndPos = ftell(file);
    fseek(file, lCurPos, 0);
    return lEndPos;
}

std::vector<gr_complex> qa_helpers::readComplexBinFile(std::string filename) {
	std::vector<gr_complex> binData;
	FILE *file = NULL;
	const char* filePath = filename.c_str();
	if ((file = fopen(filePath, "rb")) == NULL) {
		std::cout << "Could not open specified file" << std::endl;
		return binData;
	}

	long fileSize = qa_helpers::getFileSize(file);
	long numElem = fileSize/sizeof(gr_complex);
	binData.resize(numElem);
	fread(&binData[0], sizeof(gr_complex), numElem, file);
	fclose(file);

	return binData;
}

void qa_helpers::writeComplexFile(std::string filename, std::vector<gr_complex> vec) {
	std::ofstream fout;
	std::string fname = filename;
	fout.open(fname.c_str());
	for(int ii=0; ii<vec.size(); ii++) {
		fout << vec[ii].real() << "\n" << vec[ii].imag() << "\n";
	}
	fout.close();
}

void qa_helpers::writeFloatFile(std::string filename, std::vector<float> vec) {
	std::ofstream fout;
	std::string fname = filename;
	fout.open(fname.c_str());
	for(int ii=0; ii<vec.size(); ii++) {
		fout << vec[ii] << "\n";
	}
	fout.close();
}

void qa_helpers::writeIntFile(std::string filename, std::vector<int> vec) {
	std::ofstream fout;
	std::string fname = filename;
	fout.open(fname.c_str());
	for(int ii=0; ii<vec.size(); ii++) {
		fout << vec[ii] << "\n";
	}
	fout.close();
}

void qa_helpers::writeUCharFile(std::string filename, std::vector<unsigned char> vec) {
	std::ofstream fout;
	std::string fname = filename;
	fout.open(fname.c_str());
	for(int ii=0; ii<vec.size(); ii++) {
		fout << (int)vec[ii] << "\n";
	}
	fout.close();
}
