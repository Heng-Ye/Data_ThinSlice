#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <vector>
#include <string>
#include <algorithm>    // copy
#include <iterator>     // ostream_operator
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <iostream>
#include <fstream>

//using namespace std;
using namespace boost;

void readdata() {

    string data("pos_cut.csv");

    ifstream in(data.c_str());
    if (!in.is_open()) return 1;

    typedef tokenizer< escaped_list_separator<char> > Tokenizer;
    vector< string > vec;
    string line;
    vector<int> RUN;
    vector<double> mu_Z;
    vector<double> sigma_Z;
    vector<double> mu_Y;
    vector<double> sigma_Y;
    vector<double> mu_X;
    vector<double> sigma_X;

    cout<<"Format: run#,mu_Zst,sigma_Zst,mu_Yst,sigma_Yst,mu_Xst,sigma_Xst"<<endl;

    while (getline(in,line)) {
        Tokenizer tok(line);
        vec.assign(tok.begin(),tok.end());

        // vector now contains strings from one row, output to cout here
        copy(vec.begin(), vec.end(), ostream_iterator<string>(cout, ","));
        //cout << "\n" << endl;
	//cout<<stoi(vec[0])<<" "<<stof(vec[1])<<" "<<stof(vec[2])<<endl;
	RUN.push_back(stoi(vec[0]));

	mu_Z.push_back(stof(vec[1]));
	double tmp_sigma_Z=stof(vec[2]); if (tmp_sigma_Z<0) tmp_sigma_Z=-tmp_sigma_Z;
	sigma_Z.push_back(tmp_sigma_Z);

	mu_Y.push_back(stof(vec[3]));
	double tmp_sigma_Y=stof(vec[4]); if (tmp_sigma_Y<0) tmp_sigma_Y=-tmp_sigma_Y;
	sigma_Y.push_back(tmp_sigma_Y);

	mu_X.push_back(stof(vec[5])); 
	double tmp_sigma_X=stof(vec[6]); if (tmp_sigma_X<0) tmp_sigma_X=-tmp_sigma_X;
	sigma_X.push_back(tmp_sigma_X);
        cout << "\n----------------------" << endl;
    }
    
     	
    for (size_t j=0; j<RUN.size(); ++j) {
	cout<<"RUN:"<<RUN.at(j)<<" "<<mu_Z.at(j)<<" "<<sigma_Z.at(j)<<" "<<mu_Y.at(j)<<" "<<sigma_Y.at(j)<<" "<<mu_X.at(j)<<" "<<sigma_X.at(j)<<endl;
    }

}
