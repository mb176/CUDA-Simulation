#include <vector>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>


#include "vector_calculus.h"

using namespace std;

void add(const vector<double>& a,const vector<double>& b, vector<double>& c){
    if(a.size()==b.size() && b.size()==c.size()){
        for(int i =0; i!= a.size(); i++){
            c[i]=a[i]+b[i];
        }
    }else{
        throw domain_error("Vector sizes don't match" );
    }
    return;
}

void subtract(const vector<double>& a,const vector<double>& b, vector<double>& c){
    if(a.size()==b.size() && b.size()==c.size()){
        for(int i =0; i!= a.size(); i++){
            c[i]=a[i]-b[i];
        }
    }else{
        throw domain_error("Vector sizes don't match" );
    }
    return;
}

void mult(vector<double>& v, double a){
    for(vector<double>::iterator it = v.begin(); it!= v.end(); it++){
        *it = *it*a;
    }
    return;
}

double norm(vector<double>& v){
    double x= 0;
    for(vector<double>::iterator it = v.begin(); it!= v.end(); it++){
        x += (*it)*(*it);
        return sqrt(x);
    }
}

double norm2(vector<double>& v){
    double x= 0;
    for(vector<double>::iterator it = v.begin(); it!= v.end(); it++){
        x += (*it)*(*it);
        return x;
    }
}

void range(vector<double>& vec, double min, double max){
    uint size = vec.size();
    if(size==0)
        throw domain_error("Vector is empty");
    if(size == 1){
      vec[0] = (max - min)/2; //Give average if there is only one element
    } else {
      double intervall = (max-min)/(size-1);
      for(int i=0; i!= size; i++){
          vec[i] = min + i*intervall;
      }
    }
}
//Writes data into PATH in the form column[row] with " " as separation
ofstream write(const string PATH, const vector<vector<double> >& data){
    if(data.size()==0)
        throw domain_error("Data is empty");
    ofstream outfile(PATH);
    for(int vec_idx = 0; vec_idx!=data.size(); vec_idx++){//vec_idx vectors have been written in file
        for(int n = 0; n!= data[vec_idx].size(); n++){
            outfile << data[vec_idx][n] << " ";
        }
        outfile << endl;
    }
    return outfile;
}


ifstream read(const string PATH, vector<vector<double> >& data){
    if(data.size()!=0)
        cout << "Data is not empty, still proceeding";
    ifstream infile(PATH);
    string line;
    double value;
    vector<double> vec;
    while(getline(infile, line)){
        istringstream row(line);
        vec.clear();
        while(row>>value){
            vec.push_back(value);
        }
        data.push_back(vec);
    }
    return infile;
}

void print(const std::vector<double> vec){
    for(double i: vec){
        cout << i <<endl;
    }
}

void print(const std::vector<std::vector<double> > vec){
    for(vector<double> i: vec){
        for(double j : i)
        cout << j << " , ";
    }
    cout <<endl;
}
