#include <cmath>
#include <vector>
#include <stdexcept>

#include "fitting.h"

using namespace std;

//Fits f(x)=a*x+b to the data, returns {a,b,error}
vector<double> linear_fit(const vector<double>& xdata, const vector<double>& ydata){
    double size = xdata.size();
    if(size!= ydata.size())
        throw domain_error("x and y data doesn't match");
    if(size <= 1)
        throw domain_error("To few datapoints");

    double sumx(0),sumxx(0),sumxy(0),sumy(0);
    for(int i =0; i!= size; i++){//Invariant: i datapoints have been evaulated;
        sumx += xdata[i];
        sumxx += xdata[i]*xdata[i];
        sumxy += xdata[i]*ydata[i];
        sumy += ydata[i];
    }

    //Parameter for f(x)=a*x+b
    double a = (size*sumxy - sumx *sumy)/(size*sumxx-sumx*sumx); //Parameter for f(x)=a*x+b
    double b = (sumxx*sumy-sumx*sumxy)/(size*sumxx-sumx*sumx);

    //Fehler berechnen
    double sum(0),error(0);
    for(int i=0; i!= size; i++){
        error += pow(ydata[i]-a*xdata[i]+b,2);
        sum += ydata[i]*ydata[i];
    }
    error *=1/sum;
    vector<double> output={a,b,error};
    return output;
}

//Fits f(x)=b*x^a to the data, returns {a,ln(b),error}
vector<double> exponent_fit(const vector<double>& xdata, const vector<double>& ydata){
    double size = xdata.size();
    if(size!= ydata.size())
        throw domain_error("x and y data doesn't match");
    if(size <= 1)
        throw domain_error("To few datapoints");
    vector<double> logXdata, logYdata;
    for(int i =0; i!=size; i++){
        logXdata.push_back(log(xdata[i]));
        logYdata.push_back(log(ydata[i]));
    }
    return linear_fit(logXdata,logYdata);
}
