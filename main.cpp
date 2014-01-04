#include "pmtree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <sys/time.h>

class Timer
{
    timeval timer[2];
    public:
    timeval start()
    {
        gettimeofday(&this->timer[0], NULL);
        return this->timer[0];
    }
    timeval stop()
    {
        gettimeofday(&this->timer[1], NULL);
        return this->timer[1];
    }
    int duration() const
    {
        int secs(this->timer[1].tv_sec - this->timer[0].tv_sec);
        int usecs(this->timer[1].tv_usec - this->timer[0].tv_usec);
        if(usecs < 0)
        {
            --secs;
            usecs += 1000000;
        }
        return static_cast<int>(secs * 1000 + usecs / 1000.0 + 0.5);
    }
};

int main() {
    /*
       vector<double> times {306, 455, 1010, 210, 883, 1022, 310, 361, 218, 166,
       170, 654, 728i, 71, 567, 144, 613, 707, 61, 88 };
       vector<int> status {1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
       vector<vector<double> > adjustments;
       vector<vector<double> > treatments{{0,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,0}};
       vector<vector<double> > splits{{1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,0,0,0}};

*/
    vector<double> times;
    vector<int> status;
    vector<vector<double> > adjustments;
    vector<vector<double> > treatments;
    vector<vector<double> > splits;


    string line;

    ifstream infile("times.txt");
    while(getline(infile, line))
    {
        istringstream iss(line);
        double i;
        while( iss >> i ) 
            times.push_back(i);
    }
    /*
       for(vector<double>::iterator it=times.begin(); it!=times.end();it++){
       cout<<*it<<endl;
       }

       cout<<times.size()<<endl;
       */

    ifstream infile2("status.txt");
    while(getline(infile2, line))
    {
        stringstream ss(line);
        int i;
        while( ss >> i ) 
            status.push_back(i);
    }

    ifstream infile3("adjustments.txt");
    while(getline(infile3, line))
    {
        stringstream ss(line);
        double i;
        vector<double> v;
        while( ss >> i ) 
            v.push_back(i);
        adjustments.push_back(v);
    }
    ifstream infile4("treatments.txt");

    while(getline(infile4, line))
    {
        stringstream ss(line);
        double i;
        vector<double> v;
        while( ss >> i ) 
            v.push_back(i);
        treatments.push_back(v);
    }

    ifstream infile5("splits.txt");

    while(getline(infile5, line))
    {
        stringstream ss(line);
        double i;
        vector<double> v;
        while( ss >> i ) 
            v.push_back(i);
        splits.push_back(v);
    }
    //print dimensions of inputs. so far so good.
    /*cout<<times.size()<<" "<<status.size()<<" "<<adjustments.size()<<" "<<
      adjustments[0].size()<<" "<<treatments.size()<<" "
      <<treatments[0].size()<< " "<<splits.size()<<" "<<splits[0].size()<<
      endl<<flush; 
      */ 
Timer tm;

    tree pmtree(times, status, adjustments, treatments, splits);
tm.start();
    pmtree.build_tree(10, 10);
tm.stop();
pmtree.print_tree(true);
pmtree.print_tree(false);    
cout<<tm.duration()<<endl;
    /*
       pmtree.print_tree(true);
       pmtree.print_tree(false);
       pmtree.build_tree(2, 2);
       pmtree.print_tree(true);
       pmtree.print_tree(false);
       pmtree.build_tree(1, 1);
       pmtree.print_tree(true);
       pmtree.print_tree(false);
       */
    return 0;
}
