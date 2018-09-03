#ifndef CPPWATCH_H
#define CPPWATCH_H

#include <ctime>

class cppWatch{
 public:
  cppWatch();
  ~cppWatch(){}

  double totalInt;
  double currentInt;
  std::clock_t c_start;
  std::time_t t_start;

  void start(){c_start = std::clock(); t_start = std::time(NULL); return;};
  void stop();
  double total(){return totalInt;}
  double current(){return currentInt;}
  void clear(){totalInt = 0; currentInt = 0; return;}
};

cppWatch::cppWatch(){totalInt = 0; currentInt = 0; return;}

void cppWatch::stop()
{
  std::clock_t c_end = std::clock();
  //std::time_t t_end = std::time(NULL);

  //  totalInt += t_end  - t_start;
  //  currentInt = t_end  - t_start;
  totalInt += c_end - c_start; 
  currentInt = c_end - c_start;
  return;
}

#endif
