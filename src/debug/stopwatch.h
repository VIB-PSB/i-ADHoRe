#ifndef __STOPWATCH
#define __STOPWATCH
#include <time.h>

class Stopwatch {
    public:
    	clock_t begin, end;
    	Stopwatch () {
        	start ();
    	}
    	void start (){
      	  begin = clock();
    	}
    	void stop (){
     	   end = clock();
    	}
    	double time () {
        	return (double) (end-begin) / CLOCKS_PER_SEC;
    	}
};

#endif
