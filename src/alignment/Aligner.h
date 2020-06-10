#ifndef __ALIGNER_H
#define __ALIGNER_H

#include <vector>
#include "../GeneList.h"
#include "../Profile.h"

class Aligner
{
public:
    /**
     * Align a vector of genelists
     * @param segments Vector of genelists to be aligned (input/output)
     */
    virtual void align(std::vector<GeneList*>& segments) = 0;

    /**
     * Virtual destructor
     */
    virtual ~Aligner() {}
};

#endif
