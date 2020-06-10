#ifndef __SETTINGS_H
#define __SETTINGS_H

//#define DEBUG

#include "ListFile.h"


#include "debug/FileException.h"
#include <sstream>
#include <string>
#include <iostream>
#include "headers.h"

enum AlignmentMethod {
    NeedlemanWunsch,
    GreedyGraphbased,
    GreedyGraphbased2,
    GreedyGraphbased3,
    GreedyGraphbased4
};

enum MultHypCor {
    None,
    Bonferroni,
    FDR
};

enum ClusterType 
{
    Collinear,
    Cloud,
    Hybrid
};

enum FilterMethod
{
    Binomial,
    BinomialCorr,
};


class Settings {

public:

    ///////////////
    //CONSTRUCTOR//
    ///////////////

    /*
    *constructs the settings object and reads every property from the specified settings file
    */
    Settings(const string& settingsFile);


    //////////////////
    //PUBLIC METHODS//
    //////////////////

    /**
    * Write all relevant parameters to the output buffer
    */
    void displaySettings() const;

    /*
    *returns a list with listfiles read from the settings file
    */
    const list<ListFile>& getListFiles() const {
        return listfiles;
    }

    /*
    *calculates the 10 or fewer integers ranging from 3 to gap size
    */
    void getGapSizes(int* gaps) const;
    void getCloudGapSizes(int* cloudgapsizes) const;

    /*
    *returns the gap size
    */
    int getGapSize() const {
        return gap_size;
    }

    int getCloudGapSize() const{
        return cloud_gap_size;
    }

    int getMaxGapsInAlignment() const {
        return max_gaps_in_alignment;
    }

    MultHypCor getMultHypCorMethod() const {
        return mulHypCor;
    }

    /*
    *returns the cluster gap size
    */
    int getClusterGap() const {
        return cluster_gap;
    }

    int getCloudClusterGap() const{
        return cloud_cluster_gap;
    }

    /*
    *returns the tandem gap size
    */
    int getTandemGap() const {
        return tandem_gap;
    }

    /*
    *returns the q-value
    */
    double getQValue() const {
        return q_value;
    }

    /*
    *returns the number of anchorpoints
    */
    int getAnchorPoints() const {
        return anchorpoints;
    }

    /*
    *returns the probability cutoff
    */
    double getProbCutoff() const {
        return prob_cutoff;
    }

    /*
    *returns a C-string of the blast table filename
    */
    const string& getBlastTable() const {
        return blast_table;
    }

    /*
    *returns a C-string of the path where the output needs to be
    */
    const string& getOutputPath() const {
        return output_path;
    }

    /*
    *returns true if the user wants only level 2 multiplicons calculated
    */
    bool level2Only() const {
        return level_2_only;
    }

    /*
    *returns true if the user uses blast pairs or families
    */
    bool useFamily() const {
        return use_family;
    }

    /*
    *returns true if the user wants to calculate duplicated and collinear portions
    */
    bool writeStatistics() const {
        return write_statistics;
    }
    /*
    *returns the alignment method specified in the settings file, default is NeedlemanWunsch
    */
    AlignmentMethod getAlignmentMethod() const {
        return alignment_method;
    }

    /**
     * Get the number of threads
     */
    int getNumThreads() const {
        return nThreads;
    }

    bool getCompareAligners() const {
        return compareAligners;
    }

    int getFlushOutput() const {
        return flush_output;
    }

    ClusterType getClusterType() const
    {
        return clusterType;
    }

    FilterMethod getCloudFilterMethod() const
    {
        return cloudFiltermethod;
    }

    bool showGHM(int listX, int listY) const;

    bool showAlignedProfiles() const
    {
        return visualizeAlignment;
    }

    bool verboseOutput() const
    {
        return verbose_output;
    }

    bool isBruteforce() const
    {
        return bruteForceSynthenyMode;
    }

private:
    ///////////////////
    //PRIVATE METHODS//
    ///////////////////

    /*
    *returns a boolean saying the specified string 'buffer' starts the same
    *as the string 'start'
    */
    //bool startsWith(const char* buffer, const char* start, int maxbuffer, int maxstart, int& next);
    bool startsWith(const string& buffer, const string& start, unsigned int& next);

    /*
    *reads the string from the buffer into the target
    */
    void readFromBuffer(char* target, const char* buffer, const int maxtarget, const int maxbuffer);
    void readFromBuffer(string& target, const string& buffer);

    /*
    *reads a line from the settings file and is placed into the buffer
    */
    void readLine(ifstream& fin, char* buffer, const int maxbuffer);
    void readLine(ifstream& fin, string& buffer);


    void readPairsFromString(const string& pairs);

    int findGeneListName(string lname, string gname) const;
    //////////////
    //ATTRIBUTES//
    //////////////

    //static const int MAX = 50;
    //static const int MAXBUFFER = 256;

    list<ListFile> listfiles;
    string blast_table;
    string output_path;
    int gap_size;
    int cluster_gap;
    int max_gaps_in_alignment;
    int tandem_gap;
    int flush_output;
    double q_value;
    int anchorpoints;
    double prob_cutoff;
    bool level_2_only;
    bool use_family;
    bool write_statistics;
    AlignmentMethod alignment_method;
    MultHypCor mulHypCor;
    int nThreads;
    bool compareAligners;

    //cloud related parameters
    int cloud_gap_size;
    int cloud_cluster_gap;
    ClusterType clusterType;
    bool visualizeGHM;
    bool visualizeAlignment;
    FilterMethod cloudFiltermethod;
    bool verbose_output;
    bool bruteForceSynthenyMode;

    map<int, set<int> > GHMPairsToVisualize;

};

#endif
