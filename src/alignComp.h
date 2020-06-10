#ifndef ALIGN_COMP
#define ALIGN_COMP

#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>

typedef unsigned int uint;

class AlignScore
{
    public:
        /**
         * Add the partial results of a certain alignment
         * @param level Level of the generated alignment
         * @param numAP Number of initial anchorpoints in the alignment
         * @param numAlAP Number of aligned anchorpoints in the alignment
         * @param numHom Number of initial homologs in the alignment
         * @param numAlHom Number of aligned homologs in the alignment
         * @param lengthAl Lenght of the alignment
         * @param time Time it took to do the alignment
         */
        void addScore(uint level, uint numAP, uint numAlAP, uint numHom,
                      uint numAlHom, uint lengthAl, double time);

        /**
         * Reset all information
         */
        void reset();

        /**
         * Return the number of levels for which we have information
         */
        uint getNumLevels() {
            return numAP.size();
        }

        /**
         * Return the number of AP for a certain level
         * @param level Level
         */
        uint getNumAP(uint level) {
            return level < numAP.size() ? numAP.at(level) : 0;
        }

        /**
         * Return the total number of AP over all levels
         * @return The total number of AP
         */
        uint getTotalNumAP() {
            return std::accumulate(numAP.begin(), numAP.end(), 0u);
        }

        /**
         * Return the number of aligned AP for a certain level
         * @param level Level
         */
        uint getNumAlAP(uint level) {
            return level < numAlAP.size() ? numAlAP.at(level) : 0;
        }

        /**
         * Return the total number of aligned AP over all levels
         * @return The total number of aligned AP
         */
        uint getTotalNumAlAP() {
            return std::accumulate(numAlAP.begin(), numAlAP.end(), 0u);
        }

        /**
         * Return the number of homologs for a certain level
         * @param level Level
         */
        uint getNumHom(uint level) {
            return level < numHom.size() ? numHom.at(level) : 0;
        }

        /**
         * Return the total number of homologs over all levels
         * @return The total number of homologs
         */
        uint getTotalNumHom() {
            return std::accumulate(numHom.begin(), numHom.end(), 0u);
        }

        /**
         * Return the number of aligned homologs for a certain level
         * @param level Level
         */
        uint getNumAlHom(uint level) {
            return level < numAlHom.size() ? numAlHom.at(level) : 0;
        }

        /**
         * Return the total number of aligned homologs over all levels
         * @return The total number of aligned homologs
         */
        uint getTotalNumAlHom() {
            return std::accumulate(numAlHom.begin(), numAlHom.end(), 0u);
        }

        /**
         * Return the total length of profiles for a certain level
         * @param level Level
         */
        uint getLengthAl(uint level) {
            return level < lengthAl.size() ? lengthAl.at(level) : 0;
        }

        /**
         * Return the total number of profiles over all levels
         * @return The total length of the profiles of
         */
        uint getTotalLengthAl() {
            return std::accumulate(lengthAl.begin(), lengthAl.end(), 0u);
        }

        /**
         * Return the number of profiles in a certain level
         * @param level Level
         */
        uint getNumProfiles(uint level) {
            return level < numProfiles.size() ? numProfiles.at(level) : 0;
        }

        /**
         * Return the total number of profiles over all levels
         * @return The total number of profiles
         */
        uint getTotalNumProfiles() {
            return std::accumulate(numProfiles.begin(), numProfiles.end(), 0u);
        }

        /**
         * Return the alignment time to generate the profiles at a certain level
         * @param level Level
         */
        double getTime(uint level) {
            return level < timeProfiles.size() ? timeProfiles.at(level) : 0;
        }

        /**
         * Return the total alignment time over all profile levels
         * @return The total alignment time
         */
        double getTotalTime() {
            return std::accumulate(timeProfiles.begin(),
                                   timeProfiles.end(), 0.0);
        }

    private:
        std::vector<uint> numAP;
        std::vector<uint> numAlAP;
        std::vector<uint> numHom;
        std::vector<uint> numAlHom;
        std::vector<uint> lengthAl;
        std::vector<uint> numProfiles;
        std::vector<double> timeProfiles;
};

#endif
