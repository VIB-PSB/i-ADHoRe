#include "Settings.h"

Settings::Settings(const string& settingsFile)
        : gap_size(0), cluster_gap(0),cloud_gap_size(0),cloud_cluster_gap(0),
        tandem_gap(0), q_value(0.0), anchorpoints(0), prob_cutoff(0.0), level_2_only(false),
        use_family(false), alignment_method(NeedlemanWunsch), nThreads(1),
        mulHypCor(Bonferroni), compareAligners(false), max_gaps_in_alignment(0),
        flush_output(1000), clusterType(Collinear),visualizeGHM(false),cloudFiltermethod(Binomial),
        visualizeAlignment(false), verbose_output(true), bruteForceSynthenyMode(false)
{
    string genomename, listname, filename;

    ifstream fin (settingsFile.c_str());

    // any errors opening the file?
    if (!fin) throw FileException ("Error opening the settings file: " +
                                       settingsFile);

    string buffer;

    write_statistics = false;

    while (fin.good()) {
        readLine(fin, buffer);
        if (buffer.length() == 0) continue;

        unsigned int next;

        if (startsWith(buffer, "compareAligners", next)) {
            buffer.erase(0, next);
            compareAligners = true;
        }
        else if (startsWith(buffer, "number_of_threads", next)) {
            nThreads = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "genome", next)) {
            buffer.erase(0, next);
            readFromBuffer(genomename, buffer);
        }
        else if (startsWith(buffer, "blast_table", next)) {
            buffer.erase(0, next);
            readFromBuffer(blast_table, buffer);
        }
        else if (startsWith(buffer, "output_path", next)) {
            buffer.erase(0, next);
            readFromBuffer(output_path, buffer);
            // add a backslash to the directory if necessary
            if (!output_path.empty())
                if (output_path[output_path.length() - 1] != '/')
                    output_path.append("/");

        }
        else if (startsWith(buffer, "flush_output", next)) {
            flush_output = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "gap_size", next)) {
            gap_size = atoi(&buffer[next]);
        }
        else if (startsWith(buffer,"cloud_gap_size",next)){
            cloud_gap_size=atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "cluster_gap", next)) {
            cluster_gap = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "cloud_cluster_gap", next)) {
            cloud_cluster_gap = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "max_gaps_in_alignment", next)) {
            max_gaps_in_alignment = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "tandem_gap", next)) {
            tandem_gap = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "q_value", next)) {
            q_value = atof(&buffer[next]);
        }
        else if (startsWith(buffer, "anchor_points", next)) {
            anchorpoints = atoi(&buffer[next]);
        }
        else if (startsWith(buffer, "prob_cutoff", next)) {
            prob_cutoff = atof(&buffer[next]);
        }
        else if (startsWith(buffer, "verbose_output", next)) {
            buffer.erase(0, next);
            string boolean;
            readFromBuffer(boolean, buffer);
            if (boolean == "true") {
                verbose_output = true;
            }
            else if (boolean == "false") {
                verbose_output = false;
            }
            else {
                throw FileException ("ERROR: verbose_output "
                                     "attribute is not correct (should be \"true\""
                                     " or \"false\")");
            }
        }
        else if (startsWith(buffer, "level_2_only", next)) {
            buffer.erase(0, next);
            string boolean;
            readFromBuffer(boolean, buffer);
            if (boolean == "true")
                level_2_only = true;
            else if (boolean == "false")
                level_2_only = false;
            else
                throw FileException ("ERROR: The level_2_only attribute is not "
                                     "correct (should be \"true\" or \"false\")");
        }
        else if (startsWith(buffer, "write_stats", next)) {
            buffer.erase(0, next);
            string boolean;
            readFromBuffer(boolean, buffer);
            if (boolean == "true") {
                write_statistics = true;
            }
            else if (boolean == "false") {
                write_statistics = false;
            }
            else {
                throw FileException ("ERROR: The write_stats "
                                     "attribute is not correct (should be \"true\""
                                     " or \"false\")");
            }
        }
        else if (startsWith(buffer, "table_type", next)) {
            buffer.erase(0, next);
            string type;
            readFromBuffer(type, buffer);
            if (type == "family") {
                use_family = true;
            }
            else if (type == "pairs") {
                use_family = false;
            }
            else {
                throw FileException ("ERROR: The blast_table_type "
                                     "attribute is not correct (should be \"family\""
                                     " or \"pairs\")");
            }
        }
        else if (startsWith(buffer, "alignment_method", next)) {
            buffer.erase(0, next);
            string alignment_str;
            readFromBuffer(alignment_str, buffer);
            if (alignment_str == "nw") {
                alignment_method = NeedlemanWunsch;
            }
            else if (alignment_str == "gg") {
                alignment_method = GreedyGraphbased;
            }
            else if (alignment_str == "gg2") {
                alignment_method = GreedyGraphbased4;
            }
            //NOTE old aligners not available in release of i-ADHoRE 3.0
            /*else if (alignment_str == "gg3") {
                alignment_method = GreedyGraphbased3;
            }*/
            else if (alignment_str == "gg4") {
                alignment_method = GreedyGraphbased4;
            }
            else {
                throw FileException ("ERROR: Alignment method "
                                     "specified is not correct (should be one of "
                                     "these: \"nw\", \"gg2\")");
            }
        }
        else if (startsWith(buffer, "multiple_hypothesis_correction", next)) {
            buffer.erase(0, next);
            string multCorStr;
            readFromBuffer(multCorStr, buffer);
            if (multCorStr == "none") {
                mulHypCor = None;
            }
            else if (multCorStr == "bonferroni") {
                mulHypCor = Bonferroni;
            }
            else if (multCorStr == "FDR") {
                mulHypCor = FDR;
            }
            else {
                throw FileException ("ERROR: Multiple hypothesis "
                                     "testing specified is not correct (should be one of "
                                     "these: \"none\", \"bonferroni\", \"FDR\")");
            }
        }

        else if (startsWith(buffer, "cluster_type", next)) {

            buffer.erase(0, next);
            string clustertype_str;
            readFromBuffer(clustertype_str, buffer);
            if (clustertype_str == "colinear" or clustertype_str == "collinear") //backward compatibility (previous cases used colLinear)
                clusterType=Collinear;
            else if (clustertype_str == "cloud") 
                clusterType=Cloud;
            else if (clustertype_str == "hybrid") 
                clusterType=Hybrid;
            else {
                throw FileException ("ERROR: Cluster type "
                                     "specified is not correct (should be \"collinear\" or \"cloud\" or \"hybrid\")");
            }
        }

        else if (startsWith(buffer, "visualizeGHM", next))
        {

            buffer.erase(0, next);
            string vis_str;
            readFromBuffer(vis_str, buffer);
            if (vis_str == "true")
                visualizeGHM=true;
            else if (vis_str == "false")
                visualizeGHM=false;
            else
                throw FileException ("ERROR: visualizeGHM should be 'true' or 'false'");
        }

        else if (startsWith(buffer, "visualizeAlignment", next))
        {

            buffer.erase(0, next);
            string vis_str;
            readFromBuffer(vis_str, buffer);
            if (vis_str == "true")
                visualizeAlignment=true;
            else if (vis_str == "false")
                visualizeAlignment=false;
            else
                throw FileException ("ERROR: visualizeAlignment should be 'true' or 'false'");
        }

        else if (startsWith(buffer, "cloud_filter_method", next))
        {

            buffer.erase(0, next);
            string cfm_str;
            readFromBuffer(cfm_str, buffer);
            if (cfm_str == "binomial") cloudFiltermethod=Binomial;
            else if (cfm_str == "binomial_corr") cloudFiltermethod=BinomialCorr;
            else {
                throw FileException ("ERROR: cloud_filter_method should be 'binomial', 'binomial_corr' or 'density'");
            }
        }

        else if (startsWith(buffer,"visGPairs",next)) {
            buffer.erase(0, next);
            string pairString;
            readFromBuffer(pairString, buffer);
            readPairsFromString(pairString);
        }

        else if (startsWith(buffer,"bruteForceSynthenyMode",next)) {
            buffer.erase(0, next);
            string bfString;
            readFromBuffer(bfString, buffer);
            if (bfString=="on")
                bruteForceSynthenyMode=true;
            else if (bfString=="off")
                bruteForceSynthenyMode=false;
            else {
                throw FileException ("ERROR: bruteForceSynthenyMode should be 'on' or 'off'");
            }
        }

        else {
            unsigned int i = 0;
            //reading listName
            while (i < buffer.length() && buffer[i] != ' ') {
                listname.push_back(buffer[i]);
                i++;
            }
            while (i < buffer.length() && buffer[i] == ' ') {
                i++;
            }
            //reading fileName
            while (i < buffer.length() && buffer[i] != '\n') {
                filename.push_back(buffer[i]);
                i++;
            }

            if (filename.empty()){
                cout << "Problematic listname=" << listname << endl;
                throw FileException ("ERROR: No Filename next to listname! (Also check in ini file for invalid parameters)");
            }

            listfiles.push_back(ListFile (genomename, listname, filename));

            listname.clear();
            filename.clear();
        }
    }

    fin.close();

    //check whether the necessary options are filled in
    if (genomename.empty()) {
        throw FileException ("ERROR: Genome name not found in settings file");
    }
    if (listfiles.empty()) {
        throw FileException ("ERROR: Genelist files not found in settings file");
    }
    if (output_path.empty()) {
        throw FileException ("ERROR: Output path not found in settings file");
    }
    if (blast_table.empty()) {
        throw FileException ("ERROR: BLAST table not found in settings file");
    }

    if (clusterType!=Cloud){
        if (gap_size <= 0) {
    //for picoplaza it was requested that this setting could be set at 1
            throw FileException ("ERROR: Gap size not correct in settings file");
        }
        if (cluster_gap <= 0 || cluster_gap < gap_size) {
            throw FileException ("ERROR: Cluster gap size not correct in settings"
                                    " file (should be >= gap_size)");
        }
    }
    if (clusterType!=Collinear){
        level_2_only=true;
        if (cloud_gap_size <= 1) {
            throw FileException ("ERROR: Cloud Gap size not correct in settings file");
        }
        if (cloud_cluster_gap <= 1 || cloud_cluster_gap < cloud_gap_size) {
            throw FileException ("ERROR: Cloud Cluster gap size not correct in settings"
                                " file (should be >= cloud_gap_size)");
        }
    }

    if (max_gaps_in_alignment < 1) {
        cerr << "WARNING: Maximum allowed number of gaps in the alignment "
             "not specified.  Setting to cluster_gap." << endl;
        max_gaps_in_alignment = cluster_gap;
    }
    if (tandem_gap <= 1) {
        // This raises an error, but continious using default value
        // for max compatibility with older ini files
        cerr << "WARNING: Tandem gap size not correct in settings file. "
             "Using default (gap_size / 2)" << endl;
        tandem_gap = gap_size / 2;
    }
    if (q_value <= 0.0) {
        throw FileException ("q_value should be > 0.0");
    }
    if (anchorpoints <= 2) {
        throw FileException ("ERROR: Minimum number of anchorpoints not "
                             "correct in settings file (should be >= 3)");
    }
    if (prob_cutoff <= 0.0) {
        throw FileException ("ERROR: Probability cutoff not correct in "
                             "settings file (should be > 0)");
    }

    if (max_gaps_in_alignment < cluster_gap) {
        cerr << "WARNING: max_gaps_in_alignment is smaller than "
                "cluster_gap.  This may reject L2 multiplicons." << endl;
    }

    if (max_gaps_in_alignment < gap_size) {
        cerr << "WARNING: max_gaps_in_alignment is smaller than "
                "gap_size.  This may reject L2 multiplicons." << endl;
    }

    if ((clusterType != Collinear) && (!level_2_only))
        throw FileException ("ERROR: For Synthenic Cloud Search and Hybrid search no "
                             "profiles have been implemented! => level_2_only should be true!");
}

void Settings::displaySettings() const
{
    cout << endl;
    cout << "************* i-ADHoRe parameters *************" << endl;

    cout << "\tNumber of genelists = "     << listfiles.size() << endl;

    cout << "\tBlast table = "             << blast_table             << endl;
    cout << "\tOutput path = "             << output_path             << endl;
    cout << "\tGap size = "                << gap_size                << endl;
    cout << "\tCluster gap size = "        << cluster_gap             << endl;
    cout << "\tCloud gap size = "          << cloud_gap_size          << endl;
    cout << "\tCloud cluster gap size = "  << cloud_cluster_gap       << endl;

    cout << "\tMax gaps in alignment = "   << max_gaps_in_alignment   << endl;
    cout << "\tTandem gap = "              << tandem_gap              << endl;
    cout << "\tFlush output = "            << flush_output            << endl;
    cout << "\tQ-value = "                 << q_value                 << endl;
    cout << "\tAnchorpoints = "            << anchorpoints            << endl;
    cout << "\tProbability cutoff = "      << prob_cutoff             << endl;

    cout <<  "\tCloud filtering method = ";
    switch (cloudFiltermethod) {
        case Binomial:
            cout << "Binomial";
            break;
        case BinomialCorr:
            cout << "BinomialCorr";
            break;

    }
    cout << endl;

    cout << "\tLevel 2 only = ";
    if (level_2_only)
        cout << "true" << endl;
    else
        cout << "false" << endl;

    cout << "\tUse family = ";
    if (use_family)
        cout << "true" << endl;
    else
        cout << "false" << endl;

    cout << "\tWrite statistics = ";
    if (write_statistics)
        cout << "true" << endl;
    else
        cout << "false" << endl;

    cout << "\tAlignment method = ";
    switch (alignment_method){
        case NeedlemanWunsch:
            cout << "NeedlemanWunsch";
            break;
        case GreedyGraphbased:
            cout << "GreedyGraphbased";
            break;
        case GreedyGraphbased2:
            cout << "GreedyGraphbased2";
            break;
        case GreedyGraphbased3:
            cout << "GreedyGraphbased3";
            break;
        case GreedyGraphbased4:
            cout << "GreedyGraphbased4";
            break;
    }
    cout << endl;

    cout << "\tMultiple hypothesis correction = ";
    switch (mulHypCor){
        case None:
            cout << "None";
            break;
        case Bonferroni:
            cout << "Bonferroni";
            break;
        case FDR:
            cout << "FDR";
            break;
    }
    cout << endl;

    cout << "\tNumber of threads = " << nThreads << endl;

    cout << "\tCompare aligners = ";
    if (compareAligners) cout << "true" << endl;
    else cout << "false" << endl;

    switch(clusterType){
        case Collinear:
            cout << "\tCollinear searches only" << endl;
            break;
        case Cloud:
            cout << "\tSynthenic cloud searches only" << endl;
            break;
        case Hybrid:
            cout << "\tCombined search (first collinear, then cloud)" << endl;
            break;
    }

    #ifdef HAVE_PNG
    cout << "\tVisualize GHM.png = ";
    #else
    cout << "\tVisualize GHM.bmp = ";
    #endif
    if (visualizeGHM)
        cout << "true" << endl;
    else
        cout << "false" << endl;

    cout << "\tVisualize Alignment = ";
    if (visualizeAlignment)
        cout << "true" << endl;
    else
        cout << "false" << endl;
    cout << "\tVerbose output = ";
    if (verboseOutput())
        cout << "true" << endl;
    else
        cout << "false" << endl;

    cout << "************ END i-AdDHoRe parameters *********" << endl;
    cout << endl;
}

bool Settings::startsWith(const string& buffer, const string& start, unsigned int& next) {
    unsigned int i = 0;
    while (i < start.length() && i < buffer.length() && buffer[i] == start[i]) {
        i++;
    }
    if (i == buffer.length()) {
        next = i;
        return true;
    }
    else if (i == start.length()) {
        while (buffer[i] != '=') {
            i++;
        }
        i++;
        while (buffer[i] == ' ') {
            i++;
        }
        next = i;
        return true;
    }
    else {
        next = 0;
        return false;
    }
}

void Settings::readFromBuffer(char* target, const char* buffer, const int maxtarget, const int maxbuffer) {
    int i = 0;
    while (i < maxtarget && i < maxbuffer && buffer[i] != '\n') {
        target[i] = buffer[i];
        i++;
    }
    if (i < maxtarget)
        target[i] = 0;
}

void Settings::readFromBuffer(string& target, const string& buffer) {
    target.clear();
    unsigned int i = 0;
    while (i < buffer.length() && buffer[i] != '\n') {
        target.push_back(buffer[i]);
        i++;
    }
}

void Settings::readLine(ifstream& fin, char* buffer, const int maxbuffer) {
    char c = fin.get();
    while (fin.good() && (c == ' ' || c == '\n')) {
        c = fin.get();
    }
    if (fin.good()) {
        buffer[0] = c;
        c = fin.get();
        int i = 1;
        while (fin.good() && i < maxbuffer && c != '\n') {
            buffer[i] = c;
            c = fin.get();
            i++;
        }
        if (i != maxbuffer) {
            buffer[i] = 0;
        }
    }
    else {
        buffer[0] = 0;
    }
}

void Settings::readLine(ifstream& fin, string& buffer) {
    buffer.clear();
    char c = fin.get();
    while (fin.good() && (c == ' ' || c == '\n')) {
        c = fin.get();
    }
    if (fin.good()) {
        buffer.push_back(c);
        c = fin.get();
        int i = 1;
        while (fin.good() && c != '\n') {
            buffer.push_back(c);
            c = fin.get();
            i++;
        }
    }
}

void Settings::getGapSizes(int* gapsizes) const {
    int nr_gaps = 0;
    double step = (log(gap_size) / log(3) - 1 ) / 9;
    for (int i = 0; i < 10; i++) {
        int value = (int) round(pow(3,(1 + (i * step))));
        if (i == 0 || value != gapsizes[nr_gaps-1]) {
            gapsizes[nr_gaps] = value;
            nr_gaps++;
        }
    }

    //end of the gaps in the array
    if (nr_gaps < 10) gapsizes[nr_gaps] = -100;
}

void Settings::getCloudGapSizes(int* cloudgapsizes) const {
    int nr_gaps = 0;
    double step = (log(cloud_gap_size) / log(3) - 1 ) / 9;
    for (int i = 0; i < 10; i++) {
        int value = (int) round(pow(3,(1 + (i * step))));
        if (i == 0 || value != cloudgapsizes[nr_gaps-1]) {
            cloudgapsizes[nr_gaps] = value;
            nr_gaps++;
        }
    }

    //end of the gaps in the array
    if (nr_gaps < 10) cloudgapsizes[nr_gaps] = -100;
}

bool Settings::showGHM(int listX, int listY) const
{
    int list1, list2;

    if (listX<listY){
        list1=listX; list2=listY;
    } else {
        list2=listX; list1=listY;
    }

    if (visualizeGHM)
        return true;
    else {

        map<int, set<int> >::const_iterator itMap=GHMPairsToVisualize.find(list1);

        if (itMap!=GHMPairsToVisualize.end()){

            set<int> pairsWithX=itMap->second;
            set<int>::const_iterator itSet=pairsWithX.find(list2);

            if (itSet!=pairsWithX.end())
                return true;
        }
    }
    return false;
}



/*void Settings::readIDsFromString(const string& s) //NOTE can be removed
{
    int ind1=0;
    int ind2=s.find_first_of(" ",ind1);
    bool success=true;

    while (ind2!=string::npos){
        stringstream ss(s.substr(ind1,ind2-ind1),stringstream::in);
        int id;

        if (ss>>id)
            ;
        else
            success = false;

        mIDsToVisualize.insert(id);

        ind1=ind2+1;
        ind2=s.find_first_of(" ",ind1);


    }
    if (!success) throw FileException ("ERROR: One of the MultiplionIDs in the input files appears not to be an integer!");
}*/


void Settings::readPairsFromString(const string& pairs)
{
    //NOTE pairlist has following layout: list1a genome1a li1b, pair2a pair2b, etc
    //genelists
    int indStart=0;
    int indSpace1=pairs.find_first_of(" ",indStart);
    int indSpace2=pairs.find_first_of(" ",indSpace1+1);
    int indSpace3=pairs.find_first_of(" ",indSpace2+1);
    int indComma=pairs.find_first_of(",",indSpace3+1);

    while (indStart!=string::npos ){
        string firstList=pairs.substr(indStart,indSpace1-indStart);
        string firstGenome=pairs.substr(indSpace1+1,indSpace2-indSpace1-1);
        string secondList=pairs.substr(indSpace2+1,indSpace3-indSpace2-1);
        string secondGenome;
        if (indComma==string::npos)
            secondGenome=pairs.substr(indSpace3+1);
        else
            secondGenome=pairs.substr(indSpace3+1,indComma-indSpace3-1);

        if (secondList.empty() or firstList.empty()) {
            throw FileException ("ERROR: One of the GHM listnames in GHMpairsToVisualize is empty!");
        }
        else if (secondGenome.empty() or firstGenome.empty()) {
            throw FileException ("ERROR: One of the GHM genomenames in GHMpairsToVisualize is empty!");
        }
        int indFirst=findGeneListName(firstList,firstGenome);
        int indSecond=findGeneListName(secondList,secondGenome);

        if (indFirst<indSecond){

            if (GHMPairsToVisualize.find(indFirst)!=GHMPairsToVisualize.end())
                GHMPairsToVisualize[indFirst].insert(indSecond);
            else {
                set<int> firstEl;
                firstEl.insert(indSecond);
                GHMPairsToVisualize.insert(pair<int,set<int> >(indFirst, firstEl));
            }
        } else {

            if (GHMPairsToVisualize.find(indSecond)!=GHMPairsToVisualize.end())
                GHMPairsToVisualize[indSecond].insert(indFirst);
            else {
                set<int> firstEl;
                firstEl.insert(indFirst);
                GHMPairsToVisualize.insert(pair<int,set<int> >(indSecond, firstEl));
            }
        }

        indStart=(indComma!=string::npos)? indComma+2: string::npos;
        indSpace1=pairs.find_first_of(" ",indStart);
        indSpace2=pairs.find_first_of(" ",indSpace1+1);
        indSpace3=pairs.find_first_of(" ",indSpace2+1);
        indComma =pairs.find_first_of(",",indSpace3+1);
    }
}

int Settings::findGeneListName(string lname, string gname) const
{
    list<ListFile>::const_iterator it=listfiles.begin();
    int index=0;
    for (; it!=listfiles.end(); it++) {
        string listname=it->getListName();
        string genomename=it->getGenomeName();
        if ( (lname.compare(listname) ==0)
            and (gname.compare(genomename)==0))
                return index;

        index++;
    }

    throw FileException ("ERROR: Genelist not found!");
    return -1;
}
