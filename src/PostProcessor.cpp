#include "PostProcessor.h"
#include "GeneList.h"
#include "ListElement.h"
#include "AlignmentDrawer.h"


PostProcessor::PostProcessor(int mID, const string& path, DataSet* data, bool useFamily) : multipliconID(mID), outputPath(path)
{
    segmentInfo=getSegmentInfo();
    vector<vector<listElInfo> > listElementsInfo=getListElements(segmentInfo);

    dataset=data;
    useFamilies=useFamily;

    int maxPos=0;
    for (int i=0; i<listElementsInfo.size(); i++){
        int indexLastListEl=listElementsInfo[i].size()-1;
        maxPos=max(listElementsInfo[i].at(indexLastListEl).pos,maxPos);
    }

    //cout << "Maxpos= " << maxPos << endl;

    for (int i=0; i<segmentInfo.size(); i++){
        createSegment(segmentInfo[i],listElementsInfo[i],maxPos);
    }

    createLinks();
}

PostProcessor::~PostProcessor()
{
    for (int i=0; i<segments.size(); i++){
        delete segments[i];
    }
    segments.clear();
}

vector<segInfo> PostProcessor::getSegmentInfo()
{
    vector<segInfo> res;
    string filename=outputPath+segmentFile;
    ifstream ifs(filename.c_str(),ifstream::in);
    int maxLineLength=250;
    char buffer[maxLineLength];
    string str; 

    bool stop=false;
    bool firstFound=false;

    int segId=-1;
    int mID=-1;
    string gName;
    string lName;
    string firstEl;
    string lastEl;

    do {

        getline(ifs,str);
        stringstream ss(str);

        ss >> segId;
        ss >> mID;

        if (mID==multipliconID) {

            ss >> gName;
            ss >> lName;
            ss >> firstEl;
            ss >> lastEl;

            segInfo s;
            s.id=segId;
            s.genomeName=gName;
            s.listName=lName;
            s.firstGene=firstEl;
            s.lastGene=lastEl;
            res.push_back(s);

            firstFound=true;
        } else {

            if (firstFound) stop=true;
        }

        segId=-1; //reset
        mID=-1;

    }  while (!ifs.eof() and !stop);

    ifs.close();
    return res;
}

vector<vector<listElInfo> > PostProcessor::getListElements(const vector<segInfo>& sInfo)
{
    vector<vector<listElInfo> >  res;
    int nSegments=sInfo.size();

    for (int i=0; i<nSegments; i++){
        vector<listElInfo> temp;
        res.push_back(temp);
    }

    string filename=outputPath+listElementsFile;
    int sIDFirst=sInfo[0].id;
    int sIDLast=sInfo[nSegments-1].id;

    ifstream ifs(filename.c_str(),ifstream::in);
    int maxLineLength=250;
    char buffer[maxLineLength];
    string str;

    int segId=-1;
    int lEId=-1;
    string geneName;
    int genePos=-1;
    bool firstFound=false; bool stop=false;

    do {

        getline(ifs,str);
        stringstream ss(str);

        ss >> lEId;
        ss >> segId;

        if (segId>=sIDFirst and segId<=sIDLast) {

            ss >> geneName;
            ss >> genePos;

            assert(genePos>=0);
            assert(geneName.length()>0);

            listElInfo l;
            l.gene=geneName;
            l.pos=genePos;

            res[segId-sIDFirst].push_back(l);

            firstFound=true;
        } else {

            if (firstFound) stop=true;
        }

        segId=-1; //reset
        lEId=-1;
        geneName.clear();
        genePos=-1;

    }  while (!ifs.eof() and !stop);

    ifs.close();

    return res;
}

void PostProcessor::createSegment(segInfo& sI, vector< listElInfo >& lEs, int maxPos)
{
    GeneList* geneList=dataset->getGeneList(sI.listName,sI.genomeName);

    if (geneList==NULL){
        cerr << "Genelist not found" << endl;
        return;
    }

    //create <genename, listEl map>
    map<string,ListElement*> geneListElMap;
    for (int i=0; i<geneList->getSize(); i++){
        ListElement& element=geneList->getLe(i);
        string key=element.getGene().getID();
        geneListElMap.insert(pair<string,ListElement*>(key,&element));
    }

    vector<ListElement*> assembledListEls; //assembled list of listelements without gaps
    //but correct order as in aligned seg
    for (int i=0; i<lEs.size(); i++){
        map<string,ListElement*>::const_iterator itEl=geneListElMap.find(lEs[i].gene);
        assert(itEl!=geneListElMap.end());
        assembledListEls.push_back(itEl->second);
    }

    GeneList* segment;
    //copy the segment into a genelist
    segment= new GeneList(geneList->getListName(),geneList->getGenomeName(), assembledListEls);

    //introduce the gaps made in the alignment
    int lEsSize=lEs.size();
    int posPreviousElement=-1;
    for (int i=0; i<lEsSize; i++){
        int currentPos=lEs[i].pos;
        int nGaps=currentPos-posPreviousElement-1;

        for (int j=0; j<nGaps; j++)
            segment->introduceGap(posPreviousElement+1);

        posPreviousElement=currentPos;
    }
    int nGapsAtEnd=maxPos-lEs[lEsSize-1].pos;
    for (int j=0; j<nGapsAtEnd; j++)
            segment->pasteGap();

    segments.push_back(segment);
}

void PostProcessor::createLinks()
{
    if (!useFamilies) {
        for (int i=1; i<segments.size(); i++){
            for (int j=0; j<segments[i]->getSize(); j++){ //iterate over all listelements in ith segment
                //cout << "segY "<< i << " " << j << endl;
                ListElement& lEY=segments[i]->getLe(j);
                if (lEY.isGap())
                    continue;
                if (!lEY.getGene().hasPairs())
                    continue;

                for (int k=0; k<i; k++){ //iterate over all previous segments
                    for (int l=0; l<segments[k]->getSize(); l++){
                       // cout << "segX "<< k << " " << l << endl;
                        ListElement& lEX=segments[k]->getLe(l);
                        if (lEX.isGap())
                            continue;
                        if (!lEX.getGene().isIndirectPairWith(lEY.getGene()))
                            continue;
                        // store the homolog gene IDs
                        homologs.insert(Link(k,i,lEX.getNumID(), lEY.getNumID()));
                       // cout << "link inserted" << endl;
                        lEX.hasHomolog = true;
                        lEY.hasHomolog = true;
                    }
                }
            }

        }
    } else {

        map<string, vector<pair<ListElement*, int> > > pairMap;

        //insert all listElements in pairMap
        for (int i=0; i<segments.size(); i++){
            for (int j=0; j<segments[i]->getSize(); j++){
                ListElement &le = segments[i]->getLe(j);
                if (le.isGap())
                    continue;
                pairMap[le.getGene().getFamily()].push_back(pair<ListElement*, int>(&le, i));

            }
        }

        for (int i=1; i<segments.size(); i++) {
            for (int j=0; j<segments[i]->getSize(); j++){ //iterate over all listelements in ith segment
                ListElement& lEY=segments[i]->getLe(j);
                if (lEY.isGap())
                    continue;

                map<string, vector<pair<ListElement*, int> > >::iterator it=pairMap.find(lEY.getGene().getFamily());
                vector<pair<ListElement*, int> > &group = it->second;

                //cout << "familySize=" << group.size() << endl;

                for (int k=0; k<group.size(); k++) {
                    ListElement &lEX = *group[k].first;
                    int segX = group[k].second;

                    if (segX>=i) //only links with previous segments!
                        continue;

                    // to have compatibility with the pairwise case:
                    if (lEX.getNumID() == lEY.getNumID())
                        continue;

                    homologs.insert(Link(segX, i,
                                        lEX.getNumID(), lEY.getNumID()));

                    //cout << "linkInfo " << lEX.getGene().getID() << " <-> " << segments[i]->getLe(j).getGene().getID() << endl;

                    lEX.hasHomolog = true;
                    lEY.hasHomolog = true;
                }
            }
        }
    }
}

void PostProcessor::visualizeMultiplicon(int tandemG)
{
    AlignmentDrawer drawer(segments);

    drawer.buildColorMatrixPostP(dataset, tandemG,homologs);
    drawer.visualizeConflicts=true;

    string tempfilename=outputPath+"AlignmentMultiplicon";
    char buffer[50];
    sprintf(buffer,"%i.svg",multipliconID);
    string filename=tempfilename+string(buffer);
    drawer.generateAlignmentSVG(filename);
}

void PostProcessor::printMultiplicon()
{
    for (int i=0; i<segments.size(); i++){
        for (int j=0; j<segments[i]->getSize(); j++){
            if (segments[i]->getLe(j).isGap()) cout << "g";
            else cout << segments[i]->getLe(j).getGene().getID();
            cout << " ";
        }
        cout << endl;
    }
}
