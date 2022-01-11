#include "AlignmentDrawer.h"
#include "ListElement.h"

using namespace std;

Conflict::Conflict(const Conflict& c)
{
    indexOfSeg1=c.indexOfSeg1;
    posInSeg1=c.posInSeg1;
    indexOfSeg2=c.indexOfSeg2;
    posInSeg2=c.posInSeg2;
}

Conflict& Conflict::operator=(const Conflict& c)
{
    indexOfSeg1=c.indexOfSeg1;
    posInSeg1=c.posInSeg1;
    indexOfSeg2=c.indexOfSeg2;
    posInSeg2=c.posInSeg2;

    return (*this);
}


AlignmentDrawer::AlignmentDrawer(const Profile* prof) : numColorsNeeded(0), visualizeConflicts(false), svg(NULL),
border(0.5), squareWidth(1.0),squareHeight(1.0),horSquareDist(0.5), vertSquareDist(1.0), badProfile(false)
{
    if (prof==NULL)
        cerr << "AlignmentDrawer: profile not initialized!" << endl;

    segments=prof->getSegments();
    profPtr=prof;

}

AlignmentDrawer::AlignmentDrawer(const vector<GeneList*>& segs) : numColorsNeeded(0), visualizeConflicts(false), svg(NULL),
border(0.5), squareWidth(1.0),squareHeight(1.0),horSquareDist(0.5), vertSquareDist(1.0), badProfile(false)
{
    profPtr=NULL;
    segments=segs;
}

AlignmentDrawer::~AlignmentDrawer()
{
    profPtr=NULL;
    colorMatrix.clear();
    conflicts.clear();
}

void AlignmentDrawer::initializeColorMatrix()
{
    assert(segments.size()!=0);
    colorMatrix.resize(segments.size());

    int maxSegSize=segments[0]->getSize();

    for (int i=0; i<segments.size(); i++) {
        if (segments[i]->getSize() != maxSegSize)
        {
            badProfile=true;
            int s=segments[i]->getSize();
            maxSegSize=(s>maxSegSize) ? s : maxSegSize;
        }
    }

    for (int i=0; i<colorMatrix.size(); i++) {
            colorMatrix[i].resize(maxSegSize);
    }

    if (!badProfile) {
        for (int i=0; i<colorMatrix.size(); i++) {
            vector<ListElement*>::const_iterator geneIt=segments[i]->getLEBegin();
            for (int j=0; geneIt<segments[i]->getLEEnd(); geneIt++) {
                assert(j<maxSegSize);
                const ListElement& listIJ=*(*geneIt);
                if (listIJ.isGap())
                    colorMatrix[i][j]=-1;
                else
                    colorMatrix[i][j]=0;
                j++;
            }
        }
    } else {
        for (int i=0; i<colorMatrix.size(); i++) {
            for (int j=0; j<colorMatrix[0].size(); j++) {
                colorMatrix[i][j]=0;
            }
        }

    }
}

bool AlignmentDrawer::buildColorMatrixOverlappingCode(DataSet* data, int tandemG, set< Link >::const_iterator start, set< Link >::const_iterator stop)
{
    initializeColorMatrix();

    if (badProfile){
        cout << "badprofile (all boxes will be black! => segmentlength differs among the segs of alignment)" << endl;
        return false;
    }

    set<Link>::const_iterator linkIt=start;

    int colorID=1;

    //(i,j) matrix position of first gene in link
    //(k,l) matrix position of second gene in link
    for (int i=0; i<colorMatrix.size(); ) {

        if (linkIt->segmentX > i) {
            i++; continue;
        }
        vector<ListElement*>::const_iterator geneIt=segments[i]->getLEBegin();

        for (int j=0; geneIt!=segments[i]->getLEEnd(); ) {

            const Gene& genIJ=(*geneIt)->getGene(); //gene at ith segment and jth place

            if (colorMatrix[i][j]<0) { //ignore gaps
                geneIt++; j++;
                continue;

            } else {

                const Gene& genLinkFirst=data->getGene(linkIt->geneXID); //first gene in link
                const Gene& genLInkSecond=data->getGene(linkIt->geneYID); //second gene in link

                //check if link is connected to gene at (i,j) in color matrix
                if (genLinkFirst==genIJ) {

                    int k=linkIt->segmentY;
                    int l=findGenePosInSegment(k, genLInkSecond);

                    if (j!=l)
                        if (abs(j-l) < tandemG) //ignore tandem duplicates
                            conflicts.push_back(Conflict(i,j,k,l));

                    if (colorMatrix[i][j]==0) {
                        if (colorMatrix[k][l]==0) {
                            //two blanco AP
                            colorMatrix[i][j]=colorID;
                            colorMatrix[k][l]=colorID;
                            colorID++;
                        } else {
                            //second AP is coloured
                            colorMatrix[i][j]=colorMatrix[k][l];
                        }
                    } else {
                        if (colorMatrix[k][l]==0) {
                            //first AP is coloured
                            colorMatrix[k][l]=colorMatrix[i][j];
                        } else {
                            //both are coloured
                            if (colorMatrix[i][j]!=colorMatrix[k][l]){
                                cerr << "Warning: Dubious indirect gene relationship between ";
                                cerr << genLinkFirst.getID() << " and " << genLInkSecond.getID();
                                cerr <<  ", closest genes get same color in alignment" << endl;


                                if (j==l){ //overwrite color otherwise do nothing (first link gets priority)
                                    colorMatrix[k][l]=colorMatrix[i][j];
                                }
                            }
                        }
                    }
                    linkIt++;
                    if (linkIt==stop) //terminate algorithm
                    {
                        numColorsNeeded=colorID-1;
                        return true;
                    }

                }
                else {
                    j++;
                    geneIt++;
                }
            }
        }
    }
    return true;
}


bool AlignmentDrawer::buildColorMatrixPostP(DataSet* data, int tandemG, set< Link >& homologs)
{
    return buildColorMatrixOverlappingCode(data,tandemG, homologs.begin(), homologs.end());
}



bool AlignmentDrawer::buildColorMatrix(DataSet* data, int tandemG)
{
    if (profPtr==NULL)
    {
        cerr << "ProfilePtr not initialized, you can use other buildColorMatrix if you provide the links" << endl;
        return false;
    }

    set<Link>::const_iterator start=profPtr->getHomBegin();
    set<Link>::const_iterator stop=profPtr->getHomEnd();
    return buildColorMatrixOverlappingCode(data,tandemG,start,stop);
}

int AlignmentDrawer::findGenePosInSegment(int indexSeg, const Gene& g)
{
    vector<ListElement*>::const_iterator geneFindIt=segments[indexSeg]->getLEBegin();
    int l=0; //index of second gene in segY
    for (; geneFindIt!=segments[indexSeg]->getLEEnd(); geneFindIt++) {

        if ((*geneFindIt)->getGene()==g) {

            break;
        }
        l++;
    }
    assert(geneFindIt!=segments[indexSeg]->getLEEnd()); //error if gene is not found!
    return l;
}

void AlignmentDrawer::generateAlignmentSVG(const std::string filename)
{
    double xDim,yDim;

    setBorder(); //shift border to have room for genome and listname
    setBoxSize(); //adapt box size to length of gene names

    //canvas dimensions
    yDim=border+(squareHeight+vertSquareDist)*colorMatrix.size();
    xDim=border+(squareWidth+horSquareDist)*colorMatrix[0].size();

    svg=new SvgWriter(filename,xDim,yDim);

    double coUpperCornerX=border;
    double coUpperCornerY=border;

    generateColorVector();

    for (int i=0; i<colorMatrix.size(); i++) {
        coUpperCornerX=border;
        for (int j=0; j<colorMatrix[i].size(); j++) {

            const int colorIndex=colorMatrix[i][j];
            switch (colorIndex) {

                case -1:
                    svg->drawRectangle(coUpperCornerX, coUpperCornerY, squareWidth, squareHeight,white);
                    break;
                case  0:
                    svg->drawRectangle(coUpperCornerX, coUpperCornerY, squareWidth, squareHeight,black);
                    break;
                default:
                    svg->drawRectangle(coUpperCornerX, coUpperCornerY, squareWidth, squareHeight,colors[colorIndex-1]);
                    //NOTE colorindex starts from 1 =>correction
                    break;
            }
            svg->drawDot(coUpperCornerX+squareWidth*0.5,coUpperCornerY+squareHeight*0.5,0.2*squareHeight);
            coUpperCornerX+=squareWidth+horSquareDist;
        }
        coUpperCornerY+=squareHeight+vertSquareDist;
    }

    //add gene names
    coUpperCornerX=border;
    coUpperCornerY=border;

    for (int i=0; i<colorMatrix.size(); i++) {
        coUpperCornerX=border;
        for (int j=0; j<colorMatrix[0].size(); j++) {

            const int colorIndex=colorMatrix[i][j];
            string geneName;
            switch (colorIndex) {

                case -1:
                    ;
                    break;
                case  0:
                    geneName=segments[i]->getLe(j).getGene().getID();
                    svg->drawText(coUpperCornerX, coUpperCornerY-0.1*squareHeight,geneName,6);
                    break;
                default:
                    geneName=segments[i]->getLe(j).getGene().getID();
                    svg->drawText(coUpperCornerX, coUpperCornerY-0.1*squareHeight,geneName,6);
                    break;
            }
            coUpperCornerX+=squareWidth+horSquareDist;
        }
        coUpperCornerY+=squareHeight+vertSquareDist;
    }

    //add genome and listnames
    coUpperCornerX=border;
    coUpperCornerY=border;

    for (int i=0; i<colorMatrix.size(); i++) {
        svg->drawText(0.5,coUpperCornerY+0.5*squareHeight,segments[i]->getGenomeName()+" "+segments[i]->getListName(),25);
        coUpperCornerY+=squareHeight+vertSquareDist;
    }


    if (visualizeConflicts) drawConflicts();
    svg->finish();
    delete svg;
}

void AlignmentDrawer::generateColorVector()
{
    Color tempColor;
    colors.resize(numColorsNeeded);
    for (int i=0; i<numColorsNeeded; i++) {
        tempColor.getRandomColor();
        colors[i]=tempColor;
    }
}

void AlignmentDrawer::drawConflicts()
{
    if (badProfile) {return;}
    if (svg!=NULL) {
        for (int i=0; i<conflicts.size(); i++) {
            const Conflict& conflict=conflicts[i];
            double startX,startY;
            double endX,endY;
            coordCenterAPij(conflict.indexOfSeg1,conflict.posInSeg1,startX,startY);
            coordCenterAPij(conflict.indexOfSeg2,conflict.posInSeg2,endX,endY);
            svg->drawLine(startX,startY,endX,endY);
        }
    }
    else cerr << "ERROR: SVGPTR not initialized!!" << endl;

}

void AlignmentDrawer::coordCenterAPij(int i, int j,double& x, double& y) const
{
    x=border+(squareWidth+horSquareDist)*j+squareWidth*0.5;
    y=border+(squareHeight+vertSquareDist)*i+squareHeight*0.5;
}

void AlignmentDrawer::setBorder()
{
    int maxChars=0;
    for (int i=0; i<segments.size(); i++){
        int charsInGenomeName=segments[i]->getGenomeName().size();
        int charsInListName=segments[i]->getListName().size();
        maxChars=max(charsInListName+1+charsInGenomeName,maxChars);
    }
    border+=(maxChars/5+1);
}

void AlignmentDrawer::setBoxSize()
{
    int maxChars=0;
    for (int i=0; i<segments.size(); i++){
        for (int j=0; j<segments[i]->getSize(); j++){
            int charsInGeneName=segments[i]->getLe(j).getGene().getID().length();
            maxChars=max(charsInGeneName,maxChars);
        }
    }
    squareHeight=maxChars*0.125;
    squareWidth=squareHeight;
}
