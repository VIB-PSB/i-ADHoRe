#include "GHMProfile.h"

#include "AnchorPoint.h"
#include "BaseCluster.h"
#include "Profile.h"
#include "ListElement.h"
#include "Gene.h"
#include "Multiplicon.h"

#include <cassert>

GHMProfile::GHMProfile(const Profile& xObject, const GeneList& yObject)
: GHM(*xObject.getSegments()[0], yObject), x_object(xObject) {
    level = xObject.getLevel() + 1;
}

void GHMProfile::buildMatrix ()
{
    // reset the number of points in the GHM
    count_points[OPP_ORIENT] = 0;
    count_points[SAME_ORIENT] = 0;

    // reset the GHM datastructure
    matrix.clear();
    matrix.resize(2);

    for (unsigned int i = 0; i < x_object.getSegments().size(); i++) {

        const vector<ListElement*>& xList = x_object.getSegments()[i]->getRemappedElements();
        const vector<ListElement*>& yList = y_object.getRemappedElements();

        for (unsigned int x = 0; x < xList.size(); x++) {
            const ListElement &xElement = *xList[x];

            if (xElement.isGap()) continue;
            if (!xElement.getGene().hasPairs()) continue;

            std::queue<int> queue;
            xElement.matchingPositions(yList, queue);

            while (!queue.empty()) {
                // get and remove first element
                unsigned int y = queue.front();
                queue.pop();

                const ListElement &yElement = *yList[y];

                Orientation orient = (xElement.getOrientation() ==
                                      yElement.getOrientation()) ?
                    SAME_ORIENT : OPP_ORIENT;

                if (matrix[orient][x].find(y) == matrix[orient][x].end()) {
                    count_points[orient]++;
                    matrix[orient][x].insert(y);
                }
            }
        }
    }
}

void GHMProfile::getMultiplicons(vector<Multiplicon*>& mps) const
{
    for (unsigned int i = 0; i < multiplicons.size(); i++) {
        mps.push_back(multiplicons[i]);
    }
}

void GHMProfile::setMultiplicons(bool useFamily, int minHomologs)
{
    vector<const vector<ListElement*>* > x_elements;
    for (unsigned int i = 0; i < x_object.getSegments().size(); i++) {
        x_elements.push_back(&x_object.getSegments()[i]->getRemappedElements());
    }

    const vector<ListElement*>& y_elements = y_object.getRemappedElements();

    vector<Multiplicon*> multipliconsCopy;
    multipliconsCopy.reserve(multiplicons.size());

    for (unsigned int i = 0; i < multiplicons.size(); i++) {
        multiplicons[i]->setBounds();

        const vector<BaseCluster*>& baseclusters = multiplicons[i]->getBaseClusters();
        for (unsigned int j = 0; j < baseclusters.size(); j++) {
            baseclusters[j]->setBounds();
            baseclusters[j]->filterDuplicates();

            vector<AnchorPoint> false_anchorpoints;

            multiset<AnchorPoint>::const_iterator e = baseclusters[j]->getAPBegin();
            for ( ; e != baseclusters[j]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();
                const Gene &gY = y_elements[y]->getGene();

                bool first_match_encountered = false;
                for (unsigned int n = 0; n < x_elements.size(); n++) {
                    const vector<ListElement*>* list = x_elements[n];
                    const ListElement *le = list->at(x);
                    if (list->at(x)->isGap()) continue;
                    const Gene &gX = le->getGene();

                    if (!gX.isPairWith(gY)) continue;

                    if (first_match_encountered) {
                        false_anchorpoints.push_back(AnchorPoint(le->getNumID(), y_elements[y]->getNumID(), x, y, false));
                    } else {
                        // this const_cast is allowed, because setting the gene IDs
                        // does not invalidate the multiset iterator
                        const_cast<AnchorPoint&>(*e).setGeneIDs(le->getNumID(),
                                                                y_elements[y]->getNumID());
                        first_match_encountered = true;
                    }
                }

                if (!first_match_encountered) {
                    cerr << "\t\tAnchorpoint without matching x-gene!!!" << endl;
                    exit(1);
                }
            }
            for (unsigned int k = 0; k < false_anchorpoints.size(); k++) {
                baseclusters[j]->addAnchorPoint(false_anchorpoints[k]);
            }
        }

        multiplicons[i]->extractXObject(x_object);
        multiplicons[i]->extractYObject(y_object);

        if (multiplicons[i]->createHomologs(useFamily, minHomologs))
            multipliconsCopy.push_back(multiplicons[i]);
        else
            delete multiplicons[i];
    }

    multiplicons = multipliconsCopy;
}


#ifdef HAVE_PNG
void GHMProfile::visualizeGHM(const std::string& output_path) const
{
    char lev[20];
    sprintf (lev, "_Level%d", x_object.getLevel());
    char multip[20];
    sprintf(multip, "x=%d..%d",x_object.getMultiplicon().getLowestX(),x_object.getMultiplicon().getHighestX());
    string filename = "GHMProfile"+string(lev)+"_"+string(multip)+"_"+y_object.getListName()+".png";
    cout << "Visualize: " << filename << endl;

    map<int, set<int> >::iterator itX;
    set<int>::const_iterator itY;

    const vector<ListElement*>& yList = y_object.getRemappedElements();

    Grafix png(x_object.getSize() + 4 - (x_object.getSize() % 4), yList.size());

    for (int orient = 0; orient < 2; orient++) {
        const map<int, set<int> > &mat = matrix[orient];
        map<int, set<int> >::const_iterator itX = mat.begin();
        map<int, set<int> >::const_iterator endX = mat.end();

        png.setDrawingColor(white);
        // plot the dots which have never been in any cluster
        for ( ; itX != endX; itX++) {
            int x = itX->first;
            const set<int> setY = itX->second;

            set<int>::const_iterator itY = setY.begin();
            set<int>::const_iterator endY = setY.end();

            for ( ; itY != endY; itY++) {
                int y = *itY;
                png.putPixel(x, y);
            }
        }

        vector<BaseCluster*>::const_iterator it;
        it = filteredBC[orient].begin();
        while (it != filteredBC[orient].end()) {
            png.setDrawingColor(red);
            int xmin=(*it)->getLowestX();
            int xmax=(*it)->getHighestX();
            int ymin=(*it)->getLowestY();
            int ymax=(*it)->getHighestY();
            png.drawBox(xmin,xmax,ymin,ymax);

            png.setDrawingColor(white);
            multiset<AnchorPoint>::const_iterator e = (*it)->getAPBegin();
            for ( ; e != (*it)->getAPEnd(); e++) {
                    int x = e->getX();
                    int y = e->getY();
                    png.putPixel(x, y);
                }
            it++;
        }
    }

    vector<Multiplicon*>::const_iterator it = multiplicons.begin();
    for (int j=0; j<multiplicons.size(); j++) {
        vector<BaseCluster*> BCs=multiplicons[j]->getBaseClusters();

        for (int i=0; i<BCs.size(); i++){
            png.setDrawingColor(green);
            int xmin=BCs[i]->getLowestX();
            int xmax=BCs[i]->getHighestX();
            int ymin=BCs[i]->getLowestY();
            int ymax=BCs[i]->getHighestY();
            png.drawBox(xmin,xmax,ymin,ymax);

            png.setDrawingColor(blue);
            multiset<AnchorPoint>::const_iterator f = BCs[i]->getAPBegin();
            double upL,downL;
            double upR,downR;
            int xL,xR;
            xL=f->getX();
            BCs[i]->intervalBounds(xL,upL,downL);
            f++;
            for ( ; f != BCs[i]->getAPEnd(); f++) {
                xR=f->getX();
                BCs[i]->intervalBounds(xR,upR,downR);
                png.drawCircle(xR,upR,1.0);
                png.drawCircle(xR,downR,1.0);
                upL=upR;
                downL=downR;
                xL=xR;
            }

            multiset<AnchorPoint>::const_iterator e = BCs[i]->getAPBegin();
            png.setDrawingColor(yellow);
            for ( ; e != BCs[i]->getAPEnd(); e++) {
                int x = e->getX();
                int y = e->getY();

                png.drawPixel(x, y);
            }
        }
    }
    string outputName=output_path+filename;
    png.saveCanvasPng(const_cast<char*>(outputName.c_str()));
}
#endif
