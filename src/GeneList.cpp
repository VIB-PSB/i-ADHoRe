#include "GeneList.h"

#include "Gene.h"
#include "ListElement.h"

#include "debug/FileException.h"
#include <cassert>

GeneList::GeneList(const string& listName, const string& genomeName,
                   const string& fileName) :
        is_segment(false), listname(listName), genomename(genomeName), id(-1)
{
    ifstream fin (fileName.c_str());

    if (!fin) throw FileException("Could not open gene list file: " +
                                      fileName + "!");

    int coordinate = 0;
    string ID;

    while (fin >> ID) {
        char c = ID[ID.length()-1];  // final character in ID

        if (c != '+' && c != '-') {
            cerr << "No orientation +/- in gene list." << endl;
            exit(EXIT_FAILURE);
        }

        bool orientation = (c == '+');
        ID.erase(ID.length()-1);    // erase final character

        Gene gene (ID, genomename, coordinate, orientation);
        elements.push_back(new ListElement(gene, orientation, false));
        coordinate++;
    }

    fin.close();
}

GeneList::GeneList(const string& listName, const string& genomeName, const vector< ListElement* >& segmentFromFile)
    : is_segment(true), listname(listName), genomename(genomeName), id(-1)
{
    for (int i=0; i<segmentFromFile.size(); i++){
        remapped_elements.push_back(new ListElement(*segmentFromFile[i]));
    }
}

GeneList::GeneList(int size) : is_segment(true), id(-1)
{
    Gene gene;
    for (int i = 0; i < size; i++)
        remapped_elements.push_back(new ListElement(gene, 0, false));
}

GeneList::GeneList(const GeneList& genelist, int begin, int end) :
    is_segment(true), listname(genelist.listname),
    genomename(genelist.genomename), id(genelist.getID())
{
    if (end < begin) return;

    for (int i = begin; i <= end; i++)
        remapped_elements.push_back(new ListElement (*genelist.getRemappedElements()[i]));
}

GeneList::~GeneList ()
{
    if (isSegment()) {
        for (int i = 0; i < remapped_elements.size(); i++)
            delete remapped_elements[i];
    } else {
        for (int i = 0; i < elements.size(); i++)
            delete elements[i];
    }
}

void GeneList::remapTandems(int gapSize, bool useFamily)
{
    if (useFamily)
        remapTandemsFamily(gapSize);
    else
        remapTandemsPairs(gapSize);
}

void GeneList::remapTandemsPairs(int gapSize)
{
    // map positive genes upon positive genes (back to front)
    for (int i = elements.size() - 1; i > 0; i--) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getOrientation() == false) continue;

        for (int j = i - 1; j >= 0 && i - j <= gapSize + 1; j--) {
            if (elements[j]->getOrientation() == false) continue;
            if (elements[i]->getGene().isIndirectPairWith(elements[j]->getGene())) {
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    // map negative genes upon negative genes (front to back)
    for (int i = 0; i < elements.size() - 1; i++) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getOrientation() == true) continue;

        for (int j = i + 1; j < elements.size() && j - i <= gapSize + 1; j++) {
            if (elements[j]->getOrientation() == true) continue;
            if (elements[i]->getGene().isIndirectPairWith(elements[j]->getGene())) {
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    // map opposite genes (back to front)
    for (int i = elements.size() - 1; i > 0; i--) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getGene().isRemapped()) continue;

        for (int j = i - 1; j >= 0 && i - j <= gapSize + 1; j--) {
            if (elements[i]->getGene().isIndirectPairWith(elements[j]->getGene())) {
                // gene i and gene j will have a different orientation
                // otherwise, they would have been remapped above
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    //calculate remapped coordinate for each gene
    int substract = 0;
    for (unsigned int i = 0; i < elements.size(); i++) {
        if (elements[i]->getGene().isRemapped()) {
            substract++;
        }
        else {
            remapped_elements.push_back(elements[i]);
        }

        elements[i]->getGene().setRemappedCoordinate(substract);
    }
}

void GeneList::remapTandemsFamily(int gapSize)
{
    // map positive genes upon positive genes (back to front)
    for (int i = elements.size() - 1; i > 0; i--) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getOrientation() == false) continue;

        for (int j = i - 1; j >= 0 && i - j <= gapSize + 1; j--) {
            if (elements[j]->getOrientation() == false) continue;
            if (elements[i]->getGene().isPairWith(elements[j]->getGene())) {
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    // map negative genes upon negative genes (front to back)
    for (int i = 0; i < elements.size() - 1; i++) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getOrientation() == true) continue;

        for (int j = i + 1; j < elements.size() && j - i <= gapSize + 1; j++) {
            if (elements[j]->getOrientation() == true) continue;
            if (elements[i]->getGene().isPairWith(elements[j]->getGene())) {
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    // map opposite genes (back to front)
    for (int i = elements.size() - 1; i > 0; i--) {
        if (!elements[i]->getGene().hasPairs()) continue;
        if (elements[i]->getGene().isRemapped()) continue;

        for (int j = i - 1; j >= 0 && i - j <= gapSize + 1; j--) {
            if (elements[i]->getGene().isPairWith(elements[j]->getGene())) {
                // gene i and gene j will have a different orientation
                // otherwise, they would have been remapped above
                elements[i]->getGene().remapTo(elements[j]->getGene());
                break;
            }
        }
    }

    //calculate remapped coordinate for each gene
    int substract = 0;
    for (unsigned int i = 0; i < elements.size(); i++) {
        if (elements[i]->getGene().isRemapped()) {
            substract++;
        }
        else {
            remapped_elements.push_back(elements[i]);
        }

        elements[i]->getGene().setRemappedCoordinate(substract);
    }
}

int GeneList::getNumberOfMaskedElements() const {
    int count = 0;

    for (unsigned int i = 0; i < remapped_elements.size(); i++) {
        if (remapped_elements[i]->isMasked()) {
            count++;
        }
    }
    return count;
}

void GeneList::introduceGaps(int begin, int end, int number_of_gaps)
{
    assert(begin <= end);
    assert(end <= remapped_elements.size());
    if (number_of_gaps < 1) return;

    vector<ListElement*>::iterator it;

    // special case
    if (begin == end) {
        if (end == remapped_elements.size())
            it = remapped_elements.end();
        else
            it = std::find (remapped_elements.begin(),
                            remapped_elements.end(),
                            remapped_elements[end]);
        for (int i = 0; i < number_of_gaps; i++, it++) {
            ListElement *gap = new ListElement();
            it = remapped_elements.insert(it, gap);
        }
        return;
    }

    // divide the gaps evenly among the available positions
    int number_of_elements = end - begin;
    double avg_gaps =  (double)number_of_gaps / (double)number_of_elements;

    it = std::find (remapped_elements.begin(), remapped_elements.end(),
                    remapped_elements[begin+1]);

    double gapCount = 0.5;
    for (int i = 0; i < number_of_elements; i++, it++) {
        gapCount += avg_gaps;
        while (gapCount > 1.0) {
            gapCount -= 1.0;
            ListElement *gap = new ListElement();
            it = remapped_elements.insert(it, gap);
            it++;
        }
    }
}

void GeneList::introduceGap(int position)
{
    vector<ListElement*>::iterator it = find (remapped_elements.begin(),
                                              remapped_elements.end(),
                                              remapped_elements[position]);

    //insert gap before position
    remapped_elements.insert(it, new ListElement());
}

void GeneList::pasteGap()
{
    remapped_elements.push_back(new ListElement());
}

void GeneList::invertSection(int begin, int end)
{
    vector<ListElement*>::iterator it_begin =
        remapped_elements.begin() + begin;
    vector<ListElement*>::iterator it_end =
        remapped_elements.begin() + end + 1;

    for (vector<ListElement*>::iterator it = it_begin; it != it_end; it++)
        (*it)->invertOrientation();

    reverse(it_begin, it_end);
}

void GeneList::mask(int begin, int end) {
    for (int i = begin; i <= end; i++) {
        remapped_elements[i]->setMasked(true);
    }
}

void GeneList::write() {
    cerr << "--> ";
    for (unsigned int i = 0; i < remapped_elements.size(); i++) {
        if (remapped_elements[i]->isGap()) {
            cerr << setw(3) << "0";
        }
        else {
            cerr << setw(3) << "1";
        }
    }
}
void GeneList::write2() {
    cout << "-->\t";
    for (unsigned int i = 0; i < remapped_elements.size(); i++) {
        if (remapped_elements[i]->isGap()) {
            cout <<"gap\t";
        }
        else {
            cout << remapped_elements[i]->getGene().getID() << "\t";
        }
    }
    cout << endl;
}

const string GeneList::getGeneName(int i) {
    string out = "not found";
    if (i >= 0 && i < (int)remapped_elements.size())
    {
        if (remapped_elements[i]->isGap()) {
            out = "gap";
            return out;
        }
        else {
            return remapped_elements[i]->getGene().getID();
        }
    }

    return out;
}

void GeneList::removeGaps()
{
    vector<ListElement*>::iterator it = remapped_elements.begin();
    while (it != remapped_elements.end()) {
        if ((*it)->isGap()) {
            delete (*it);
            it = remapped_elements.erase(it);
        }
        else {
            it++;
        }
    }
}


