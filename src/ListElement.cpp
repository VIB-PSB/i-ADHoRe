#include "ListElement.h"

ListElement::ListElement (Gene& _gene, bool _orientation, bool _masked)
        : gene (_gene), orientation(_orientation), masked(_masked), gap(false),
        numID(-1), hasHomolog(false), hasAP(false), hasAlHomolog(false),
        hasAlAP(false) { }

ListElement::ListElement() : orientation(false), masked(false), gap(true),
        numID(-1), hasHomolog(false), hasAP(false), hasAlHomolog(false),
        hasAlAP(false) { }

const Gene& ListElement::getGene() const {
    return gene;
}

Gene& ListElement::getGene() {
    return gene;
}

void ListElement::matchingPositions(const vector<ListElement*>& list,
                                    queue<int>& q) const
{
    for (unsigned int i = 0; i < list.size(); i++) {
        if (list[i]->isGap()) continue;
        if (list[i]->isMasked()) continue;

        if (gene.isPairWith(list[i]->getGene()))
            q.push(i);
    }
}

bool ListElement::isMasked() const {
    return masked;
}

bool ListElement::getOrientation() const {
    return orientation;
}
