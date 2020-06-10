#include "Gene.h"

Gene::Gene(const string& _ID, const string& genomeName, const int _coordinate,
           const bool _orientation) : ID(_ID), genomename(genomeName),
        coordinate(_coordinate), orientation(_orientation), is_tandem(false),
        is_tandem_representative(false), remapped(false), has_gf(false),
        pairs(NULL)
{

}

Gene::Gene() : pairs(NULL)
{

}

bool Gene::isPairWith(const Gene& gene) const
{
    if (has_gf) {
        if (gf_id == gene.getFamily() && ID != gene.getID())
            return true;
        else
            return false;
    } else {
        if (hasPairs() && pairs->find(string(gene.getID())) != pairs->end())
            return true;
    }

    return false;
}


bool Gene::isIndirectPairWith(const Gene& gene) const
{
    // if either of the genes has no pairs
    if ( (!hasPairs()) || (!gene.hasPairs()) )
        return false;

    // is it a direct pair?
    if (isPairWith(gene)) return true;

    // is it an indirect pair
    hash_set<string,stringhash>::const_iterator it = pairs->begin();
    for ( ; it != pairs->end(); it++)
        if (gene.pairs->find(*it) != gene.pairs->end())
            return true;

    return false;
}

void Gene::remapTo(Gene& gene)
{
    is_tandem = true;
    is_tandem_representative = false;
    tandem_representative = &gene;
    remapped = true;

    gene.is_tandem = true;
    gene.is_tandem_representative = true;
    gene.tandem_representative = &gene;
    gene.remapped = false;
}
