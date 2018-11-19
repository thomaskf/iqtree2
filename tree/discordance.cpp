//
//  discordance.cpp
//  tree
//
//  Created by Minh Bui on 24/9/18.
//

#include "phylosupertree.h"

void PhyloTree::computeSiteConcordance(map<string,string> &meanings) {
    BranchVector branches;
    getInnerBranches(branches);
    for (auto it = branches.begin(); it != branches.end(); it++) {
        computeSiteConcordance((*it), params->site_concordance, randstream);
        Neighbor *nei = it->second->findNeighbor(it->first);
        double sCF;
        GET_ATTR(nei, sCF);

        string sup_str = convertDoubleToString(sCF);
        Node *node = it->second;
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += "[&sCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",!sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }
    meanings.insert({"sCF", "Site concordance factor (%) averaged over " + convertIntToString(params->site_concordance) +  " quartets"});
    meanings.insert({"sDF1", "Site discordance factor (%) for alternative quartet 1"});
    meanings.insert({"sDF2", "Site discordance factor (%) for alternative quartet 2"});
    meanings.insert({"sN", "Number of informative sites averaged over " + convertIntToString(params->site_concordance) +  " quartets"});
}

void Alignment::computeQuartetSupports(IntVector &quartet, vector<size_t> &support) {
    for (auto pat = begin(); pat != end(); pat++) {
        if (!pat->isInformative()) continue;
        bool informative = true;
        for (int j = 0; j < quartet.size(); j++)
            if (pat->at(quartet[j]) >= num_states) {
                informative = false;
                break;
            }
        if (!informative) continue;
        if (pat->at(quartet[0]) == pat->at(quartet[1]) && pat->at(quartet[2]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[2]))
            support[0] += pat->frequency;
        if (pat->at(quartet[0]) == pat->at(quartet[2]) && pat->at(quartet[1]) == pat->at(quartet[3]) && pat->at(quartet[0]) != pat->at(quartet[1]))
            support[1] += pat->frequency;
        if (pat->at(quartet[0]) == pat->at(quartet[3]) && pat->at(quartet[1]) == pat->at(quartet[2]) && pat->at(quartet[0]) != pat->at(quartet[1]))
            support[2] += pat->frequency;
    }
}

void SuperAlignment::computeQuartetSupports(IntVector &quartet, vector<size_t> &support) {
    for (int part = 0; part < partitions.size(); part++) {
        IntVector part_quartet;
        for (auto i = quartet.begin(); i != quartet.end(); i++) {
            if (taxa_index[*i][part] >= 0)
                part_quartet.push_back(taxa_index[*i][part]);
            else
                break;
        }
        if (part_quartet.size() != quartet.size())
            continue;
        if (Params::getInstance().site_concordance_partition) {
            vector<size_t> part_support;
            part_support.resize(3, 0);
            partitions[part]->computeQuartetSupports(part_quartet, part_support);
            for (int j = 0; j < 3; j++) {
                support[part*3+3+j] += part_support[j];
                support[j] += part_support[j];
            }
        } else
            partitions[part]->computeQuartetSupports(part_quartet, support);
    }
}

void PhyloTree::computeSiteConcordance(Branch &branch, int nquartets, int *rstream) {
    vector<IntVector> taxa;
    taxa.resize(4);

    // extract the taxa from the two left subtrees
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        getTaxaID(taxa[id], (*it)->node, branch.first);
        id++;
        if (id > 2)
            outError(__func__, " only work with bifurcating tree");
    }

    // extract the taxa from the two right subtrees
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        getTaxaID(taxa[id], (*it)->node, branch.second);
        id++;
        if (id > 4)
            outError(__func__, " only work with bifurcating tree");
    }
    
    double sCF = 0.0; // concordance factor
    double sDF1 = 0.0;
    double sDF2 = 0.0;
    double sN = 0.0;
    size_t sum_sites = 0;
    int i;
    vector<size_t> support;
    support.resize(3, 0);
    // reserve size for partition-wise concordant/discordant sites
    if (Params::getInstance().site_concordance_partition && isSuperTree())
        support.resize(((PhyloSuperTree*)this)->size()*3+3, 0);
    for (i = 0; i < nquartets; i++) {
        int j;
        // get a random quartet
        IntVector quartet;
        quartet.resize(taxa.size());
        for (j = 0; j < taxa.size(); j++) {
            quartet[j] = taxa[j][random_int(taxa[j].size(), rstream)];
        }
        support[0] = support[1] = support[2] = 0;
        aln->computeQuartetSupports(quartet, support);
        size_t sum = support[0] + support[1] + support[2];
        sum_sites += sum;
        if (sum > 0) {
            sCF += ((double)support[0]) / sum;
            sDF1 += ((double)support[1]) / sum;
            sDF2 += ((double)support[2]) / sum;
        }
    }
    sN = (double)sum_sites / nquartets;
    sCF = round(sCF / nquartets * 10000)/100;
    sDF1 = round(sDF1 / nquartets * 10000)/100;
    sDF2 = round(sDF2 / nquartets * 10000)/100;
    Neighbor *nei = branch.second->findNeighbor(branch.first);
    PUT_ATTR(nei, sCF);
    PUT_ATTR(nei, sN);
    PUT_ATTR(nei, sDF1);
    PUT_ATTR(nei, sDF2);
    // insert key-value for partition-wise con/discordant sites
    for (i = 3; i < support.size(); i++) {
        string key;
        switch (i%3) {
            case 0: key = "sC"; break;
            case 1: key = "sD1"; break;
            case 2: key = "sD2"; break;
        }
        key = key + convertIntToString(i/3);
        nei->putAttr(key, (double)support[i]/nquartets);
    }
}

/**
 assign branch supports to a target tree
 */
void PhyloTree::computeGeneConcordance(MTreeSet &trees, map<string,string> &meanings) {
    StrVector names;
    getTaxaName(names);
    StringIntMap name_map;
    for (auto stri = names.begin(); stri != names.end(); stri++)
        name_map[*stri] = stri - names.begin();
    BranchVector branches;
    vector<Split*> subtrees;
    extractQuadSubtrees(subtrees, branches, root->neighbors[0]->node);
    IntVector decisive_counts; // number of decisive trees
    decisive_counts.resize(branches.size(), 0);
    IntVector supports[3]; // number of trees supporting 3 alternative splits
    supports[0].resize(branches.size(), 0);
    supports[1].resize(branches.size(), 0);
    supports[2].resize(branches.size(), 0);
    StrVector occurence_trees; // list of tree IDs where each split occurs
    if (verbose_mode >= VB_MED)
        occurence_trees.resize(branches.size());
    int treeid, taxid;
    for (treeid = 0; treeid < trees.size(); treeid++) {
        MTree *tree = trees[treeid];
        StrVector taxname;
        tree->getTaxaName(taxname);
        // create the map from taxa between 2 trees
        Split taxa_mask(leafNum);
        for (StrVector::iterator it = taxname.begin(); it != taxname.end(); it++) {
            if (name_map.find(*it) == name_map.end())
                outError("Taxon not found in full tree: ", *it);
            taxa_mask.addTaxon(name_map[*it]);
        }
        // make the taxa ordering right before converting to split system
        taxname.clear();
        int smallid;
        for (taxid = 0, smallid = 0; taxid < leafNum; taxid++)
            if (taxa_mask.containTaxon(taxid)) {
                taxname.push_back(names[taxid]);
                tree->findLeafName(names[taxid])->id = smallid++;
            }
        ASSERT(taxname.size() == tree->leafNum);
        
        SplitGraph sg;
        //NodeVector nodes;
        tree->convertSplits(sg);
        SplitIntMap hash_ss;
        for (auto sit = sg.begin(); sit != sg.end(); sit++)
            hash_ss.insertSplit((*sit), 1);
        
        // now scan through all splits in current tree
        int id, qid;
        for (id = 0, qid = 0; qid < subtrees.size(); id++, qid += 4)
        {
            bool decisive = true;
            int i;
            for (i = 0; i < 4; i++) {
                if (!taxa_mask.overlap(*subtrees[qid+i])) {
                    decisive = false;
                    break;
                }
            }
            if (!decisive) continue;
            
            decisive_counts[id]++;
            for (i = 0; i < 3; i++) {
                Split this_split = *subtrees[qid]; // current split
                this_split += *subtrees[qid+i+1];
                Split *subsp = this_split.extractSubSplit(taxa_mask);
                if (subsp->shouldInvert())
                    subsp->invert();
                if (hash_ss.findSplit(subsp)) {
                    supports[i][id]++;
                    if (verbose_mode >= VB_MED && i == 0)
                        occurence_trees[id] += convertIntToString(treeid+1) + " ";
                }
                delete subsp;
            }
        }
        
    }
    
    for (int i = 0; i < branches.size(); i++) {
        Neighbor *nei = branches[i].second->findNeighbor(branches[i].first);
        double gCF = round((double)supports[0][i]/decisive_counts[i] * 10000)/100;
        double gDF1 = round((double)supports[1][i]/decisive_counts[i] * 10000)/100;
        double gDF2 = round((double)supports[2][i]/decisive_counts[i] * 10000)/100;
        int gN = decisive_counts[i];
        PUT_ATTR(nei, gCF);
        PUT_ATTR(nei, gDF1);
        PUT_ATTR(nei, gDF2);
        PUT_ATTR(nei, gN);
        
        stringstream tmp;
        tmp.precision(3);
        tmp << (double)supports[0][i]/decisive_counts[i]*100;
        if (verbose_mode >= VB_MED)
            tmp << "%" << decisive_counts[i];
        
        Node *node = branches[i].second;
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']')
                node->name += "[&CF=" + tmp.str() + "]";
            else
                node->name = node->name.substr(0, node->name.length()-1) + ",!CF=" + tmp.str() + "]";
        } else {
            if (!node->name.empty())
                node->name.append("/");
            node->name.append(tmp.str());
        }
        if (verbose_mode >= VB_MED) {
            cout << node->name << " " << occurence_trees[i] << endl;
        }
    }
    for (vector<Split*>::reverse_iterator it = subtrees.rbegin(); it != subtrees.rend(); it++)
        delete (*it);

    meanings.insert({"gCF", "Gene concordance factor (%)"});
    meanings.insert({"gDF1", "Gene discordance factor (%) for alternative split 1"});
    meanings.insert({"gDF2", "Gene discordance factor (%) for alternative split 2"});
    meanings.insert({"gN", "Number of trees decisive for the branch"});
}

/**
 compute quartet internode certainty, similar to Zhou et al (biorxiv)
 */
void PhyloTree::computeQuartetConcordance(MTreeSet &trees) {
    
    outError("Not working yet, need consent from Zhou et al.");
    BranchVector branches;
    getInnerBranches(branches);
    
    for (auto treeit = trees.begin(); treeit != trees.end(); treeit++) {
        
    }
    
    for (auto it = branches.begin(); it != branches.end(); it++) {
        Node *node = it->second;
        double sup = computeQuartetConcordance(*it, trees);
        string sup_str = convertDoubleToString(sup);
        
        if (Params::getInstance().newick_extended_format) {
            if (node->name.empty() || node->name.back() != ']') {
                node->name += "[&qCF=" + sup_str + "]";
            } else
                node->name = node->name.substr(0, node->name.length()-1) + ",!sCF=" + sup_str + "]";
        } else {
            if (!node->name.empty())
                node->name += "/";
            node->name += sup_str;
        }
    }

}

double PhyloTree::computeQuartetConcordance(Branch &branch, MTreeSet &trees) {
    vector<IntVector> taxa;
    taxa.resize(4);
    if (branch.first->degree() != 3 || branch.second->degree() != 3)
        outError(__func__, " only work with bifurcating tree");

    // extract the taxa from the two left subtrees
    int id = 0;
    FOR_NEIGHBOR_DECLARE(branch.first, branch.second, it) {
        getTaxaID(taxa[id], (*it)->node, branch.first);
        id++;
    }
    
    // extract the taxa from the two right subtrees
    FOR_NEIGHBOR(branch.second, branch.first, it) {
        getTaxaID(taxa[id], (*it)->node, branch.second);
        id++;
    }
    
    double sum_support = 0.0;
    int num_quartets = Params::getInstance().site_concordance; // TODO: change name
    for (int i = 0; i < num_quartets; i++) {
        int j;
        // get a random quartet
        IntVector quartet;
        quartet.resize(taxa.size());
        for (j = 0; j < taxa.size(); j++) {
            quartet[j] = taxa[j][random_int(taxa[j].size())];
        }
        
        int quartetCF[3] = {0, 0, 0};
        for (auto tree = trees.begin(); tree != trees.end(); tree++) {
        }
        int sum = quartetCF[0] + quartetCF[1] + quartetCF[2];
        if (sum > 0)
            sum_support += ((double)quartetCF[0]) / sum;
    }
    return sum_support / num_quartets;
}
