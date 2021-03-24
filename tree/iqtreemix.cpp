//
//  iqtreemix.cpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#include "iqtreemix.h"

IQTreeMix::IQTreeMix() : IQTree() {
    patn_freqs = NULL;
    patn_isconst = NULL;
    ptn_like_cat = NULL;
    _ptn_like_cat = NULL;
}

IQTreeMix::IQTreeMix(Alignment *aln) : IQTree(aln) {
    patn_freqs = NULL;
    patn_isconst = NULL;
    ptn_like_cat = NULL;
    _ptn_like_cat = NULL;
}

IQTreeMix::~IQTreeMix() {
    size_t i;
    for (i=0; i<size(); i++) {
        delete (at(i));
    }
    if (ptn_like_cat != NULL) {
        delete[] ptn_like_cat;
    }
    if (_ptn_like_cat != NULL) {
        delete[] _ptn_like_cat;
    }
    if (patn_freqs != NULL) {
        delete[] patn_freqs;
    }
    if (patn_isconst != NULL) {
        delete[] patn_isconst;
    }
}

// return how many char c inside the infile
int checkCharInFile(char* infile, char c) {
    ifstream fin;
    string aline;
    size_t i,k;
    k=0;
    fin.open(infile);
    while (getline(fin,aline)) {
        for (i=0; i<aline.length(); i++) {
            if (aline[i] == c)
                k++;
        }
    }
    fin.close();
    return k;
}

/**
        initialization
 */
void IQTreeMix::init(Alignment* aln, Params* params, Checkpoint* chkpt) {
    // read the trees inside the user input file
    cout << "call IQTreeMix::init" << endl;
    ifstream in;
    int numTree, i;

    // check how many trees inside the user input file
    numTree = checkCharInFile(params->user_file, ';');
    cout << "Number of input trees: " << numTree << endl;
    
    ASSERT(numTree);

    clear();
    in.open(params->user_file);
    for (i=0; i<numTree; i++) {
        IQTree* iqtree_new = new IQTree(aln);
        iqtree_new->setParams(params);
        iqtree_new->setCheckpoint(chkpt);
        iqtree_new->readUserTree(params->SSE, in);
        if (!params->compute_ml_tree_only) {
            iqtree_new->setRootNode(params->root);
        }
        iqtree_new->initSettings(*params);
        push_back(iqtree_new);
    }
    in.close();
    
    // initialize the tree weights
    weights.clear();
    double w = 1.0 / numTree;
    for (i=0; i<numTree; i++) {
        weights.push_back(w);
    }
    
    ptn_like_cat = new double[size() * aln->getNPattern()];
    _ptn_like_cat = new double[size() * aln->getNPattern()];
    patn_freqs = new int[aln->getNPattern()];
    patn_isconst = new int[aln->getNPattern()];
    
    // get the pattern frequencies
    aln->getPatternFreq(patn_freqs);
    
    // get whether the pattern is constant
    for (i=0; i<aln->getNPattern(); i++) {
        patn_isconst[i] = aln->at(i).isConst();
    }
}

void separateModel(string model_name, vector<string>& submodel_names) {
    cout << "[separateModel] model_name: " << model_name << endl << flush;
    size_t brac_start, brac_end, comma_pos, str_pos, i;
    string s;
    
    submodel_names.clear();
    
    // remove the last "+T"
    if (model_name.length() >= 2 && model_name.substr(model_name.length()-2,2)=="+T") {
        model_name = model_name.substr(0, model_name.length()-2);
    }
    brac_start = model_name.find("MIX{");
    if (brac_start != string::npos) {
        brac_end = model_name.find("}", brac_start+4);
    }
    if (brac_start == string::npos || brac_end == string::npos) {
        submodel_names.push_back(model_name);
        return;
    }
    str_pos = brac_start+4;
    comma_pos = model_name.find(",", str_pos);
    while (comma_pos != string::npos && comma_pos < brac_end) {
        submodel_names.push_back(model_name.substr(str_pos, comma_pos - str_pos));
        str_pos = comma_pos+1;
        comma_pos = model_name.find(",", str_pos);
    }
    if (str_pos < brac_end) {
        submodel_names.push_back(model_name.substr(str_pos, brac_end - str_pos));
    }
    if (brac_end+1 < model_name.length()) {
        s = model_name.substr(brac_end+1);
        for (i=0; i<submodel_names.size(); i++) {
            submodel_names[i].append(s);
        }
    }
    // show the model names
    cout << "[separateModel] submodel names:" << endl;
    for (i=0; i<submodel_names.size(); i++) {
        cout << submodel_names[i] << endl;
    }
}

void IQTreeMix::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    size_t i;
    vector<string> submodel_names;

    cout << "call IQTreeMix::initializeModel" << endl;
    cout << "[IQTreeMix::initializeModel] model_name: " << model_name << endl << flush;
    IQTree::initializeModel(params, model_name, models_block);
    ASSERT(getModelFactory()->model->isTreeMixture());
    separateModel(model_name, submodel_names);
    ASSERT(submodel_names.size()==size());
    for (i=0; i<size(); i++) {
        cout << "submodel " << i+1 << " : " << submodel_names[i] << endl;
        at(i)->initializeModel(params, submodel_names[i], models_block);
    }
}

double IQTreeMix::computeLikelihood(double *pattern_lh) {
    cout << "call IQTreeMix::computeLikelihood" << endl;
    double* pattern_lh_tree;
    size_t i,j,ptn,t;
    size_t nptn,ntree;
    double logLike = 0.0;
    double subLike;
    double score;
    
    nptn = aln->getNPattern();
    ntree = size();

    
    // show the tree weights
    cout << "[IQTreeMix::computeLikelihood] weights: ";
    for (t=0; t<ntree; t++) {
        if (t>0)
            cout << ", ";
        cout << weights[t];
    }
    cout << endl;
    
    // compute likelihood for each tree
    pattern_lh_tree = _ptn_like_cat;
    for (t=0; t<ntree; t++) {
        score = at(t)->computeLikelihood(pattern_lh_tree);
        // cout << "[IQTreeMix::computeLikelihood] Tree " << t+1 << " : " << score << endl;
        pattern_lh_tree += nptn;
    }
    
    // reorganize the array
    i=0;
    for (t=0; t<ntree; t++) {
        j=t;
        for (ptn=0; ptn<nptn; ptn++) {
            ptn_like_cat[j] = exp(_ptn_like_cat[i]);
            i++;
            j+=ntree;
        }
    }
    
    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        subLike = 0.0;
        for (t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        if (pattern_lh) {
            pattern_lh[ptn] = subLike;
        }
        // cout << ptn << "\t" << log(subLike) << "\t" << patn_freqs[ptn] << endl;
        logLike += log(subLike) * (double) patn_freqs[ptn];
    }
    cout << "[IQTreeMix::computeLikelihood] log-likelihood: " << logLike << endl;

    return logLike;
}

/**
        compute pattern likelihoods only if the accumulated scaling factor is non-zero.
        Otherwise, copy the pattern_lh attribute
        @param pattern_lh (OUT) pattern log-likelihoods,
                        assuming pattern_lh has the size of the number of patterns
        @param cur_logl current log-likelihood (for sanity check)
        @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
 */
void IQTreeMix::computePatternLikelihood(double *pattern_lh, double *cur_logl,
                                         double *pattern_lh_cat, SiteLoglType wsl) {
    size_t i,ptn,t;
    size_t nptn,ntree;
    double subLike;

    computeLikelihood();
    nptn = aln->getNPattern();
    ntree = size();

    // compute the total likelihood
    i=0;
    for (ptn=0; ptn<nptn; ptn++) {
        subLike = 0.0;
        for (t=0; t<ntree; t++) {
            subLike += ptn_like_cat[i] * weights[t];
            i++;
        }
        pattern_lh[ptn] = subLike;
    }
}

int IQTreeMix::ensureNumberOfThreadsIsSet(Params *params) {
    size_t i;
    // int k = IQTree::ensureNumberOfThreadsIsSet(params);
    int k = 1;
    for (i=0; i<size(); i++) {
        k &= at(i)->ensureNumberOfThreadsIsSet(params);
    }
    return k;
}

void IQTreeMix::initializeAllPartialLh() {
    size_t i;
    // IQTree::initializeAllPartialLh();
    for (i=0; i<size(); i++) {
        at(i)->initializeAllPartialLh();
    }
}

void IQTreeMix::deleteAllPartialLh() {
    size_t i;
    // IQTree::deleteAllPartialLh();
    for (i=0; i<size(); i++) {
        at(i)->deleteAllPartialLh();
    }
}

void IQTreeMix::clearAllPartialLH(bool make_null) {
    size_t i;
    // IQTree::clearAllPartialLH(make_null);
    for (i=0; i<size(); i++) {
        at(i)->clearAllPartialLH(make_null);
    }
}

/**
        optimize all branch lengths of the tree
        @param iterations number of iterations to loop through all branches
        @return the likelihood of the tree
 */
double IQTreeMix::optimizeAllBranches(int my_iterations, double tolerance, int maxNRStep) {
    size_t i;
    for (i=0; i<size(); i++) {
        at(i)->optimizeAllBranches(my_iterations, tolerance, maxNRStep);
    }
    return computeLikelihood();
}

void IQTreeMix::showTree() {
    size_t i;
    for (i=0; i<size(); i++) {
        cout << "Tree " << i+1 << ": ";
        at(i)->printTree(cout);
        cout << endl;
    }
}
