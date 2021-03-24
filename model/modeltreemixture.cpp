/*
 * modeltreemixture.cpp
 *
 *  Created on: Dec 2, 2020
 *      Author: thomas
 */

#include "modeltreemixture.h"

using namespace std;

const double MIN_MIXTURE_PROP = 0.001;

/**
    constructor
    @param tree associated tree for the model
*/
ModelTreeMixture::ModelTreeMixture(IQTreeMix* treemix) : ModelMixture(treemix), ModelMarkov(treemix) {
    prop = NULL;
    fix_prop = true;
    optimizing_submodels = false;
    mixtrees = treemix;
}

ModelTreeMixture::ModelTreeMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
        StateFreqType freq, string freq_params, IQTreeMix* treemix, bool optimize_weights)
        : ModelMixture(treemix), ModelMarkov(treemix)
{
    int nptn;
    int i;
    PhyloTree *tree;
    int nmix = size();
    prop = NULL;
    fix_prop = true;
    optimizing_submodels = false;
    optimize_steps = 0;
    mixtrees = treemix;
    initMixture(orig_model_name, model_name, model_list, models_block, freq, freq_params, treemix, optimize_weights);
}

ModelTreeMixture::~ModelTreeMixture() {
}

double ModelTreeMixture::optimizeParameters(double gradient_epsilon) {

    double score = 0.0;

    if (!phylo_tree->getModelFactory()->unobserved_ptns.empty()) {
        outError("Mixture model +ASC is not supported yet. Contact author if needed.");
    }
    
    score = optimizeWithEM(gradient_epsilon);
    cout << "[ModelTreeMixture::optimizeParameters] finished optimizeWithEM()" << endl;

    // now rescale Q matrices to have proper interpretation of branch lengths

    double sum;
    int i, ncategory = size();
    for (i = 0, sum = 0.0; i < ncategory; i++) {
        sum += mixtrees->weights[i] * at(i)->total_num_subst;
    }

    if (fabs(sum-1.0) > 1e-6 && !isFused()) {
        for (i = 0; i < ncategory; i++)
            at(i)->total_num_subst /= sum;
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
    }
    return score;
}


double ModelTreeMixture::optimizeWithEM(double gradient_epsilon) {
    cout << "ModelTreeMixture::optimizeWithEM()" << endl << flush;
    
    size_t nptn, ptn, c, i;
    size_t nmix = size();
    // vector<PhyloTree*> trees;
    // vector<ModelFactory*> models;
    // vector<RateHeterogeneity*> site_rates;
    PhyloTree *tree;
    ModelFactory *model_fac;
    IQTreeMix* iqtreeMix;
    double* pattern_mix_lh;
    bool converged1, converged2;
    double logl_epsilon = 0.001;

    ASSERT(nmix);
    
    nptn = phylo_tree->aln->getNPattern();

    double *new_prop = aligned_alloc<double>(nmix);
    pattern_mix_lh = new double[nmix * nptn];
    
    double prev_score = -DBL_MAX, score;

//    int num_steps = 100000; //SC
    int step;

    converged1=converged2=false;

    // EM algorithm loop described in Wang, Li, Susko, and Roger (2008)
    for (step = 0; step < optimize_steps; step++) {
        cout << "step " << step << endl;
        
        phylo_tree->initializeAllPartialLh();
        score = phylo_tree->computeLikelihood();
        phylo_tree->clearAllPartialLH();

        converged1 = (score < prev_score + gradient_epsilon);
        
        if (converged1 && converged2) {
            break;
        }
        
        prev_score = score;
        
        memcpy(pattern_mix_lh, mixtrees->ptn_like_cat, nptn * nmix * sizeof(double));
        // multiply pattern_mix_lh by tree weights
        i = 0;
        for (ptn = 0; ptn < nptn; ptn++) {
            for (c = 0; c < nmix; c++) {
                pattern_mix_lh[i] *= mixtrees->weights[c];
                i++;
            }
        }

        memset(new_prop, 0, nmix*sizeof(double));

        // E-step
        // cout << "E-step" << endl << flush;
        // decoupled weights (prop) from pattern_mix_lh to obtain L_ci and compute pattern likelihood L_i
        for (ptn = 0; ptn < nptn; ptn++) {
            double *this_lk_cat = pattern_mix_lh + ptn*nmix;
            // double lk_ptn = phylo_tree->ptn_invar[ptn];
            double lk_ptn = 0.0;
            for (c = 0; c < nmix; c++) {
                lk_ptn += this_lk_cat[c];
            }
            ASSERT(lk_ptn != 0.0);
            lk_ptn = mixtrees->patn_freqs[ptn] / lk_ptn;

            // transform pattern_mix_lh into posterior probabilities of each category
            for (c = 0; c < nmix; c++) {
                this_lk_cat[c] *= lk_ptn;
                new_prop[c] += this_lk_cat[c];
            }
        }
        
        // M-step
        for (c = 0; c < nmix; c++) {
            new_prop[c] = new_prop[c] / phylo_tree->getAlnNSite();
            if (new_prop[c] < 1e-10) new_prop[c] = 1e-10;
        }
        
        // update the weights for rate model
        // and check for convergence
        converged2 = true;
        if (isFused() && !phylo_tree->getRate()->getFixParams()) {
            double new_pinvar = 0.0;
            for (c = 0; c < nmix; c++) {
                converged2 = converged2 && (fabs(mixtrees->weights[c] - new_prop[c]) < 1e-4);
                mixtrees->weights[c] = new_prop[c];
                prop[c] = new_prop[c];
                new_pinvar += new_prop[c];
            }
            new_pinvar = fabs(1.0 - new_pinvar);
            if (new_pinvar > 1e-6) {
                converged2 = converged2 && (fabs(phylo_tree->getRate()->getPInvar() - new_pinvar) < 1e-4);
                phylo_tree->getRate()->setPInvar(new_pinvar);
                phylo_tree->computePtnInvar();
            }
        } else if (!fix_prop) {
            for (c = 0; c < nmix; c++) {
                // check for convergence
                converged2 = converged2 && (fabs(mixtrees->weights[c]-new_prop[c]) < 1e-4);
                mixtrees->weights[c] = new_prop[c];
                prop[c] = new_prop[c];
            }
        }

        // optimize model one by one
        for (c = 0; c < nmix; c++) if (at(c)->getNDim() > 0) {
            tree = at(c)->phylo_tree; // trees[c];
            model_fac = at(c)->phylo_tree->model_factory; // models[c];
            
            // initialize likelihood
            tree->initializeAllPartialLh();
            // copy posterior probability into ptn_freq
            tree->computePtnFreq();
            
            double *this_lk_cat = pattern_mix_lh+c;
            for (ptn = 0; ptn < nptn; ptn++)
                tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
            cout << "[ModelTreeMixture::optimizeWithEM] start optimizing the model " << c+1 << ":" << endl;
            model_fac->optimizeParameters(gradient_epsilon);
            // at(c)->optimizeParameters(gradient_epsilon);
            cout << "[ModelTreeMixture::optimizeWithEM] done optimizing the model " << c+1 << ":" << endl;
        }

        
        phylo_tree->initializeAllPartialLh();
        score = phylo_tree->computeLikelihood();
        
        // update the ptn_freq array
        for (c = 0; c < nmix; c++) {
            tree = at(c)->phylo_tree;
            // initialize likelihood
            tree->initializeAllPartialLh();
            // copy posterior probability into ptn_freq
            tree->computePtnFreq();
            double *this_lk_cat = pattern_mix_lh+c;
            for (ptn = 0; ptn < nptn; ptn++)
                tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
        }
    }

    phylo_tree->initializeAllPartialLh();
    score = phylo_tree->computeLikelihood();
    
    // update the ptn_freq array
    for (c = 0; c < nmix; c++) {
        tree = at(c)->phylo_tree;
        // initialize likelihood
        tree->initializeAllPartialLh();
        // copy posterior probability into ptn_freq
        tree->computePtnFreq();
        double *this_lk_cat = pattern_mix_lh+c;
        for (ptn = 0; ptn < nptn; ptn++)
            tree->ptn_freq[ptn] = this_lk_cat[ptn*nmix];
    }
    
    // show the information
    cout << "[ModelTreeMixture::optimizeWithEM] log-likelihood: " << score << endl;
    cout << "[ModelTreeMixture::optimizeWithEM] weights: ";
    for (c=0; c<nmix; c++) {
        if (c>0)
            cout << ", ";
        cout << mixtrees->weights[c];
    }
    cout << endl;
    for (c=0; c<nmix; c++) {
        cout << "[ModelTreeMixture::optimizeWithEM] Tree " << c+1 << ": " << at(c)->phylo_tree->getTreeString() << endl;
        cout << "[ModelTreeMixture::optimizeWithEM] Model " << c+1 << ": " << endl;
        at(c)->phylo_tree->model_factory->model->writeInfo(cout);
        cout << endl;
    }

    aligned_free(new_prop);
    delete[] pattern_mix_lh;

    return score;
}

/**
    write information
    @param out output stream
*/
void ModelTreeMixture::writeInfo(ostream &out) {
    int c;
    int nmix = size();
    for (c=0; c<nmix; c++) {
        out << "Tree-mixture model " << c+1 << ":" << endl;
        at(c)->phylo_tree->model_factory->model->writeInfo(out);
        if (at(c)->phylo_tree->getRate()) {
            at(c)->phylo_tree->getRate()->writeInfo(out);
        }
    }
}
