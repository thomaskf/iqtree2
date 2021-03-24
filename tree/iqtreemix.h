//
//  iqtreemix.hpp
//  tree
//
//  Created by Thomas Wong on 14/12/20.
//

#ifndef iqtreemix_h
#define iqtreemix_h

#include <stdio.h>
#include <vector>
#include "iqtree.h"

class IQTreeMix : public IQTree, public vector<IQTree*> {
public:

    /**
            default constructor
     */
    IQTreeMix();

    IQTreeMix(Alignment *aln);

    /**
            destructor
     */
    virtual ~IQTreeMix();

    /**
            initialization
     */
    void init(Alignment* aln, Params* params, Checkpoint* chkpt);
    
    void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    virtual double computeLikelihood(double *pattern_lh = NULL);

    /**
            compute pattern likelihoods only if the accumulated scaling factor is non-zero.
            Otherwise, copy the pattern_lh attribute
            @param pattern_lh (OUT) pattern log-likelihoods,
                            assuming pattern_lh has the size of the number of patterns
            @param cur_logl current log-likelihood (for sanity check)
            @param pattern_lh_cat (OUT) if not NULL, store all pattern-likelihood per category
     */
    virtual void computePatternLikelihood(double *pattern_lh, double *cur_logl = NULL,
            double *pattern_lh_cat = NULL, SiteLoglType wsl = WSL_RATECAT);

    virtual int ensureNumberOfThreadsIsSet(Params *params);
    
    virtual void initializeAllPartialLh();

    virtual void deleteAllPartialLh();

    virtual void clearAllPartialLH(bool make_null = false);

    /**
            optimize all branch lengths of the tree
            @param iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);

    virtual void showTree();
    
    /**
            pattern frequencies
     */
    int* patn_freqs;
    
    /**
            whether pattern is constant
     */
    int* patn_isconst;
    
    /**
            weights of trees
     */
    vector<double> weights;
    
    /**
            pattern likelihoods for all trees
     */
    double* ptn_like_cat;
    
private:
    
    /**
            immediate array for pattern likelihoods during computation
     */
    double* _ptn_like_cat;
};

#endif /* iqtreemix_h */
