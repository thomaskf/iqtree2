/*
 * modeltreemixture.h
 *
 *  Created on: Dec 2, 2020
 *      Author: thomas
 */

#ifndef MODELTREEMIXTURE_H_
#define MODELTREEMIXTURE_H_

#include "modelmixture.h"
#include "tree/iqtreemix.h"

/**
 * tree mixture model
 */
class ModelTreeMixture : public ModelMixture {
public:
    
    /**
        constructor
        @param model_name model name, e.g., JC, HKY.
        @param freq state frequency type
        @param trees associated phylogenetic trees
    */
    ModelTreeMixture(string orig_model_name, string model_name, string model_list, ModelsBlock *models_block,
            StateFreqType freq, string freq_params, IQTreeMix *treemix, bool optimize_weights);

    /**
        constructor
        @param trees associated trees for the model
    */
    ModelTreeMixture(IQTreeMix *treemix);

    virtual ~ModelTreeMixture();

    /**
     * @return the corresponding tree for the i-th model
     * For this model mixture, it should return the same tree.
     * For tree-mixture model, which inherits this class,
     * it will return different tree for different model
     */
    virtual PhyloTree* getTree(int i) {
        if (mixtrees == nullptr || i >= mixtrees->size()) {
            outError("The number of trees inputed is less than the number of models specified");
        }
        return mixtrees->at(i);
    }

    /**
        optimize rate parameters using EM algorithm
        @param gradient_epsilon
        @return log-likelihood of optimized parameters
    */
    double optimizeWithEM(double gradient_epsilon);

    /**
        optimize model parameters
        @return the best likelihood
    */
    virtual double optimizeParameters(double gradient_epsilon);

    /**
     * @return TRUE if this is a tree-mixture model, FALSE otherwise
     */
    virtual bool isTreeMixture() { return true; }
    
    /**
        write information
        @param out output stream
    */
    virtual void writeInfo(ostream &out);
    
private:
    IQTreeMix *mixtrees;
};

#endif /* MODELTREEMIXTURE_H_ */

