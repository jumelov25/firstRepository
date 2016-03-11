/**
 * VOCALIZE - Speech and Language Technology Solutions
 *
 *      "Viterbi Search for TTS" 
 *
 * @author VOCALIZE Team
 * @version 9.7
 *
 * (All Rights Reserved)
 */

/**
 * Standard include files
 */
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cmath>
#include <assert.h>

#include <log4cxx/logger.h>


static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("unit_selection.viterbi_search.viterbi"));
static log4cxx::LoggerPtr logger_valuation(log4cxx::Logger::getLogger("unit_selection.viterbi_search.viterbi.valuation"));


using namespace std;

/**
 * Viterbi Search specific include files
 */
#include "unit_selection/viterbi_search/viterbi.h"

/**
 * Print matrix MD ???
 */
void printmd (float **mat, int m, int n/*, char *str*/) {
    int i, j;

    printf("\n");
    for (i = 0; i < m; i++) {
        printf( "linha %d total linhas %d totalcolunas %d\n", i, m, n );
        for (j = 0; j < n; j++) {
            printf ("%e ", mat[i][j]);
        }
        printf ("\n");
    }
    printf ("\n");
}


/**
 * Print matrix MI ???
 */
void printmi (int **mat, int m, int n) {
    int i, j;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf ("%2d ", mat[i][j]);
        }
        printf ("\n");
    }
    printf ("\n");
}


/**
 * Maximum in a vector: It return the maximum of an Array of Integer
 */
int maxiv (int *vector, int n) {
    int i, in, max;

    in = 0;
    max = vector[0];
    for (i = 1; i < n; i++) {
        if (vector[i] > max) {
            max = vector[i];
            in = i;
        }
    }
    return max;
}


/**
 * Maximum in a vector: It return the maximum of an Array of Integer and the position of this maximum
 */
void maxdv (float *vector, int n, float *max, int *in) {
    int i;

    *in = 0;
    *max = vector[0];
    for (i = 1; i < n; i++) {
        if (vector[i] > *max) {
            *max = vector[i];
            *in = i;
        }
    }
}


/**
 * Allocating the Array of Array of integer Indexes. It is used for backtracking
 */
int **allocMi (int m, int n) {
    int i, **mat;

    mat = (int **) calloc (m, sizeof (int *));
    if (mat == NULL)
        printf ("Erro: calloc(mat)\n");
    for (i = 0; i < m; i++) {
        mat[i] = (int *) calloc (n, sizeof (int));
        if (mat[i] == NULL)
            printf ("Erro: calloc(mat[])\n");
        memset (mat[i], -1, n*sizeof(int));
    }
    return mat;
}


/**
 * Deallocating the Array of Array of integers
 */
void deallocMi (int **mat, int m, int n) {
    int i;

    for (i = 0; i < m; i++)
        free (mat[i]);
    free (mat);
}


/**
 * Allocating the Array of Array of Maximum Likelihoods
 */
float **allocMd (int m, int n) {
    int i;
    float **mat;

    mat = (float **) calloc (m, sizeof (float *));
    if (mat == NULL)
        printf ("Erro: calloc(mat)\n");
    for (i = 0; i < m; i++) {
        mat[i] = (float *) calloc (n, sizeof (float));
        if (mat[i] == NULL)
            printf ("Erro: calloc(mat[])\n");
    }

    return mat;
}


/**
 * Deallocating the Array of Array of Maximum Likelihoods
 */
void deallocMd (float **mat, int m, int n) {
    int i;

    for (i = 0; i < m; i++)
        free (mat[i]);
    free (mat);
}


/**
 * Allocating the Array of Array of Maximum Likelihoods
 */
float **allocCP (int nCol, int *nEst) {
    float **cp;
    int i, aux;

    aux = maxiv (nEst, nCol);
    cp  = (float **) calloc (aux, sizeof (float *));
    if (cp == NULL)
        printf ("Erro: calloc(cp)\n");
    for (i = 0; i < aux; i++) {
        cp[i] = (float *) calloc (nCol, sizeof (float));
        if (cp[i] == NULL)
            printf ("Erro: calloc(cp[])\n");
    }

    return cp;
}


/**
 * Deallocating the Array of Array of Maximum Likelihoods
 */
void deallocCP (float **cp, int nCol, int *nEst) {
    int i;

    for (i = 0; i < maxiv (nEst, nCol); i++)
        free (cp[i]);
    free (cp);
}


/**
 * Viterbi Search
 */
UttViterbiOutput* viterbi (UttViterbiInputLangFeats* inputLangFeat,
                           UttViterbiInputLeafSelection* inputCluster,
                           UttViterbiInputCCost* inputCCost, 
                           UttViterbiInputPPCost* inputPPCost,
                           UttViterbiModel* model)
{
    string classes;

    int     i, j, k, in;
    int     **MatIn;
    float   sum, max;
    float   beta;
    float   **ProbCP;
    float   **ProdProb;
    float   **MaxProdProb;
    int     nCol;

    /* // Unused variables
    bool boundarySyll;
    bool boundaryWord;
    bool boundaryPhW;
    bool boundaryPhP;
    bool boundaryIP;
    */

    UttViterbiOutput* uttOut = new UttViterbiOutput();
    nCol    = (int)inputCluster->leafs.size();
    ProbCP  = allocCP (nCol, inputPPCost->NExem);

    /// Setting the Phonetic Prosodic Cost
    for (i = 0; i < nCol; i++){
        string buffer = "";
        char valueTmp[20];
        LOG4CXX_TRACE(logger_valuation, "viterbi_num_samples = " << inputPPCost->NExem[i]);
        for (j = 0; j < inputPPCost->NExem[i]; j++) {
            ProbCP[j][i] = inputPPCost->prob[i][j];    // Phonetic prosodic costs given by regression models.
            sprintf(valueTmp, "%f", ProbCP[j][i]);
            buffer = buffer + string(valueTmp) + " ";
        }
        LOG4CXX_TRACE(logger_valuation, "viterbi_prob_input = " << buffer);
    }

    /// From the first to the previous last Synthesis Unit
    for (i = 0; i < nCol-1; i++) {

        /// Linguistic Class associated to the PolyPhones to be Concatenated
        classes = inputCluster->leafs[i].rightClass + "_" + inputCluster->leafs[i+1].leftClass;

        /// Linguistic Class associated to the PolyPhones to be Concatenated
        for (j = 0; j < inputCCost->CCM[i].n; j++) {

            max = inputCCost->CCM[i].p[0][j];
            for (k = 1; k < inputCCost->CCM[i].m; k++) {
                if (inputCCost->CCM[i].p[k][j] > max) {
                    max = inputCCost->CCM[i].p[k][j];
                }
            }

            /// Setting the cCost selective pressure weight
            beta = (float)(1.0*(*model->ccostSelectivePressure)[classes] / ((float)(max)+0.000000000000000000001));

            /// Transforming the cCost distance into cCost Probabilities using the SoftMax function
            sum = 0.0;
            for (k = 0; k < inputCCost->CCM[i].m; k++) {
                inputCCost->CCM[i].p[k][j] = (float)(expf(-1*beta*inputCCost->CCM[i].p[k][j]) + 0.000000000001/(inputCCost->CCM[i].m));
                sum = sum + inputCCost->CCM[i].p[k][j];
         }

            for (k = 0; k < inputCCost->CCM[i].m; k++) {
            inputCCost->CCM[i].p[k][j] /= (float)(sum + 0.000000000000000001);
                inputCCost->CCM[i].p[k][j]  = log10(inputCCost->CCM[i].p[k][j]); 
            }

        }
    }

    /// Evaluating the Maximum Likelihood Matrix: It uses cCost and pPCost (in log10 domain) simultaneouly
    MatIn       = allocMi (maxiv (inputPPCost->NExem, nCol), nCol-1);
    MaxProdProb = allocCP (nCol, inputPPCost->NExem);
    for (i = 0; i < maxiv (inputPPCost->NExem, nCol); i++){
        MaxProdProb[i][0] = ProbCP[i][0];
    }

    for (i = 0; i < nCol-1; i++) {
        ProdProb = allocMd (inputCCost->CCM[i].m, inputCCost->CCM[i].n);

        for (j = 0; j < inputCCost->CCM[i].n; j++) {
            for (k = 0; k < inputCCost->CCM[i].m; k++) {
                ProdProb[k][j] = MaxProdProb[j][i] + inputCCost->CCM[i].p[k][j] + ProbCP[k][i+1];
            }
        }

        for (j = 0; j < inputPPCost->NExem[i+1]; j++){
            maxdv (ProdProb[j], inputPPCost->NExem[i], &MaxProdProb[j][i+1], &MatIn[j][i]);
        }

        deallocMd (ProdProb, inputCCost->CCM[i].m, inputCCost->CCM[i].n);
    }

    /// Identifying the Winner Sample for each Synthesis Unit
    in  = 0;
    max = MaxProdProb[0][nCol-1];
    for (i = 1; i < inputPPCost->NExem[nCol-1]; i++) {
        if (MaxProdProb[i][nCol-1] > max) {
            max = MaxProdProb[i][nCol-1];
            in  = i;
        }
    }

    /// Performing the Backtracking Decodification
    uttOut->seq         = (int *)calloc(nCol, sizeof(int));
    uttOut->seqlen      = nCol;
    uttOut->seq[nCol-1] = in;

    stringstream vitSeqStr; 
    for (i = nCol-2; i >= 0; i--){
        uttOut->seq[i] = MatIn[uttOut->seq[i+1]][i];
        vitSeqStr << uttOut->seq[i] << " ";
    }
    LOG4CXX_DEBUG(logger, "Viterbi.Sequence = " << vitSeqStr.str());

    float sumPPcost=0.0;
    stringstream PPcostStr; 
    for (int pp=0; pp<nCol; ++pp){
        float ff = 1 - pow(10,ProbCP[uttOut->seq[pp]][pp]);
        PPcostStr << ff << " ";
        sumPPcost = sumPPcost + ff;
    }
    LOG4CXX_TRACE(logger, "Viterbi.WinSeq_PPcost:" << PPcostStr.str());
    LOG4CXX_TRACE(logger, "Viterbi.WinSeq_PPcost_SUM:" << sumPPcost);

    float sumCCcost=0.0;
    stringstream CCcostStr; 
    for (int pp=0; pp<nCol-1;++pp){
        float ff = 1 - pow(10,inputCCost->CCM[pp].p[ uttOut->seq[pp+1] ][ uttOut->seq[pp] ]);
        CCcostStr << ff << " ";
        sumCCcost = sumCCcost + ff;
    }
    LOG4CXX_TRACE(logger, "Viterbi.WinSeq_CCcost:" << CCcostStr.str());
    LOG4CXX_TRACE(logger, "Viterbi.WinSeq_CCcost_SUM:" << sumCCcost);

    /// Deallocating all the memory used 
    deallocCP (MaxProdProb, nCol, inputPPCost->NExem);
    deallocCP (ProbCP,      nCol, inputPPCost->NExem);
    deallocMi (MatIn, maxiv (inputPPCost->NExem, nCol), nCol-1);

    return uttOut;
}

