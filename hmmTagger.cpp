/**
 * VOCALIZE - Speech and Language Technology Solutions
 *
 *      "hmmTagger for TTS Unit Selection" 
 *
 * @author VOCALIZE Team
 * @version 9.7
 *
 * (All Rights Reserved)
 */

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <limits>

#include <log4cxx/logger.h>

#include "frontend/hmm_tagger/Transition.h"
#include "frontend/hmm_tagger/Tag.h"
#include "frontend/hmm_tagger/Lexicon.h"
#include "frontend/hmm_tagger/WordTree.h"
#include "frontend/hmm_tagger/HMMTagger.h"
#include "frontend/hmm_tagger/Train.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("frontend.hmm_tagger.HMMTagger"));

using namespace std;

// Definição dos vetores R^3 e K^3 usados no algoritmo
typedef map< int, map< string, map < string, double > > > DeltaMap;
typedef map< int, map< string, map < string, string > > > PsiMap;

/**
 * Construtor. Instancia as variáveis locais.
 * @param tr objeto contendo matrizes n-grama
 * @param tv objeto contendo lista de etiquetas usadas
 * @param lex objeto contendo entradas lexicais e ambiguidades sintáticas
 * @param sr objeto contendo árvore de sufixos para palavras desconhecidas
 */
HMMTagger::HMMTagger( Transition & tr, TagPool & tv, Lexicon & lex, WordTree & st ):
        trans_( tr ), lexicon_( lex ), suffixTree_( st ) {
   tags_ = tv.getTags();
}

/**
 * Destrutor
 */
HMMTagger::~HMMTagger() {
}

/**
 * Retorna uma lista com todas as possíveis etiquetas para
 * uma palavra e suas respectivas probabilidades.
 * Se a palavra não consta no léxico, ela é estimada a
 * partir da árvore de sufixos.
 * 
 * @param word palavra a ser consultada
 * @return mapa com etiquetas possível e suas probabilidades
 */
map< string, Lexicon::Form > HMMTagger::getTagsForWord( string & word ) {

    if( lexicon_.existOnLexicon( word ) ) {
        return lexicon_.getForms( word );
    }
    else {
        // Repensar em como utilizar o retorno das palavras desconhecidas
        // Montar um mapa de Forms não é rápido.... (TODO)

        string inverted = word;
        Train::invertWord( inverted );
        map< string, double > & prob = suffixTree_.recoverData( inverted );
        map< string, Lexicon::Form > mapForms;
        
        map< string, double >::iterator i;
        for( i = prob.begin(); i != prob.end(); i++ ) {
            Lexicon::Form f;
            f.probability = i->second;
            f.count = 0;
            mapForms[ i->first ] = f;
        }
        
        return mapForms;
    }
}

/**
 * Retorna o valor de transição de uma dada sequencia de três
 * etiquetas, temporalmente equivalente a ordem dos parametros.
 * Atualmente implementada com n-gramas.
 *
 * @param word_i primeira etiqueta da sequencia
 * @param word_j segunda etiqueta da sequencia
 * @param word_k terceira etiqueta da sequencia
 * @return probabilidade de transição pela sequencia
 */
double HMMTagger::transitionValue( const string & tag_i, const string & tag_j, const string & tag_k ) {
    if( trans_.doesExistOnTrigram( tag_i, tag_j, tag_k ) ) {
        return trans_.getTrigram( tag_i, tag_j, tag_k );
    }
    if( trans_.doesExistOnBigram( tag_j, tag_k ) ) {
        return trans_.getBigram( tag_j, tag_k ) + log( 0.5 );
    }
    return trans_.getUnigram( tag_k ) + log( 0.25 ) ;
}

/**
 * Algoritmo de viterbi para HMM de segunda ordem, utilizando
 * n-gramas para calculo de transição de estados e algumas simplificações.
 * 
 * @param wordSequence um vetor contendo as palavras de uma sentença
 * @return um vetor de etiquetas co-indexado com as palavras de entrada
 */
vector< string > HMMTagger::viterbi( vector< string > & wordSequence ) {

    int numWords = wordSequence.size();
    const double MINUS_INF = - numeric_limits< double >::infinity();
    DeltaMap delta;
    PsiMap psi;

    map< string, Lexicon::Form > probTagWord;
    map< string, Lexicon::Form >::iterator ptwIt;
    double coefB = 0.0;
    vector< string >::iterator i, j, k;

    LOG4CXX_TRACE( logger,  "----- Initialization -----" );

    /*
     * First step: Initialization
     *
     * Como a primeira palavra da sentença é o marcador <s>, as probabilidades
     * são 1 para delta(i,j) onde Tj = <s>, e zero o contrário.
     */

    for( j = tags_.begin(); j != tags_.end(); j++ ) {
        if( *j == Train::START_TAG ) {
            for( i = tags_.begin(); i != tags_.end(); i++ ) {
                delta[ 0 ][ *i ][ *j ] = 0; 
             }
        }
    }

    LOG4CXX_TRACE( logger, "----- Maximum Probability -----" );

    /*
     * Second step: Calculate maximum probability for word Wt
     */
    int t;
    for( t = 1; t < numWords; t++ ) {

        probTagWord = getTagsForWord( wordSequence[ t ] );
        LOG4CXX_TRACE( logger, "The word is << " << wordSequence[ t ] << " >>" );

        for( k = tags_.begin(); k != tags_.end(); k++ ) {

            bool findMaximum = true;

            // If it is the last word, then only the tag </s> has probability different from zero
            if( t == numWords - 1 ) {
                if( *k != Train::END_TAG ) {
                    findMaximum = false;
                    coefB = MINUS_INF;
                }
                else {
                    coefB = 0;
                }
            }
            else {
                // Retrieving tag for current word. If not existing, probability is zero...
                ptwIt = probTagWord.find( *k );
                if( ptwIt == probTagWord.end() ) {
                    findMaximum = false;
                    coefB = MINUS_INF;
                }
                // ...or else, calculate the probability P( W | Tk ) = P( Tk | W ) / P( Tk )
                else {
                    coefB = ptwIt->second.probability - trans_.getUnigram( *k );
                }
            }

            for( j = tags_.begin(); j != tags_.end() && findMaximum; j++ ) {

                double bestProb = MINUS_INF;
                string bestTag = "";

                for( i = tags_.begin(); i != tags_.end(); i++ ) {  

                    if( delta[ t - 1][ *i ].find( *j ) == delta[ t - 1 ][ *i ].end() ) {
                        continue;
                    }

                    double prob = delta[ t - 1 ][ *i ][ *j ] + transitionValue( *i, *j, *k );

                    if( prob > bestProb ) {                        
                        bestProb = prob;
                        bestTag = *i;
                    }
                }

                if( bestTag != "" ) {
                    delta[ t ][ *j ][ *k ] = bestProb + coefB;
                    psi[ t ][ *j ][ *k ] = bestTag;
 
                }
            }
        }
    }

    LOG4CXX_TRACE( logger, "----- Retrieve Process -----" );

    /*
     * Third step: Retrieve tags of best probability
     */
    vector< string > bestTag( numWords );

    // Retrieving T + 1, where T = numWords
    double bestProb = MINUS_INF;
    for( i = tags_.begin(); i != tags_.end(); i++ ) {
        for( j = tags_.begin(); j != tags_.end(); j++ ) {
            if( delta[ numWords - 1 ][ *i ].find( *j ) != delta[ numWords - 1 ][ *i ].end() ) {
                double prob = delta[ numWords - 1 ][ *i ][ *j ];
                if( prob > bestProb ) {
                    bestTag[ numWords - 1 ] = *j;
                    bestProb = prob;
                }
            }
        }
    }
    LOG4CXX_TRACE( logger, "bestTag[" << (numWords - 1) << "] = " << bestTag[ numWords - 1] );

    // Retrieving T
    bestProb = MINUS_INF;
    for( i = tags_.begin(); i != tags_.end(); i++ ) {
        for( j = tags_.begin(); j != tags_.end(); j++ ) {
            if( delta[ numWords - 1 ][ *i ].find( *j ) != delta[ numWords - 1 ][ *i ].end() ) {
                double prob = delta[ numWords - 1 ][ *i ][ *j ];
                if( prob > bestProb ) {
                    bestTag[ numWords - 2 ] = *i;
                    bestProb = prob;
                }
            }
        }
    }
    LOG4CXX_TRACE( logger, "bestTag[" << (numWords - 2) << "] = " << bestTag[ numWords - 2] );

    // Retrieving the remaining tags
    for( t = numWords - 3; t >= 0 ; t-- ) {
        bestTag[ t ] = psi[ t + 2 ][ bestTag[ t + 1 ] ][ bestTag[ t + 2] ];
        LOG4CXX_TRACE( logger, "bestTag[" << t << "] = " << bestTag[ t ] );
    }

    LOG4CXX_TRACE( logger, "----- Algorithm Concluded -----" );

    return bestTag;
}

/**
 * Aplica o algoritmo de viterbi a uma sequencia de palavras
 * que não contém as marcas de início e fim de sentença
 * @param wordSequence sequencia de palavras sem marcas de inicio e fim de sentença
 * @return vetor correspondente ao de entrada, com as etiquetas morfossintáticas
 */
vector< string > HMMTagger::tag( vector< string > & wordSequence ) {
    wordSequence.insert( wordSequence.begin(), Train::START_TAG );
    wordSequence.push_back( Train::END_TAG );

    vector< string > wordWithTags = viterbi( wordSequence );
    
    wordSequence.pop_back();
    wordSequence.erase( wordSequence.begin() );

    wordWithTags.pop_back();
    wordWithTags.erase( wordWithTags.begin() );
    return wordWithTags;
}

/**
 * Retorna o objeto que encapsula o léxico.
 * @return léxico usado pelo tagger
 */
Lexicon * HMMTagger::getLexicon() {
    return &lexicon_;
}
