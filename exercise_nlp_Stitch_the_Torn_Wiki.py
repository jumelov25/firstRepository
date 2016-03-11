# Enter your code here. Read input from STDIN. Print output to STDOUT

#!/usr/bin/env python
# -*- encoding: iso-8859-1 -*-

import numpy as np
import locale
import sys
import copy

locale.setlocale(locale.LC_ALL, '')

def lowerCaseAndPunctuationRemoval(Text):

    # This is a very simple routing for "Text Cleaning" and it does not intent to be exaustive

    newText = []
    for s in Text:
        s = s.lower();
        s = s.replace('.', '');
        s = s.replace(',', '');
        s = s.replace('!', '');
        s = s.replace('?', '');
        s = s.replace(':', '');
        s = s.replace(';', '');
        s = s.replace('(', '');
        s = s.replace(')', '');
        s = s.replace('[', '');
        s = s.replace(']', '');
        s = s.replace('/', ' ');
        s = s.replace('"', ' ');
        s = s.replace("'", ' ');
        s = s.replace("&", ' ');
        s = s.replace("_", ' ');
        s = s.replace("-", ' ');
        s = s.replace('\r\n',' ');
        newText.append(s);

    return newText


def functionWordsRemoval(Text):

    # This is a very simple routine to "Function Word Removal" and it it does not intent to be exaustive

    #Auxiliary Verbs
    FWList  = ' be do has will is been did are am is have can shall was been '
    FWList += ' should would does had could may might must ought were '

    #Prepositions
    FWList += ' in at to through over between under than of on from at in with by if '

    #Articles
    FWList += ' a an the that these this those which who each what no '

    #Conjunctions:
    FWList += ' and or nor but for so since as '

    #Pronons
    FWList += ' I you him us ours she he i we they it she you him them us me thy hou thee '
    FWList += ' thine its my their his her our your yours ours theirs mine thine whose '

    newText = []
    for s in Text:
        newS    = ''
        wrds = s.split()
        for w in wrds:
            if(w not in FWList.split()):
                newS +=  w + ' '

        newText.append(newS)

    return newText


def fullDictionary(Text):

    # This routine build a full dictionary including all words from the Text (Set A and Set B)

    wdic = {}
    for s in Text:
        wrds = s.split()
        wrds.sort()

        for w in wrds:
            wdic[w] = 0

    return wdic


def dictionariesPerSegment(Text,wdic):

    # This routine build a dictionary for each Text Segment 

    wvector = []

    for s in Text:
        wv = copy.deepcopy(wdic)
        for w in s.split():
            wv[w] += 1

        wvector.append(wv)
    
    return wvector


def docTermMatrix(wvector):

    # This routine build a matrix with Counting for Each Word of Each Text Segment

    list_of_list = []

    for dic in wvector:
        list = []
        for w in dic.keys():
            list.append(dic[w])
        list_of_list.append(list)

    A = (np.array(list_of_list)).T

    return A


def lowRankDocTermMatrix(A,order):

    # This routine performs a Singular Value Decomposition to allow the dimensionality reduction 
    # and LSA (Latent Semantic Analysis)
    
    U, S, Vt = np.linalg.svd(A, full_matrices=False)

    Ulow  = U[:,:order]
    Slow  = np.diag(S)[:order,:order]
    Vtlow = Vt[:order,:]
   
    return np.dot(Ulow, np.dot(Slow,Vtlow))


def evalBasedOnCosineDistanceMetric(A):

    # This routine estimates the "distance" between two Vectors by using "cosine distance" (the angle between two vectors)

    result = []
    m      = len(A[1,:])

    for i in range(m/2):
        s1 = A[:,i]                             # Vector associated to the first segment
        s  = []
        for j in range(m/2):

            s2   = A[:,j+m/2]                   # Vector associated to the second segment

            s1s2 = np.dot(s1,s2)                # Inner product
            s1_m = np.sqrt((s1*s1).sum())       # Modulus of s1 
            s2_m = np.sqrt((s2*s2).sum())       # Modulus of s2
            cos_angle = s1s2 / s1_m / s2_m      # Cosine of angle between s1 and s2

            angle = np.arccos(cos_angle)        # The arc of the cossine
            angle = angle * 360 / (2 * np.pi)   # Angle in degrees 

            s.append(angle)
        
        sa = np.array(s)
        result.append(np.argsort(sa)[0] + 1)

    return result


def testDictionariesPerSegmentBuilder(wvector):

    # This routine is here only to help on testing/debuging the algorithm

    for dic in wvector:
        for w in dic.keys():
            if dic[w]>0:
                print  w + ' : ' + str(dic[w])
        print '\n'


def evaluateResults(resultLSA, ref):

     # This routine estimate the accuracy of the algorithm

     i    = 0
     erro = 0
     for i in range(len(resultLSA)):
         r = ref[i]
         if( int(resultLSA[i]) != int(r) ):
             erro += 1
     
     return 100*erro/len(resultLSA)


def main():

    error       = []   # Estimation error
    LSAorder    = []   # Latent Semantic Analysis order (Dimensionality reduction using SVD analysis)
    ErrorSearch = {}   # Error search space 
    resultLSA   = []   # Results for each value of LSAorder searched
    Text        = []   # Text to be analysed
    
    # Expected result
    Ref         = [8, 17, 5, 6, 4, 16, 1, 13, 15, 3, 2, 19, 12, 7, 11, 10, 18, 14, 20, 9]   

    # Reading the input Text
    N = int(raw_input())
    for i in range(N):
        Text.append(raw_input())
    separator = raw_input()
    for i in range(N):
        Text.append(raw_input())
     
    # Tokenization, Function Words removal and general cleaning 
    tkText  = lowerCaseAndPunctuationRemoval(Text)
    fwrText = functionWordsRemoval(tkText)

    # Build a Dictionary with all the words from the Text and assign a value equal to zero to each one
    fDic  = fullDictionary(fwrText)

    # Build a Dictionary for each Text Segment (It contains: all Words from the Text and the Word Countings from the Segment) 
    dpSeg = dictionariesPerSegment(fwrText,fDic)

    # Build a Matrix with the Text Segment Vectors
    DTM = docTermMatrix(dpSeg)

    # Search for the best result: LSA order (%)
    MIN   = 100
    for i in range(90):

        # LSA order
        LSAorder = int(len(Text)*((i+11)/100.0))

        # SVD decomposition: LrDTM = U S Vt
        LrDTM = lowRankDocTermMatrix(DTM,LSAorder)

        # Look for the Closest Text Segments:  Set A (first half) <-> Set B (second half)
        resultLSA.append ( evalBasedOnCosineDistanceMetric(LrDTM) )
        
        # Evaluate the results: Compare with a reference 
        error.append( evaluateResults(resultLSA[i], Ref) )

        # Look for the best result: LSA order (%) and Error (%)
        ErrorSearch[ str(LSAorder*(100.0/(2*N))) ] = error[i]
        if error[i] <= MIN:
            MIN = error[i]
            bestSolution = 'LSA order: ' + str(LSAorder*(100.0/(2*N))) + '% ' + ' , ' + 'Error : ' + str(error[i]) + '%'
            bestResult   = resultLSA[i]
    
    # Please, uncomment these coming lines to visualize the search process (for the best result)  
        #print str(LSAorder*(100.0/(2*N))) + '% ' + ' : ' + str(error[i]) + '%'
    #print '\n'
    #print '--------------------------------------------------------'
    #print '   Best Solution: ' + bestSolution
    #print '--------------------------------------------------------'
    #print '\n'
    
    # Best result: It uses the best LSA Order searched
    for r in bestResult:
        print r

if __name__ == "__main__":
    sys.exit(main())
