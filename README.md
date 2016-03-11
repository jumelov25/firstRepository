In this repository you can find:

    1. Solution to an NLP problem from Hackerrank:
      1.1. Algorithm implementation
      1.2. Brief Report
      
    2. Four examples of Codes written in C (by me and by my team here at VOCALIZE)
      1.1 viterbiSearch (viterbi search for TTS unit selection)
      1.2 hmmTagger (PoS Tagger based on Hidden Markov Models)
      1.3 vad_Streaming (Voice Activity Detection for Streaming operation)
      1.4 vad_Batch (Voice Activity Detection for Batch operation)


Here is the description problem of the NLP problem from Hackerrank. The name of the problem is: Stitch the Torn Wiki.  

You can find bellow the original description of the problem that was extracted from: https://www.hackerrank.com/challenges/stitch-the-torn-wiki 

Text blocks which are approximately 500 to 1000 words in length are picked up from N different Wikipedia articles. Every block of text has been picked up from a unique Wikipedia article, about a well known person or place.

Each of these text blocks is split into two parts of roughly equal length.

The first (starting) part obtained after splitting is placed in Set A which will hold all the starting blocks. The second part of the block, is placed in Set B which will contain the second part for all the text fragments which we selected. Both the Sets A and B are shuffled up, and the ordering of elements is lost.

Your task is to identify, for each text fragment (a) in Set A, which is the correct, corresponding text fragment (b) in Set B, such that both a and b were in the same text block initially.

Input Format 
An Integer N on the first line. This is followed by 2N+1 lines. 
Text fragments (numbered 1 to N) from Set A, each on a new line (so a total of N lines). 
A separator with five asterisk marks "*" which indicates the end of Set A and beginning of Set B. 
Text fragments (numbered 1 to N) from Set B, each on a new line (so a total of N lines).

Output Format
N lines, each containing one integer. 
The i-th line should contain an integer j such that the i-th element of Set A and the j-th element of Set B are a pair, i.e., both originally came from the same block of text/Wikipedia article.

Constraints 
1 <= N <= 100 
No text fragment will have more than 10000 characters.

Explanation 
The first, second and third text fragment of Set A are about Delhi, Seattle and Martin Luther respectively. 
In set B, the paragraph on Delhi, is the second text fragment. 
The paragraph on Seattle is the first text fragment in Set B. 
The paragraph on Martin Luther Kind is the third text fragment in Set B. 
So, the expected output is 2, 1, 3 respectively.

Scoring 
A sample test case with twenty paragraphs is provided to you when you Compile and Test. 
Extensive training data is not required for this challenge. The weightage for a test case will be proportional to the number of tests (Articles) which it contains. This works out to a ratio of 1:2 (Sample Test: Hidden Test). 
Score = M * (C)/N Where M is the Maximum Score for the test case. 
C = Number of correct answers in your output. 
N = Total number of Wikipedia Articles (which were split into 2N fragments and divided into Set A and Set B respectively).

