NLP: Stitch the Torn Wiki

Text blocks which are approximately 500 to 1000 words in length are picked up from N different Wikipedia articles. Every block of text has been picked up from a unique Wikipedia article, about a well known person or place.

Each of these text blocks is split into two parts of roughly equal length.

The first (starting) part obtained after splitting is placed in Set A which will hold all the starting blocks. The second part of the block, is placed in Set B which will contain the second part for all the text fragments which we selected. Both the Sets A and B are shuffled up, and the ordering of elements is lost.

Your task is to identify, for each text fragment (a) in Set A, which is the correct, corresponding text fragment (b) in Set B, such that both a and b were in the same text block initially.
