
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  135
Structure2: input/casp  Length=  135 (by which all scores are normalized)
Number of residues in common=  135
RMSD of  the common residues=   13.519

TM-score    = 0.2693  (d0= 4.32)
MaxSub-score= 0.1944  (d0= 3.50)
GDT-TS-score= 0.2444 %(d<1)=0.1481 %(d<2)=0.2074 %(d<4)=0.2519 %(d<8)=0.3704
GDT-HA-score= 0.1722 %(d<0.5)=0.0815 %(d<1)=0.1481 %(d<2)=0.2074 %(d<4)=0.2519

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1      9.1234277996   0.2033974675  -0.9789977844   0.0138855458
 2     -4.2609183657  -0.6873065571  -0.1326663464   0.7141493801
 3     -1.4981179590  -0.6973085162  -0.1547998020  -0.6998555955

Superposition in the TM-score: Length(d<5.0)= 33  RMSD=  1.97
(":" denotes the residue pairs of distance < 5.0 Angstrom)
MPLEWRAGASSDEINAIIRAVYRQVLGNDYVMSTERLTSAESLLRGGEISVRDFVRAVALSELYREKFFHNNAHNRFIELNFKHLLGRAPYDQAEVAAHAATYHSHGYDADINSYIDSAEYTESFGDNVVPYFRG
  :::::::::::::::::::::::::::::: :     ::                                                                                              
MPLEWRAGASSDEINAIIRAVYRQVLGNDYVMSTERLTSAESLLRGGEISVRDFVRAVALSELYREKFFHNNAHNRFIELNFKHLLGRAPYDQAEVAAHAATYHSHGYDADINSYIDSAEYTESFGDNVVPYFRG
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345

