
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  152
Structure2: input/casp  Length=  152 (by which all scores are normalized)
Number of residues in common=  152
RMSD of  the common residues=   22.362

TM-score    = 0.2241  (d0= 4.59)
MaxSub-score= 0.1382  (d0= 3.50)
GDT-TS-score= 0.1859 %(d<1)=0.0855 %(d<2)=0.1316 %(d<4)=0.1908 %(d<8)=0.3355
GDT-HA-score= 0.1201 %(d<0.5)=0.0724 %(d<1)=0.0855 %(d<2)=0.1316 %(d<4)=0.1908

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     65.6110335421  -0.1755561062  -0.8732475387  -0.4545533959
 2      7.7282778490  -0.8541365564  -0.0944828871   0.5113939060
 3     29.1313078679  -0.4895209869   0.4780289952  -0.7292855978

Superposition in the TM-score: Length(d<5.0)= 29  RMSD=  2.75
(":" denotes the residue pairs of distance < 5.0 Angstrom)
NLYFQGERNYNKWAESYIKYNLSNLKIEKEDLTIYFDNLQVSGNACVSIRKGKQINSFEYIIKFEWLYSKKKEGKDYFGGSVEIPDFSTFSLEENDYAINIERTDESENLRFIYDSILKKEGKEKIKECLKNFQEDLLKHDKNESNKELKIK
                    :                                                                       :::::::::::  ::::  ::::  :   :   :   :  :   ::  ::          
NLYFQGERNYNKWAESYIKYNLSNLKIEKEDLTIYFDNLQVSGNACVSIRKGKQINSFEYIIKFEWLYSKKKEGKDYFGGSVEIPDFSTFSLEENDYAINIERTDESENLRFIYDSILKKEGKEKIKECLKNFQEDLLKHDKNESNKELKIK
12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012

