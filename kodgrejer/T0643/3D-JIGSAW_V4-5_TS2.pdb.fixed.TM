
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=   73
Structure2: input/casp  Length=   73 (by which all scores are normalized)
Number of residues in common=   73
RMSD of  the common residues=    5.980

TM-score    = 0.5675  (d0= 3.00)
MaxSub-score= 0.5527  (d0= 3.50)
GDT-TS-score= 0.5822 %(d<1)=0.2603 %(d<2)=0.4932 %(d<4)=0.7123 %(d<8)=0.8630
GDT-HA-score= 0.4110 %(d<0.5)=0.1781 %(d<1)=0.2603 %(d<2)=0.4932 %(d<4)=0.7123

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     27.3868145689   0.7549620380   0.6095388800  -0.2418567241
 2     15.2854143009   0.3255044073  -0.0281516368   0.9451213500
 3     29.7242805655   0.5692795465  -0.7922561701  -0.2196610090

Superposition in the TM-score: Length(d<5.0)= 55  RMSD=  2.21
(":" denotes the residue pairs of distance < 5.0 Angstrom)
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
      :        ::::::::::::::::::::::::::::::::::    ::::::::::::::::::::
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
1234567890123456789012345678901234567890123456789012345678901234567890123

