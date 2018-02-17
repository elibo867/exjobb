
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
RMSD of  the common residues=    8.467

TM-score    = 0.5457  (d0= 3.00)
MaxSub-score= 0.5218  (d0= 3.50)
GDT-TS-score= 0.5822 %(d<1)=0.2740 %(d<2)=0.4658 %(d<4)=0.7397 %(d<8)=0.8493
GDT-HA-score= 0.4075 %(d<0.5)=0.1507 %(d<1)=0.2740 %(d<2)=0.4658 %(d<4)=0.7397

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     25.0104587406  -0.4084490387   0.7676356434   0.4938672917
 2     12.3064320955  -0.6106648565  -0.6319469273   0.4772122318
 3     -8.0821086142   0.6784230360  -0.1066705215   0.7268862249

Superposition in the TM-score: Length(d<5.0)= 56  RMSD=  2.38
(":" denotes the residue pairs of distance < 5.0 Angstrom)
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
            ::::::::::::::::::::::::::::::::::::     ::::::::::::::::::::
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
1234567890123456789012345678901234567890123456789012345678901234567890123

