
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
RMSD of  the common residues=    8.460

TM-score    = 0.5462  (d0= 3.00)
MaxSub-score= 0.5225  (d0= 3.50)
GDT-TS-score= 0.5788 %(d<1)=0.2603 %(d<2)=0.4658 %(d<4)=0.7397 %(d<8)=0.8493
GDT-HA-score= 0.4041 %(d<0.5)=0.1507 %(d<1)=0.2603 %(d<2)=0.4658 %(d<4)=0.7397

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     25.0580190111  -0.4171899679   0.7725111954   0.4787264185
 2     12.2233961991  -0.5978170708  -0.6300174596   0.4956740365
 3     -8.0197337064   0.6845197445  -0.0794005898   0.7246573437

Superposition in the TM-score: Length(d<5.0)= 56  RMSD=  2.37
(":" denotes the residue pairs of distance < 5.0 Angstrom)
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
            ::::::::::::::::::::::::::::::::::::     ::::::::::::::::::::
HHSHMLPPEQWSHTTVRNALKDLLKDMNQSSLAKECPLSQSMISSIVNSTYYANVSAAKCQEFGRWYKHFKKT
1234567890123456789012345678901234567890123456789012345678901234567890123

