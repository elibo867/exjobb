
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=   99
Structure2: input/casp  Length=   99 (by which all scores are normalized)
Number of residues in common=   99
RMSD of  the common residues=   13.542

TM-score    = 0.2493  (d0= 3.63)
MaxSub-score= 0.1889  (d0= 3.50)
GDT-TS-score= 0.2500 %(d<1)=0.1414 %(d<2)=0.1717 %(d<4)=0.2828 %(d<8)=0.4040
GDT-HA-score= 0.1793 %(d<0.5)=0.1212 %(d<1)=0.1414 %(d<2)=0.1717 %(d<4)=0.2828

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1      2.0037028283  -0.4919707078   0.4241641202  -0.7602957463
 2   -107.5417502204  -0.1977407449   0.7960252317   0.5720510715
 3     12.7164357616   0.8478581370   0.4317738177  -0.3077465675

Superposition in the TM-score: Length(d<5.0)= 30  RMSD=  2.70
(":" denotes the residue pairs of distance < 5.0 Angstrom)
KENKLDYIPEPMDLSLVDLPESLIQLSERIAENVHEVWAKARIDEGWTYGEKRDDIHKKHPCLVPYDELPEEEKEYDRNTAMNTIKMVKKLGFRIEKED
                                  :::                               ::::::::::::::::::::::: : :::  
KENKLDYIPEPMDLSLVDLPESLIQLSERIAENVHEVWAKARIDEGWTYGEKRDDIHKKHPCLVPYDELPEEEKEYDRNTAMNTIKMVKKLGFRIEKED
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789

