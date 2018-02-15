
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  232
Structure2: input/casp  Length=  232 (by which all scores are normalized)
Number of residues in common=  232
RMSD of  the common residues=    3.378

TM-score    = 0.8402  (d0= 5.65)
MaxSub-score= 0.6988  (d0= 3.50)
GDT-TS-score= 0.7166 %(d<1)=0.3707 %(d<2)=0.6681 %(d<4)=0.8448 %(d<8)=0.9828
GDT-HA-score= 0.5097 %(d<0.5)=0.1552 %(d<1)=0.3707 %(d<2)=0.6681 %(d<4)=0.8448

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     55.9156666581  -0.7524324644   0.6567414214  -0.0503596261
 2     13.3104241623   0.6351055059   0.7031230319  -0.3197796091
 3     60.9382981377  -0.1746035021  -0.2725962352  -0.9461526883

Superposition in the TM-score: Length(d<5.0)=200  RMSD=  1.81
(":" denotes the residue pairs of distance < 5.0 Angstrom)
KQKGLAIITGASQGIGAVIAAGLATDGYRVVLIARSKQNLEKVHDEIMRSNKHVQEPIVLPLDITDCTKADTEIKDIHQKYGAVDILVNAAAMFMDGSLSEPVDNFRKIMEINVIAQYGILKTVTEIMKVQKNGYIFNVASADGGIYGSTKFALLGLAESLYRELAPLGIRVTTLCPGWVNTDMAKKAGTPFKDEEMIQPDDLLNTIRCLLNLSENVCIKDIVFEMKKSIIE
  :::::::::::::::::::::::::::::::::::::::::::::::::: :::::::::::::::::::::::::::::::::::::::::::     :::::::::::::::::::::::::::::::::::::::: : ::::::::::::::::::::::::::::::::::::::::::         :::::::::::::::::::  ::::         :  
KQKGLAIITGASQGIGAVIAAGLATDGYRVVLIARSKQNLEKVHDEIMRSNKHVQEPIVLPLDITDCTKADTEIKDIHQKYGAVDILVNAAAMFMDGSLSEPVDNFRKIMEINVIAQYGILKTVTEIMKVQKNGYIFNVASADGGIYGSTKFALLGLAESLYRELAPLGIRVTTLCPGWVNTDMAKKAGTPFKDEEMIQPDDLLNTIRCLLNLSENVCIKDIVFEMKKSIIE
1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012

