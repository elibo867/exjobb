
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  301
Structure2: input/casp  Length=  305 (by which all scores are normalized)
Number of residues in common=  301
RMSD of  the common residues=    4.101

TM-score    = 0.8074  (d0= 6.41)
MaxSub-score= 0.5908  (d0= 3.50)
GDT-TS-score= 0.6164 %(d<1)=0.2033 %(d<2)=0.5279 %(d<4)=0.7934 %(d<8)=0.9410
GDT-HA-score= 0.4041 %(d<0.5)=0.0918 %(d<1)=0.2033 %(d<2)=0.5279 %(d<4)=0.7934

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     45.2687397571  -0.5307583009  -0.1659276061   0.8311219258
 2     68.5943018975  -0.7674118604  -0.3220953151  -0.5543768073
 3     42.2763866793   0.3596868951  -0.9320529156   0.0436199493

Superposition in the TM-score: Length(d<5.0)=254  RMSD=  2.29
(":" denotes the residue pairs of distance < 5.0 Angstrom)
TSAIDPVSFSLYAKDFTRFAQELGASFERYGFAVLSDYDLDQARIDAAVDSAKAFFALPVETKKQYAGVKGGARGYIPFGVETAKGADHYDLKEFWHMGRDLPPGHRFRAHMADNVWPAEIPAFKHDVSWLYNSLDGMGGKVLEAIATYLKLERDFFKPTVQDGNSVLRLLHYPPIPKDATRAGAHGDINTITLLLGAEEGGLEVLDRDGQWLPINPPPGCLVINIGDMLERLTNNVLPSTVHRVVNPPPERRGVPRYSTPFFLHFASDYEIKTLQNCVTAENPDRYPESITADEFLQQRL----
  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::         :          :::::::::::              ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: ::::::::::::::::::::   ::::::::::::::::::::::::::::::::::::::::::::::::    ::::::::::::::::::::::::::::::::  :  ::::::::::::    
TSAIDPVSFSLYAKDFTRFAQELGASFERYGFAVLSDYDLDQARIDAAVDSAKAFFALPVETKKQYAGVKGGARGYIPFGVETAKGADHYDLKEFWHMGRDLPPGHRFRAHMADNVWPAEIPAFKHDVSWLYNSLDGMGGKVLEAIATYLKLERDFFKPTVQDGNSVLRLLHYPPIPKDATRAGAHGDINTITLLLGAEEGGLEVLDRDGQWLPINPPPGCLVINIGDMLERLTNNVLPSTVHRVVNPPPERRGVPRYSTPFFLHFASDYEIKTLQNCVTAENPDRYPESITADEFLQQRLREIK
12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345

