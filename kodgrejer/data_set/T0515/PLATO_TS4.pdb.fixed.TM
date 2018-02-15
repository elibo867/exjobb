
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  348
Structure2: input/casp  Length=  348 (by which all scores are normalized)
Number of residues in common=  348
RMSD of  the common residues=  163.342

TM-score    = 0.1375  (d0= 6.79)
MaxSub-score= 0.0594  (d0= 3.50)
GDT-TS-score= 0.0826 %(d<1)=0.0374 %(d<2)=0.0517 %(d<4)=0.1006 %(d<8)=0.1408
GDT-HA-score= 0.0524 %(d<0.5)=0.0201 %(d<1)=0.0374 %(d<2)=0.0517 %(d<4)=0.1006

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     40.9060363641  -0.4599267082   0.7671747714  -0.4471132890
 2      3.1668782816   0.6236188308  -0.0793766494  -0.7776881775
 3     50.1011969562  -0.6321131046  -0.6365078300  -0.4419171929

Superposition in the TM-score: Length(d<5.0)= 20  RMSD=  2.71
(":" denotes the residue pairs of distance < 5.0 Angstrom)
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
                                        :                                                         ::::::::::::::::                              :::                                                                                                                                                                                                         
MIETPYYLIDKAKLTRNMERIAHVREKSGAKALLALKCFATWSVFDLMRDYMDGTTSSSLFEVRLGRERFGKETHAYSVAYGDNEIDEVVSHADKIIFNSISQLERFADKAAGIARGLRLNPQRLGEWDVPKVERVMDRINGFMIHNNCENKDFGLFDRMLGEIEERFGALIARVDWVSLGGGIHFTGDDYPVDAFSARLRAFSDRYGVQIYLEPGEASITKSTTLEVTVLDTLYNKNLAIVDSSIEAHMLDLLIYRETAKVLPNEGSHSYMICGKSCLAGDVFGEFRFAEELKVGDRISFQDAAGYTMVKKNWFNGVKMPAIAIRELDGSVRTVREFTYADYEQSLS
123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678

