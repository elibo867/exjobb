
 *****************************************************************************
 *                                 TM-SCORE                                  *
 * A scoring function to assess the similarity of protein structures         *
 * Based on statistics:                                                      *
 *       0.0 < TM-score < 0.17, random structural similarity                 *
 *       0.5 < TM-score < 1.00, in about the same fold                       *
 * Reference: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710     *
 * For comments, please email to: zhng@umich.edu                             *
 *****************************************************************************

Structure1: input/casp  Length=  255
Structure2: input/casp  Length=  256 (by which all scores are normalized)
Number of residues in common=  255
RMSD of  the common residues=   10.728

TM-score    = 0.7572  (d0= 5.92)
MaxSub-score= 0.6403  (d0= 3.50)
GDT-TS-score= 0.6475 %(d<1)=0.3359 %(d<2)=0.6563 %(d<4)=0.7539 %(d<8)=0.8438
GDT-HA-score= 0.4756 %(d<0.5)=0.1563 %(d<1)=0.3359 %(d<2)=0.6563 %(d<4)=0.7539

 -------- rotation matrix to rotate Chain-1 to Chain-2 ------
 i          t(i)         u(i,1)         u(i,2)         u(i,3)
 1     20.7929646680   0.0089033590   0.7516236721  -0.6595320961
 2     38.9820738663  -0.9911655964   0.0939168352   0.0936503527
 3     36.9458260659   0.1323309892   0.6528717207   0.7458197005

Superposition in the TM-score: Length(d<5.0)=196  RMSD=  1.66
(":" denotes the residue pairs of distance < 5.0 Angstrom)
PSKVALSYSKSLKAPETDSLNLPVDENGYITIFDGETFNGWRGYGKDRVPTKWTIEDGCIKFNGSGGGEAQDGDGGDLIFAHKFKNFELELEWKVAKGSNSGILYLAQEVT-KDKDGNDVLEPIYISAPEYQILDNANHPDAKLGKDNNRQSASLYDMIPAVPQNSKPFGEWNKAKIMVYKGTVVHGQNDENVLEYHLWTKQWTDMLQASKFSEDKWPLAFELLNNCGGENHEGFIGLQDHGDDVWFRNIRVKVLD
                           ::::::::::::::::::::::::::::::::::::::        : ::::::::::::::::::::::::::::::::::::          :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::        :      :::::::::::::::::::::::::
PSKVALSYSKSLKAPETDSLNLPVDENGYITIFDGETFNGWRGYGKDRVPTKWTIEDGCIKFNGSGGGEAQDGDGGDLIFAHKFKNFELELEWKVAKGSNSGILYLAQEVTSKDKDGNDVLEPIYISAPEYQILDNANHPDAKLGKDNNRQSASLYDMIPAVPQNSKPFGEWNKAKIMVYKGTVVHGQNDENVLEYHLWTKQWTDMLQASKFSEDKWPLAFELLNNCGGENHEGFIGLQDHGDDVWFRNIRVKVLD
1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456

