# guide_SOAs

               type strength  n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 16       4 10   10    perhaps
    ZT2019   ZT2019       3- 24       4 11   11        yes
                                        code
    HCT2018    SOAs2plus_regular(2, 4, m=10)
    ZT2019  OSOAs_hadamard(m=11, n=24, el=2)

---

               type strength  n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 32       4 21   22    perhaps
                                     code
    HCT2018 SOAs2plus_regular(2, 5, m=21)

---

               type strength  n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 32       4 22   22    perhaps
    ZT2019   ZT2019       3- 32       4 15   15        yes
                                        code
    HCT2018    SOAs2plus_regular(2, 5, m=22)
    ZT2019  OSOAs_hadamard(m=15, n=32, el=2)

---

               type strength  n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 32       4 20   22    perhaps
    ZT2019   ZT2019       3- 48       4 20   23        yes
                                        code
    HCT2018    SOAs2plus_regular(2, 5, m=20)
    ZT2019  OSOAs_hadamard(m=20, n=48, el=2)

---

    Code
      temp5
    Output
                 type strength  n nlevels m mmax orthogonal
      HCT2018 HCT2018       2+ 64      16 8    8    perhaps
      ZT2019   ZT2019 3- or 2+ 64      16 5    5        yes
                                            code
      HCT2018   SOAs2plus_regular(s=4, k=3, m=8)
      ZT2019  OSOAs_regular(s=4, k=3, el=2, m=5)

---

    Code
      temp6
    Output
                 type strength   n nlevels  m mmax orthogonal
      HCT2018 HCT2018       2+ 256      16 45   45    perhaps
      ZT2019   ZT2019 3- or 2+ 256      16 21   21        yes
                                             code
      HCT2018   SOAs2plus_regular(s=4, k=4, m=45)
      ZT2019  OSOAs_regular(s=4, k=4, el=2, m=21)

---

               type strength   n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 256      16 15   45    perhaps
    ZT2019   ZT2019 3- or 2+ 256      16 15   21        yes
                                           code
    HCT2018   SOAs2plus_regular(s=4, k=4, m=15)
    ZT2019  OSOAs_regular(s=4, k=4, el=2, m=15)

---

               type strength  n nlevels m mmax orthogonal
    HCT2018 HCT2018       2+ 27       9 6    6    perhaps
    ZT2019   ZT2019 3- or 2+ 27       9 4    4        yes
                                          code
    HCT2018   SOAs2plus_regular(s=3, k=3, m=6)
    ZT2019  OSOAs_regular(s=3, k=3, el=2, m=4)

---

               type strength  n nlevels m mmax orthogonal
    HCT2018 HCT2018       2+ 27       9 6    6    perhaps
    ZT2019   ZT2019 3- or 2+ 27       9 4    4        yes
                                          code
    HCT2018   SOAs2plus_regular(s=3, k=3, m=6)
    ZT2019  OSOAs_regular(s=3, k=3, el=2, m=4)

---

               type strength   n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+  81       9 20   25    perhaps
    ZT2019   ZT2019 3- or 2+ 243       9 20   40        yes
                                           code
    HCT2018   SOAs2plus_regular(s=3, k=4, m=20)
    ZT2019  OSOAs_regular(s=3, k=5, el=2, m=20)

---

               type strength  n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 81       9 25   25    perhaps
    ZT2019   ZT2019 3- or 2+ 81       9 13   13        yes
                                           code
    HCT2018   SOAs2plus_regular(s=3, k=4, m=25)
    ZT2019  OSOAs_regular(s=3, k=4, el=2, m=13)

---

               type strength   n nlevels  m mmax orthogonal
    HCT2018 HCT2018       2+ 729      81 18   18    perhaps
    ZT2019   ZT2019 3- or 2+ 729      81 10   10        yes
                                           code
    HCT2018   SOAs2plus_regular(s=9, k=3, m=18)
    ZT2019  OSOAs_regular(s=9, k=3, el=2, m=10)

---

               type strength   n nlevels m mmax orthogonal
    HCT2018 HCT2018       2+ 512      64 8   16    perhaps
    ZT2019   ZT2019 3- or 2+ 512      64 8    9        yes
                                          code
    HCT2018   SOAs2plus_regular(s=8, k=3, m=8)
    ZT2019  OSOAs_regular(s=8, k=3, el=2, m=8)

---

               type strength  n nlevels  m mmax orthogonal
    ST_fam1 ST_fam1        3 16       8  5    5         no
    ST_fam2 ST_fam2        3 16       8  4    4         no
    ST_fam3 ST_fam3       3+ 16       8  3    3        yes
    LLY         LLY        3 24       8 10   10        yes
                                                          code
    ST_fam1     SOAs_8level(n=16, m=5, constr='ShiTang_alpha')
    ST_fam2 SOAs_8level(n=16, m=4, constr='ShiTang_alphabeta')
    ST_fam3 SOAs_8level(n=16, m=3, constr='ShiTang_alphabeta')
    LLY                       OSOAs_hadamard(m=10, n=24, el=3)

---

        type strength  n nlevels  m mmax orthogonal
    LLY  LLY        3 32       8 12   14        yes
                                    code
    LLY OSOAs_hadamard(m=12, n=32, el=3)

---

               type strength  n nlevels  m mmax orthogonal
    ST_fam1 ST_fam1        3 32       8  9    9         no
    ST_fam2 ST_fam2        3 32       8  8    8         no
    ST_fam3 ST_fam3       3+ 32       8  7    7        yes
    LLY         LLY        3 32       8 14   14        yes
                                                          code
    ST_fam1     SOAs_8level(n=32, m=9, constr='ShiTang_alpha')
    ST_fam2 SOAs_8level(n=32, m=8, constr='ShiTang_alphabeta')
    ST_fam3 SOAs_8level(n=32, m=7, constr='ShiTang_alphabeta')
    LLY                       OSOAs_hadamard(m=14, n=32, el=3)

---

               type strength   n nlevels  m mmax orthogonal
    ST_fam1 ST_fam1        3  64       8 20   20         no
    ST_fam3 ST_fam3       3+ 128       8 20   31        yes
    LLY         LLY        3  48       8 20   22        yes
                                                            code
    ST_fam1      SOAs_8level(n=64, m=20, constr='ShiTang_alpha')
    ST_fam3 SOAs_8level(n=128, m=20, constr='ShiTang_alphabeta')
    LLY                         OSOAs_hadamard(m=20, n=48, el=3)

---

               type strength   n nlevels  m mmax orthogonal
    ST_fam1 ST_fam1        3 128       8 32   40         no
    ST_fam2 ST_fam2        3 128       8 32   32         no
    LLY         LLY        3  72       8 32   34        yes
                                                            code
    ST_fam1     SOAs_8level(n=128, m=32, constr='ShiTang_alpha')
    ST_fam2 SOAs_8level(n=128, m=32, constr='ShiTang_alphabeta')
    LLY                         OSOAs_hadamard(m=32, n=72, el=3)

---

        type strength   n nlevels  m mmax orthogonal
    LLY  LLY  3 or 2* 256      64 20   20        yes
                                       code
    LLY OSOAs_regular(s=4, k=4, el=3, m=20)

---

        type strength   n nlevels  m mmax orthogonal
    LLY  LLY  3 or 2* 256      64 15   20        yes
                                       code
    LLY OSOAs_regular(s=4, k=4, el=3, m=15)

---

        type strength    n nlevels   m mmax orthogonal
    LLY  LLY  3 or 2* 4096      64 100  340        yes
                                        code
    LLY OSOAs_regular(s=4, k=6, el=3, m=100)

---

        type strength   n nlevels  m mmax orthogonal
    LLY  LLY  3 or 2* 625     125 10   30        yes
                                       code
    LLY OSOAs_regular(s=5, k=4, el=3, m=10)

---

             type strength n nlevels m mmax orthogonal
    ZT2019 ZT2019       3- 8       4 3    3        yes
                                     code
    ZT2019 OSOAs_hadamard(m=3, n=8, el=2)

