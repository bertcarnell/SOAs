# ocheck

    Table of correlation matrix entries:
    
     1 
    16 
    [1] FALSE

# ocheck3

    Triples that violate 3-orthogonality:
    [1] TRUE

---

    Triples that violate 3-orthogonality:
    [1] "4,2,1"
    [1] "4,3,1"
    [1] "2,4,1"
    [1] "3,4,1"
    [1] "4,1,2"
    [1] "4,3,2"
    [1] "1,4,2"
    [1] "3,4,2"
    [1] "4,1,3"
    [1] "4,2,3"
    [1] "1,4,3"
    [1] "2,4,3"
    [1] "2,1,4"
    [1] "3,1,4"
    [1] "1,2,4"
    [1] "3,2,4"
    [1] "1,3,4"
    [1] "2,3,4"
    [1] FALSE

# soacheck2D

    pairs for which SOA property in 2D is violated:
    [1] 1 2
    [1] "1x1: A2 = 1"
    [1] 1 3
    [1] "1x1: A2 = 1"
    [1] 1 4
    [1] "1x1: A2 = 1"
    [1] 2 3
    [1] "1x1: A2 = 1"
    [1] 2 4
    [1] "1x1: A2 = 1"
    [1] 3 4
    [1] "1x1: A2 = 1"
    [1] FALSE

---

    pairs for which SOA property in 2D is violated:

---

    pairs for which SOA property in 2D is violated:
    [1] 1 2
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"
    [1] 1 3
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"
    [1] 1 4
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"
    [1] 2 3
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"
    [1] 2 4
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"
    [1] 3 4
    [1] "3x1: A2 = 1"
    [1] "2x2: A2 = 3"
    [1] "1x3: A2 = 1"

# soacheck3D

    triples for which SOA property in 3D is violated:
    [1] 1 2 3
    2x1x1:
        1     2     3 
    0.000 0.000 0.031 

# Spattern

    Code
      Spattern(nullcase, s = 8)
    Message <simpleMessage>
      s equals the number of column levels. Function GWLP from package DoE.base is used.
    Output
        1   2   3   4 
        0  42 168 301 
      attr(,"call")
      Spattern(nullcase, s = 8)
      attr(,"class")
      [1] "Spaper"

---

    Code
      Spattern(nullcase, s = 8, maxdim = 3, maxwt = 3)
    Message <simpleMessage>
      s equals the number of column levels. Function GWLP from package DoE.base is used.
    Output
        1   2   3 
        0  42 168 
      attr(,"call")
      Spattern(nullcase, s = 8, maxdim = 3, maxwt = 3)
      attr(,"class")
      [1] "Spaper"

---

    1 2 3 4 
    0 0 0 7 
    attr(,"call")
    Spattern(st1, s = 2, maxdim = 4, maxwt = 4)

---

    1 2 
    0 0 
    attr(,"call")
    Spattern(st1, s = 2, maxdim = 3, maxwt = 2, detailed = TRUE)

---

     1  2  3  4  5  6  7  8 
     0  6  0 13 24 12  0  8 
    attr(,"call")
    Spattern(nullcase4, s = 2, maxwt = NULL)

---

     1  2  3  4 
     0  6  0 13 
    attr(,"call")
    Spattern(nullcase4, s = 2, maxdim = NULL)

---

     1  2  3  4  5  6  7  8 
     0  6  0 13 24 12  0  8 
    attr(,"call")
    Spattern(nullcase4, s = 2, maxdim = NULL, maxwt = NULL)

