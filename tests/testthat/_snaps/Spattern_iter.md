# Spattern_iter

    Code
      Spattern_iter(nullcase, s = 8)
    Output
        1   2   3   4 
        0  42 168 301 
      attr(,"call")
      Spattern_iter(nullcase, s = 8)

---

    Code
      Spattern_iter(nullcase, s = 8, maxdim = 3, maxwt = 3)
    Output
        1   2   3 
        0  42 168 
      attr(,"call")
      Spattern_iter(nullcase, s = 8, maxdim = 3, maxwt = 3)

---

    1 2 3 4 
    0 0 0 7 
    attr(,"call")
    Spattern_iter(st1, s = 2, maxdim = 4, maxwt = 4)

---

    1 2 
    0 0 
    attr(,"call")
    Spattern_iter(st1, s = 2, maxdim = 3, maxwt = 2, detailed = TRUE)
    attr(,"details")
    attr(,"details")[[1]]
    attr(,"details")[[1]]$combis
           
    cols  1
    cols  2
    cols  3
    cols  8
    cols  9
    cols 10
    cols 15
    cols 16
    cols 17
    cols 22
    cols 23
    cols 24
    cols 29
    cols 30
    cols 31
    cols 36
    cols 37
    cols 38
    cols 43
    cols 44
    cols 45
    
    attr(,"details")[[1]]$combiweights
     [1] 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2
    
    attr(,"details")[[1]]$contribs
    [1] 8.835242e-29 1.767048e-28
    
    
    attr(,"details")[[2]]
    attr(,"details")[[2]]$combis
              
    cols  1  8
    cols  1 15
    cols  1 22
    cols  1 29
    cols  1 36
    cols  1 43
    cols  8 15
    cols  8 22
    cols  8 29
    cols  8 36
    cols  8 43
    cols 15 22
    cols 15 29
    cols 15 36
    cols 15 43
    cols 22 29
    cols 22 36
    cols 22 43
    cols 29 36
    cols 29 43
    cols 36 43
    
    attr(,"details")[[2]]$combiweights
     [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    
    attr(,"details")[[2]]$contribs
    [1] NA  0
    
    

---

     1  2  3  4  5  6  7  8 
     0  6  0 13 24 12  0  8 
    attr(,"call")
    Spattern_iter(nullcase4, s = 2, maxwt = NULL)

---

     1  2  3  4 
     0  6  0 13 
    attr(,"call")
    Spattern_iter(nullcase4, s = 2, maxdim = NULL)

---

     1  2  3  4  5  6  7  8 
     0  6  0 13 24 12  0  8 
    attr(,"call")
    Spattern_iter(nullcase4, s = 2, maxdim = NULL, maxwt = NULL)

