(TeX-add-style-hook "appendix"
 (function
  (lambda ()
    (LaTeX-add-labels
     "sec:authors")
    (TeX-add-symbols
     "thepage"
     "centerformat"
     "sixrm"
     "egtrm"
     "title"
     "authors"
     "inst"
     "fn"))))

