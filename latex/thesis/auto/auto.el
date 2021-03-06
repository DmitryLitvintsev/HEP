(TeX-add-style-hook "auto"
 (function
  (lambda ()
    (LaTeX-add-labels
     "eq:su3"
     "eq:su6"
     "eq:wf"
     "fig:mult"
     "tab:masses"
     "eq:hadwidth"
     "eq:elwidth"
     "eq:pet"
     "fig:argus"
     "eq:decaychain"
     "fig:lambdac"
     "fig:lamcst"
     "tab:fit"
     "fig:fakes"
     "eq:cut"
     "fig:mcsic"
     "fig:sicpl"
     "fig:sic0"
     "fig:nonres"
     "eq:fraction"
     "fig:peterson1"
     "eq:dm1"
     "eq:dm2"
     "eq:dm12")
    (TeX-run-style-hooks
     "rotating"
     "epsfig"
     "lh"
     "russian"
     "latex2e"
     "art12"
     "article"
     "12pt"
     "defs"
     "biba"))))

