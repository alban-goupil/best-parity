#! /usr/bin/env clj 

(let [M (Integer/parseUnsignedInt (first *command-line-args*))
      N (int (Math/sqrt M))
      A (vec (range N))]
;      A (into [] (range (- 1 N) N 2))]
  (doseq [q (range N)
          i (range N)]
    (printf "%2d\t%2d\n" (A i) (A q))))
