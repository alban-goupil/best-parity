#! /usr/bin/env clj

;;; Génère le mapping du DVB-T2 en 16 QAM

(def mappings
  {4   [[1] [0]]
   16  [[1 0] [1 1] [0 1] [0 0]]
   64  [[1 0 0] [1 0 1] [1 1 1] [1 1 0]
        [0 1 0] [0 1 1] [0 0 1] [0 0 0]]
   256 [[1 0 0 0] [1 0 0 1] [1 0 1 1] [1 0 1 0]
        [1 1 1 0] [1 1 1 1] [1 1 0 1] [1 1 0 0]
        [0 1 0 0] [0 1 0 1] [0 1 1 1] [0 1 1 0]
        [0 0 1 0] [0 0 1 1] [0 0 0 1] [0 0 0 0]]})


(let [order (Integer/parseUnsignedInt (first *command-line-args*))
      mapping (mappings order)
      bits2int (fn [bs] (reduce + (map * (iterate #(* 2 %) 1) bs)))]

  (->> (for [q mapping i mapping] (interleave i q))
       (map-indexed #(hash-map :const %1 :gf (bits2int %2)))
       (sort-by :gf)
       (map :const)
       (run! #(print % ""))))
