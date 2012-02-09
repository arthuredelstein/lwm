(ns lwm.core
  (:import (java.awt.geom Point2D$Double)
           (org.apache.commons.math.linear Array2DRowRealMatrix
                                           LUDecompositionImpl))
  (:gen-class
    :name valelab.LocalWeightedMean
    :init init
    :constructors {[Integer java.util.List] nil}
    :state lwmFunction
    :methods [[transform [java.awt.geom.Point2D$Double] java.awt.geom.Point2D$Double]
              ^{:static true} [randomTestPoints [Integer] java.util.List]]))

(defn nearest-neighbors [point n points]
  (take n (sort-by #(.distance point %) points)))

(defn weight-function [R]
  (when (< R 1)
    (+ 1 (* -3 R R) (* 2 R R R))))

(defn polynomial-exponents [n]
  (apply concat
    (for [j (range (inc n))]
      (for [k (reverse (range (inc j)))]
        [k (- j k)]))))

(defn power-terms [x y exponents]
  (for [[x-exp y-exp] exponents]
    (* (Math/pow x x-exp) (Math/pow y y-exp))))

(defn make-matrix [nested-lists]
  (Array2DRowRealMatrix.
    (into-array (map double-array nested-lists))))
  
(defn fit-polynomial [exponents point-pairs]
  (let [src-points (map first point-pairs)
        dest-points (map second point-pairs)
        matrix (make-matrix
                 (for [src-point src-points]
                   (power-terms (.x src-point) (.y src-point) exponents)))
        solver (.getSolver (LUDecompositionImpl. matrix))
        destX (double-array (map #(.x %) dest-points))
        destY (double-array (map #(.y %) dest-points))]
    [(.solve solver destX) (.solve solver destY)]))

(defn evaluate-polynomial [x y coeffs exponents]
  (apply + (map * coeffs (power-terms x y exponents))))
        
(defrecord control-point [point Rn polyX polyY])

(defn create-control-point [point-pair order point-pairs]
  (let [exponents (polynomial-exponents order)
        neighbors (nearest-neighbors (first point-pair)
                                     (count exponents)
                                     (map first point-pairs))
        src-point (first point-pair)
        [polyX polyY] (fit-polynomial exponents neighbors)]
    (control-point. src-point
                    (.distance src-point (first (last neighbors)))
                    polyX
                    polyY)))

(defn create-control-points [order point-pairs]
  (map #(create-control-point % order point-pairs) point-pairs))

(defn sum-vals [m kw]
  (apply + (map kw m)))

(defn weighted-mean-xy [component-maps]
  (let [ws (sum-vals component-maps :w)
        wpxs (sum-vals component-maps :wpx)
        wpys (sum-vals component-maps :wpy)]
    (Point2D$Double. (/ wpxs ws) (/ wpys ws))))

(defn make-component-maps [pt exponents control-points]
  (for [control-point control-points]
    (let [R (/ (.distance pt (:point control-point))
               (:Rn control-point))]
      (when-let [weight (weight-function R)]
        (let [eval-part #(evaluate-polynomial (.x pt) (.y pt)
                                              (% control-point)
                                              exponents)]
          {:w weight
           :wpx (* weight (eval-part :polyX))
           :wpy (* weight (eval-part :polyY))})))))

(defn generate-lwm-fn [order point-pairs]
  (let [exponents (polynomial-exponents order)
        control-points (create-control-points order point-pairs)]
    (fn [pt]
      (weighted-mean-xy
        (filter identity
                (make-component-maps pt exponents control-points))))))
        
;; java interop

(defn -init [order point-pairs]
  [[] (generate-lwm-fn order point-pairs)])

(defn -transform [this src-point]
  ((.lwmFunction this) src-point))
              
;; tests

(defn create-random-points [n]
  (repeatedly n #(Point2D$Double. (rand) (rand))))

(defn -randomTestPoints [n]
  (partition 2 (create-random-points (* 2 n))))

(defn find-neighbors [size n k]
  (let [q (create-random-points size)]
    (time (doseq [q0 (take k q)]
            (nearest-neighbors q0 n q)))))

(defn find-neighbors-other [size n k]
  (let [q (create-random-points size)]
    (time
      (let [finder (lwm.neighbors/nearest-neighbor-finder q)]
        (doseq [q0 (take k q)]
          (finder q0 n))))))
     
