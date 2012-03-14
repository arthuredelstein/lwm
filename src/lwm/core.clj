(ns lwm.core
  (:import (java.awt.geom Point2D$Double)
           (java.util ArrayList)
           (org.apache.commons.math.linear Array2DRowRealMatrix
                                           LUDecompositionImpl)
           (ags.utils KdTree$SqrEuclid))
  (:require lwm.neighbors)
  (:gen-class
    :name valelab.LocalWeightedMean
    :init init
    :constructors {[Integer java.util.Map] nil}
    :state lwmFunction
    :methods [[transform [java.awt.geom.Point2D$Double] java.awt.geom.Point2D$Double]
              ^{:static true} [randomTestPoints [Integer] java.util.List]
              ^{:static true} [findNeighbor  [java.awt.geom.Point2D$Double
                                             java.util.List]
                                 java.awt.geom.Point2D$Double]
              ^{:static true} [findNeighbors [java.awt.geom.Point2D$Double
                                             Integer
                                             java.util.List]
                                 java.util.List]]))

(defn nearest-neighbor-brute [point points]
  (apply min-key #(.distance point %) points))

(defn nearest-neighbors-brute [point n points]
  (take n (sort-by #(.distance point %) points)))

(defn point-to-array [^Point2D$Double point]
  (double-array [(.x point) (.y point)]))

(defn setup-kd-tree [points]
  (let [points (vec points)
        tree (KdTree$SqrEuclid. 2 (Integer. Integer/MAX_VALUE))]
    (dotimes [i (count points)]
      (let [point (points i)]
        (.addPoint tree (point-to-array point) i)))
    {:kd-tree tree
     :points points}))

(defn nearest-neighbors-kdtree [kd-tree point n sort?]
  (let [entries (.nearestNeighbor (:kd-tree kd-tree)
                                  (point-to-array point) n sort?)]
    (map #((:points kd-tree) (.value %1)) entries)))
 
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

(defn create-control-point [kd-tree src-point order point-map]
  (let [exponents (polynomial-exponents order)
        n (count exponents)
        neighbors (nearest-neighbors-kdtree kd-tree src-point n true)
        [polyX polyY] (fit-polynomial exponents
                                      (select-keys point-map neighbors))]
    (control-point. src-point
                    (.distance src-point (first neighbors))
                    polyX
                    polyY)))

(defn create-control-points [kd-tree order point-map]
  (let [src-points (keys point-map)]
    (map #(create-control-point
            kd-tree % order point-map) src-points)))

(defn sum-vals [m kw]
  (apply + (map kw m)))

(defn weighted-mean-xy [component-maps]
  (let [ws (sum-vals component-maps :w)
        wpxs (sum-vals component-maps :wpx)
        wpys (sum-vals component-maps :wpy)]
    (Point2D$Double. (/ wpxs ws) (/ wpys ws))))

(defn make-component-maps [pt exponents control-points]
  (for [control-point control-points]
    (do
    (let [R (/ (.distance pt (:point control-point))
               (:Rn control-point))]
      (when-let [weight (weight-function R)]
        (let [eval-part #(evaluate-polynomial (.x pt) (.y pt)
                                              (% control-point)
                                              exponents)]
          {:w weight
           :wpx (* weight (eval-part :polyX))
           :wpy (* weight (eval-part :polyY))}))))))

(defn generate-lwm-fn [order point-map]
  (def q point-map)
  (let [kd-tree (setup-kd-tree (keys point-map))
        exponents (polynomial-exponents order)
        control-points (create-control-points kd-tree order point-map)]
    (fn [pt]
      (weighted-mean-xy
        (filter identity
                (make-component-maps pt exponents control-points))))))
        
;; java interop

(defn -init [order point-map]
  [[] (generate-lwm-fn order point-map)])

(defn -transform [this src-point]
  ((.lwmFunction this) src-point))
              
(defn -findNeighbors [point n points]
  (ArrayList. (nearest-neighbors-brute point n points)))

(defn -findNeighbor [point points]
  (nearest-neighbor-brute point points))
  
  
  ;; tests

(defn create-random-points [n]
  (repeatedly n #(Point2D$Double. (rand) (rand))))

(defn -randomTestPointPairs [n]
  (into {} (map vec (partition 2 (create-random-points (* 2 n))))))


(defn find-neighbors [size n k]
  (let [q (create-random-points size)]
    (time (doseq [q0 (take k q)]
            (lwm.neighbors/nearest-neighbors q0 n q)))))

(defn find-neighbors-other [size n k]
  (let [q (create-random-points size)]
    (time
      (let [finder (lwm.neighbors/nearest-neighbor-finder q)]
        (doseq [q0 (take k q)]
          (finder q0 n))))))
