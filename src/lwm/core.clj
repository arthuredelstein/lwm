(ns lwm.core
  (:import (java.awt.geom Point2D$Double)
           (org.apache.commons.math.linear Array2DRowRealMatrix
                                           LUDecompositionImpl))
  (:gen-class
    :name LocalWeightedMean
    :init [generate-lwm-fn]
    :constructors {[Integer/TYPE List] []}
    :state lwm-function
    :methods [[transform [Point2D$Double] Point2D$Double]]))

(defn nearest-neighbors [point-pair n point-pairs]
  (take n (sort-by #(.distance (first point-pair) (first %)) point-pairs)))

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

(defn neighbor-distance [sorted-neighbor-group]
  (.distance (first sorted-neighbor-group) (last sorted-neighbor-group)))
        
(defrecord control-point [point Rn polyX polyY])

(defn create-control-point [point-pair order point-pairs]
  (let [exponents (polynomial-exponents order)
        neighbors (nearest-neighbors point-pair
                                     (count exponents) point-pairs)
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
                             (:Rn control-point))
                        eval-part #(evaluate-polynomial (.x pt) (.y pt)
                                                        (% control-point)
                                                        exponents)]
                    (when-let [weight (weight-function R)]
                      { :w weight
                       :wpx (* weight (eval-part :polyX))
                       :wpy (* weight (eval-part :polyY))}))))

(defn generate-lwm-fn [order point-pairs]
  (let [exponents (polynomial-exponents order)
        _ (println exponents)
        control-points (create-control-points order point-pairs)]
    (fn [pt]
      (weighted-mean-xy
        (filter identity
                (make-component-maps pt exponents control-points))))))
        

(defn -transform [this src-point]
  ((.lwm-function this) src-point))
              
;; tests

(defn create-random-points [n]
  (repeatedly n #(Point2D$Double. (rand) (rand))))

(defn test-nearest-neighbor-map []
  (nearest-neighbor-map 5 (create-random-points 1000)))