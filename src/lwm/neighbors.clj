(ns lwm.neighbors
    (:gen-class
    :name valelab.NearestNeighborFinder
    :init init
    :constructors {[java.util.List] nil}
    :state finderFunction
    :methods [[find [java.awt.geom.Point2D$Double int] java.util.List]]))

(defn binary-search
  ([val coll]
    (binary-search val < coll))
  ([val comparator coll]
    (let [i (java.util.Collections/binarySearch coll val comparator)]
      (if (neg? i)
        (dec (Math/abs i))
        i))))

(defn positions [point by-x by-y]
  [(binary-search point #(compare (.x %1) (.x %2)) by-x)
   (binary-search point #(compare (.y %1) (.y %2)) by-y)])        

(defn get-rect [xi0 yi0 n by-x by-y]
  {:r (get by-x (+ xi0 n))
   :l (get by-x (- xi0 n))
   :t (get by-y (+ yi0 n))
   :b (get by-y (- yi0 n))})

(defn expanse [point edge-points]
  (let [x0 (.x point)
        y0 (.y point)
        {:keys [r l t b]} edge-points]
    [(when r
        (- (.x r) x0))
     (when l
        (- x0 (.x l)))
     (when t
       (- (.y t) y0))
     (when b
       (- y0 (.y b)))]))

(defn sort-points [points]
  [(vec (sort-by #(.x %) points))
   (vec (sort-by #(.y %) points))])

(defn nearest-neighbors [point n by-x by-y]
  (let [[xi0 yi0] (positions point by-x by-y)]
    (loop [i 0 neighbor-candidates []]
      (let [edge-points (get-rect xi0 yi0 i by-x by-y)
            furthest-candidate (last neighbor-candidates)]
        (if (and (<= n (count neighbor-candidates))
                 (< (.distance point furthest-candidate)
                    (apply min (filter identity (expanse point edge-points)))))
          neighbor-candidates
          (let [all-candidates (set (concat neighbor-candidates (filter identity (vals edge-points))))
                new-candidates (take n (sort-by #(.distance point %) all-candidates))]
            (recur (inc i) new-candidates)))))))
          
(defn nearest-neighbor-finder [points]
  (let [sorted-points (sort-points points)]
    (fn [point n]
      (apply nearest-neighbors point n sorted-points))))

(defn -init [points]
  [[] (nearest-neighbor-finder points)])

(defn -find [this point numberNeighbors]
  ((.finderFunction this) point numberNeighbors))
            
        