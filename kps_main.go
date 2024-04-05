package main

import (
	"fmt"
	"math"
	"sort"
)

type Point struct {
	X float64
	Y float64
}

type Pair struct {
	First  Point
	Second Point
}

type KPS struct {
	Points []Point
}

func (kps *KPS) findMedian(a []int, n int) float64 {
	if n%2 == 0 {
		sort.Ints(a)
		return float64(a[n/2-1]+a[n/2]) / 2.0
	} else {
		sort.Ints(a)
		return float64(a[n/2])
	}
}

func (kps *KPS) swap(a *float64, b *float64) {
	temp := *a
	*a = *b
	*b = temp
}

func (kps *KPS) partition(arr []float64, l int, r int, x float64) int {
	i := l
	for i < r {
		if arr[i] == x {
			break
		}
		i++
	}
	kps.swap(&arr[i], &arr[r])

	i = l
	for j := l; j <= r-1; j++ {
		if arr[j] <= x {
			kps.swap(&arr[i], &arr[j])
			i++
		}
	}
	kps.swap(&arr[i], &arr[r])
	return i
}

func (kps *KPS) kthSmallest(arr []float64, l int, r int, k int) float64 {
	if k > 0 && k <= r-l+1 {
		n := r - l + 1

		var median []float64
		for i := 0; i < n/5; i++ {
			median = append(median, kps.findMedian(arr[l+i*5:l+i*5+5], 5))
		}
		if n%5 != 0 {
			median = append(median, kps.findMedian(arr[l+n/5*5:l+n], n%5))
		}

		var medOfMed float64
		if len(median) == 1 {
			medOfMed = median[0]
		} else {
			medOfMed = kps.kthSmallest(median, 0, len(median)-1, (len(median)+1)/2)
		}

		pos := kps.partition(arr, l, r, medOfMed)

		if pos-l == k-1 {
			return arr[pos]
		}
		if pos-l > k-1 {
			return kps.kthSmallest(arr, l, pos-1, k)
		}

		return kps.kthSmallest(arr, pos+1, r, k-pos+l-1)
	}

	return math.Inf(1)
}

func (kps *KPS) findMedianDouble(arr []float64, n int) float64 {
	sort.Float64s(arr)
	return arr[n/2]
}

func (kps *KPS) abss(a float64) float64 {
	if a < 0 {
		return 0 - a
	} else {
		return a
	}
}

func (kps *KPS) getT(p1 Point, p2 Point, points []Point, flag bool) []Point {
	var upperT []Point
	slope := (p1.Y - p2.Y) / (p1.X - p2.X)
	for i := 0; i < len(points); i++ {
		currPoint := points[i]

		if currPoint.X > p1.X && currPoint.X < p2.X {

			currSlope := (p1.Y - currPoint.Y) / (p1.X - currPoint.X)
			if !flag {
				if currSlope > slope {
					upperT = append(upperT, currPoint)
				}
			} else {
				if currSlope < slope {
					upperT = append(upperT, currPoint)
				}
			}
		}
	}
	upperT = append(upperT, p1)
	upperT = append(upperT, p2)

	return upperT
}

func (kps *KPS) getUpperBridge(points []Point, median float64) Pair {
	sort.Slice(points, func(i, j int) bool {
		return points[i].X < points[j].X
	})

	var candidates []Point
	var pairs []Pair
	if len(points)%2 == 0 {
		for i := 0; i < len(points); i += 2 {
			firstPt := points[i]
			secondPt := points[i+1]

			currPair := Pair{First: firstPt, Second: secondPt}
			pairs = append(pairs, currPair)
		}
	} else {
		candidates = append(candidates, points[0])
		for i := 1; i < len(points); i += 2 {
			firstPt := points[i]
			secondPt := points[i+1]

			currPair := Pair{First: firstPt, Second: secondPt}
			pairs = append(pairs, currPair)
		}
	}
	var slopes []float64
	for i := 0; i < len(pairs); i++ {

		p1 := pairs[i].First
		p2 := pairs[i].Second
		x1 := p1.X
		x2 := p2.X
		y1 := p1.Y
		y2 := p2.Y

		if x1 == x2 {
			if y1 > y2 {
				candidates = append(candidates, p1)
			} else {
				candidates = append(candidates, p2)
			}
			slopes = append(slopes, math.Inf(1))
		} else {
			slope := (y2 - y1) / (x2 - x1)
			slopes = append(slopes, slope)
		}
	}

	var arr []float64
	for i := 0; i < len(slopes); i++ {
		if slopes[i] != math.Inf(1) {
			arr = append(arr, slopes[i])
		}
	}

	var medianSlope float64
	if len(arr) == 1 {
		medianSlope = arr[0]
	} else {
		medianSlope = kps.kthSmallest(arr, 0, len(arr)-1, (len(arr)+1)/2)
	}

	var SMALL []Pair
	var EQUAL []Pair
	var LARGE []Pair

	for i := 0; i < len(pairs); i++ {
		p1 := pairs[i].First
		p2 := pairs[i].Second
		x1 := p1.X
		x2 := p2.X
		y1 := p1.Y
		y2 := p2.Y

		if x1 != x2 {
			slope := (y2 - y1) / (x2 - x1)
			if math.Abs(slope-medianSlope) < 0.001 {
				currPair := Pair{First: p1, Second: p2}
				EQUAL = append(EQUAL, currPair)
			} else if slope < medianSlope {
				currPair := Pair{First: p1, Second: p2}
				SMALL = append(SMALL, currPair)
			} else if slope > medianSlope {
				currPair := Pair{First: p1, Second: p2}
				LARGE = append(LARGE, currPair)
			}
		}
	}

	maxC := math.Inf(-1)
	for i := 0; i < len(points); i++ {

		x := points[i].X
		y := points[i].Y
		currC := y - medianSlope*x

		if currC > maxC {
			maxC = currC
		}
	}

	pmin := Point{X: math.Inf(1), Y: math.Inf(1)}
	pmax := Point{X: math.Inf(-1), Y: math.Inf(-1)}

	for i := 0; i < len(points); i++ {

		x := points[i].X
		y := points[i].Y

		currC := y - medianSlope*x

		if math.Abs(currC-maxC) < 0.001 {

			if x < pmin.X {
				pmin.X = x
				pmin.Y = y
			}
			if x > pmax.X {
				pmax.X = x
				pmax.Y = y
			}
		}
	}

	if pmin.X <= median && pmax.X > median {
		upperBridge := Pair{First: pmin, Second: pmax}
		return upperBridge
	} else if pmax.X <= median {
		for i := 0; i < len(EQUAL); i++ {
			pt := EQUAL[i].Second
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(LARGE); i++ {
			pt := LARGE[i].Second
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(SMALL); i++ {
			pt1 := SMALL[i].First
			pt2 := SMALL[i].Second
			candidates = append(candidates, pt1)
			candidates = append(candidates, pt2)
		}
		return kps.getUpperBridge(candidates, median)
	} else if pmin.X > median {
		for i := 0; i < len(EQUAL); i++ {
			pt := EQUAL[i].First
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(LARGE); i++ {
			pt1 := LARGE[i].First
			pt2 := LARGE[i].Second
			candidates = append(candidates, pt1)
			candidates = append(candidates, pt2)
		}
		for i := 0; i < len(SMALL); i++ {
			pt := SMALL[i].First
			candidates = append(candidates, pt)
		}
		return kps.getUpperBridge(candidates, median)
	}
	return Pair{}
}

func (kps *KPS) getLowerBridge(points []Point, median float64) Pair {
	sort.Slice(points, func(i, j int) bool {
		return points[i].X < points[j].X
	})

	var candidates []Point
	var pairs []Pair
	if len(points)%2 == 0 {
		for i := 0; i < len(points); i += 2 {
			firstPt := points[i]
			secondPt := points[i+1]

			currPair := Pair{First: firstPt, Second: secondPt}
			pairs = append(pairs, currPair)
		}
	} else {
		candidates = append(candidates, points[0])
		for i := 1; i < len(points); i += 2 {
			firstPt := points[i]
			secondPt := points[i+1]

			currPair := Pair{First: firstPt, Second: secondPt}
			pairs = append(pairs, currPair)
		}
	}

	var slopes []float64
	for i := 0; i < len(pairs); i++ {

		p1 := pairs[i].First
		p2 := pairs[i].Second
		x1 := p1.X
		x2 := p2.X
		y1 := p1.Y
		y2 := p2.Y

		if x1 == x2 {
			if y1 > y2 {
				candidates = append(candidates, p2)
			} else {
				candidates = append(candidates, p1)
			}
			slopes = append(slopes, math.Inf(1))
		} else {
			slope := (y2 - y1) / (x2 - x1)
			slopes = append(slopes, slope)
		}
	}

	var arr []float64
	for i := 0; i < len(slopes); i++ {
		if slopes[i] != math.Inf(1) {
			arr = append(arr, slopes[i])
		}
	}

	var medianSlope float64
	if len(arr) == 1 {
		medianSlope = arr[0]
	} else {
		medianSlope = kps.kthSmallest(arr, 0, len(arr)-1, (len(arr)+1)/2)
	}

	var SMALL []Pair
	var EQUAL []Pair
	var LARGE []Pair

	for i := 0; i < len(pairs); i++ {
		p1 := pairs[i].First
		p2 := pairs[i].Second
		x1 := p1.X
		x2 := p2.X
		y1 := p1.Y
		y2 := p2.Y

		if x1 != x2 {
			slope := (y2 - y1) / (x2 - x1)
			if math.Abs(slope-medianSlope) < 0.001 {
				currPair := Pair{First: p1, Second: p2}
				EQUAL = append(EQUAL, currPair)
			} else if slope < medianSlope {
				currPair := Pair{First: p1, Second: p2}
				SMALL = append(SMALL, currPair)
			} else if slope > medianSlope {
				currPair := Pair{First: p1, Second: p2}
				LARGE = append(LARGE, currPair)
			}
		}
	}

	minC := math.Inf(1)

	for i := 0; i < len(points); i++ {

		x := points[i].X
		y := points[i].Y
		currC := y - medianSlope*x

		if currC < minC {
			minC = currC
		}
	}

	pmin := Point{X: math.Inf(1), Y: math.Inf(1)}
	pmax := Point{X: math.Inf(-1), Y: math.Inf(-1)}

	for i := 0; i < len(points); i++ {

		x := points[i].X
		y := points[i].Y
		currC := y - medianSlope*x

		if math.Abs(currC-minC) < 0.001 {

			if x < pmin.X {
				pmin.X = x
				pmin.Y = y
			}
			if x > pmax.X {
				pmax.X = x
				pmax.Y = y
			}
		}
	}

	if pmin.X <= median && pmax.X > median {
		lowerBridge := Pair{First: pmin, Second: pmax}
		return lowerBridge
	} else if pmax.X <= median {
		for i := 0; i < len(EQUAL); i++ {
			pt := EQUAL[i].Second
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(LARGE); i++ {
			pt1 := LARGE[i].First
			pt2 := LARGE[i].Second
			candidates = append(candidates, pt1)
			candidates = append(candidates, pt2)
		}
		for i := 0; i < len(SMALL); i++ {
			pt := SMALL[i].Second
			candidates = append(candidates, pt)
		}
		return kps.getLowerBridge(candidates, median)
	} else if pmin.X > median {
		for i := 0; i < len(EQUAL); i++ {
			pt := EQUAL[i].First
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(LARGE); i++ {
			pt := LARGE[i].First
			candidates = append(candidates, pt)
		}
		for i := 0; i < len(SMALL); i++ {
			pt1 := SMALL[i].First
			pt2 := SMALL[i].Second
			candidates = append(candidates, pt1)
			candidates = append(candidates, pt2)
		}
		return kps.getLowerBridge(candidates, median)
	}
	return Pair{}
}

func (kps *KPS) getUpperHull(pmin Point, pmax Point, points []Point) []Point {

	var upperHull []Point
	n := len(points)
	var arr []float64
	for i := 0; i < n; i++ {
		arr = append(arr, points[i].X)
	}

	var median float64
	if n == 1 {
		median = arr[0]
	} else {
		median = kps.kthSmallest(arr, 0, n-1, (n+1)/2)
	}
	upperBridge := kps.getUpperBridge(points, median)

	pl := upperBridge.First
	pr := upperBridge.Second

	if pl.X > pr.X {
		temp := pl
		pl = pr
		pr = temp
	}

	upperHull = append(upperHull, pl)
	upperHull = append(upperHull, pr)

	if pmin != pl {
		upperTLeft := kps.getT(pmin, pl, points, false)
		left := kps.getUpperHull(pmin, pl, upperTLeft)
		upperHull = append(upperHull, left...)
	}

	if pmax != pr {
		upperTRight := kps.getT(pr, pmax, points, false)
		right := kps.getUpperHull(pr, pmax, upperTRight)
		upperHull = append(upperHull, right...)
	}

	return upperHull
}

func (kps *KPS) getLowerHull(pmin Point, pmax Point, points []Point) []Point {

	var lowerHull []Point
	n := len(points)
	var arr []float64
	for i := 0; i < n; i++ {
		arr = append(arr, points[i].X)
	}
	var median float64
	if n == 1 {
		median = arr[0]
	} else {
		median = kps.kthSmallest(arr, 0, n-1, (n+1)/2)
	}
	lowerBridge := kps.getLowerBridge(points, median)

	pl := lowerBridge.First
	pr := lowerBridge.Second

	if pl.X > pr.X {
		temp := pl
		pl = pr
		pr = temp
	}

	lowerHull = append(lowerHull, pl)
	lowerHull = append(lowerHull, pr)

	if pmin != pl {
		lowerTLeft := kps.getT(pmin, pl, points, true)
		left := kps.getLowerHull(pmin, pl, lowerTLeft)
		lowerHull = append(lowerHull, left...)
	}
	if pmax != pr {
		lowerTRight := kps.getT(pr, pmax, points, true)
		right := kps.getLowerHull(pr, pmax, lowerTRight)
		lowerHull = append(lowerHull, right...)
	}

	return lowerHull
}

func (kps *KPS) ComputeHull() []Point {

	if len(kps.Points) < 3 {
		fmt.Println("Hull doesn't exist!!")
		return nil
	}

	pminU, pminL, pmaxU, pmaxL := kps.Points[0], kps.Points[0], kps.Points[0], kps.Points[0]
	for i := 1; i < len(kps.Points); i++ {
		currPoint := kps.Points[i]
		if currPoint.X < pminL.X {
			pminL = currPoint
			pminU = currPoint
		} else if currPoint.X > pmaxL.X {
			pmaxL = currPoint
			pmaxU = currPoint
		} else if currPoint.X == pminL.X {
			if currPoint.Y > pminU.Y {
				pminU = currPoint
			} else if currPoint.Y < pminL.Y {
				pminL = currPoint
			}
		} else if currPoint.X == pmaxL.X {
			if currPoint.Y > pmaxU.Y {
				pmaxU = currPoint
			} else if currPoint.Y < pmaxL.Y {
				pmaxL = currPoint
			}
		}
	}

	upperT := kps.getT(pminU, pmaxU, kps.Points, false)
	upperHull := kps.getUpperHull(pminU, pmaxU, upperT)

	lowerT := kps.getT(pminL, pmaxL, kps.Points, true)
	lowerHull := kps.getLowerHull(pminL, pmaxL, lowerT)

	hullEdges := append(upperHull, lowerHull...)

	if pminU != pminL {
		hullEdges = append(hullEdges, pminL)
		hullEdges = append(hullEdges, pminU)
	}
	if pmaxL != pmaxU {
		hullEdges = append(hullEdges, pmaxL)
		hullEdges = append(hullEdges, pmaxU)
	}

	sort.Slice(hullEdges, func(i, j int) bool {
		if hullEdges[i].X < hullEdges[j].X {
			return true
		} else if hullEdges[i].X > hullEdges[j].X {
			return false
		} else {
			return hullEdges[i].Y < hullEdges[j].Y
		}
	})

	hull := []Point{hullEdges[0]}
	i := 1
	for i < len(hullEdges) {
		for i < len(hullEdges) && hullEdges[i] == hullEdges[i-1] {
			i++
		}

		if i < len(hullEdges) {
			hull = append(hull, hullEdges[i])
		}

		i++
	}

	return hull
}

func main() {
	kps := KPS{}
	kps.Points = []Point{{X: 0, Y: 0}, {X: 1, Y: 1}, {X: 2, Y: 2}, {X: 3, Y: 3}, {X: 4, Y: 4}, {X: 5, Y: 5}}
	hull := kps.ComputeHull()
	fmt.Println(hull)
}


