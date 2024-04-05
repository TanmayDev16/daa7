import math

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __lt__(self, other):
        return (self.x < other.x) or (self.x == other.x and self.y < other.y)

class KPS:
    def __init__(self):
        self.points = []

    def find_median(self, arr):
        arr.sort()
        return arr[len(arr) // 2]

    def swap(self, arr, a, b):
        arr[a], arr[b] = arr[b], arr[a]

    def partition(self, arr, l, r, x):
        i = l
        for j in range(l, r):
            if arr[j] == x:
                break
            self.swap(arr, j, r)

        i = l
        for j in range(l, r):
            if arr[j] <= x:
                self.swap(arr, i, j)
                i += 1
        self.swap(arr, i, r)
        return i

    def kth_smallest(self, arr, l, r, k):
        if k > 0 and k <= r - l + 1:
            n = r - l + 1
            i = 0
            median = []
            while i < n // 5:
                median.append(self.find_median(arr[l + i * 5: l + i * 5 + 5]))
                i += 1
            if i * 5 < n:
                median.append(self.find_median(arr[l + i * 5: l + i * 5 + n % 5]))
                i += 1

            med_of_med = median[i - 1] if i == 1 else self.kth_smallest(median, 0, i - 1, i // 2)

            pos = self.partition(arr, l, r, med_of_med)

            if pos - l == k - 1:
                return arr[pos]
            if pos - l > k - 1:
                return self.kth_smallest(arr, l, pos - 1, k)
            return self.kth_smallest(arr, pos + 1, r, k - pos + l - 1)

        return math.inf

    def abss(self, a):
        return abs(a)

    def get_T(self, p1, p2, points, flag):
        upper_T = []
        slope = (p1.y - p2.y) / (p1.x - p2.x)
        for curr_point in points:
            if p1.x < curr_point.x < p2.x:
                curr_slope = (p1.y - curr_point.y) / (p1.x - curr_point.x)
                if (not flag and curr_slope > slope) or (flag and curr_slope < slope):
                    upper_T.append(curr_point)
        upper_T.append(p1)
        upper_T.append(p2)
        return upper_T

    def get_upper_bridge(self, points, median):
        points.sort(key=lambda p: p.x)

        candidates = []
        pairs = []
        slopes = []

        if len(points) % 2 == 0:
            for i in range(0, len(points), 2):
                curr_pair = (points[i], points[i + 1])
                pairs.append(curr_pair)
        else:
            candidates.append(points[0])
            for i in range(1, len(points), 2):
                curr_pair = (points[i], points[i + 1])
                pairs.append(curr_pair)

        for p1, p2 in pairs:
            x1, x2 = p1.x, p2.x
            y1, y2 = p1.y, p2.y

            if x1 == x2:
                if y1 > y2:
                    candidates.append(p1)
                else:
                    candidates.append(p2)
                slopes.append(math.inf)
            else:
                slope = (y2 - y1) / (x2 - x1)
                slopes.append(slope)

        slopes = [s for s in slopes if s != math.inf]
        median_slope = self.kth_smallest(slopes, 0, len(slopes) - 1, (len(slopes) + 1) // 2)

        SMALL, EQUAL, LARGE = [], [], []

        for i, (p1, p2) in enumerate(pairs):
            x1, x2 = p1.x, p2.x
            y1, y2 = p1.y, p2.y

            if x1 != x2:
                slope = (y2 - y1) / (x2 - x1)
                if abs(slope - median_slope) < 0.001:
                    EQUAL.append((p1, p2))
                elif slope < median_slope:
                    SMALL.append((p1, p2))
                elif slope > median_slope:
                    LARGE.append((p1, p2))

        max_c = -math.inf
        for x, y in [(p.x, p.y) for p in points]:
            curr_c = (y - median_slope * x)
            if curr_c > max_c:
                max_c = curr_c

        pmin, pmax = Point(math.inf, math.inf), Point(-math.inf, -math.inf)

        for x, y in [(p.x, p.y) for p in points]:
            curr_c = y - median_slope * x
            if abs(curr_c - max_c) < 0.001:
                if x < pmin.x:
                    pmin.x, pmin.y = x, y
                if x > pmax.x:
                    pmax.x, pmax.y = x, y

        if pmin.x <= median < pmax.x:
            upper_bridge = (pmin, pmax)
            return upper_bridge
        elif pmax.x <= median:
            for pt in [p2 for p1, p2 in EQUAL] + [p2 for p1, p2 in LARGE] + [p1 for p1, p2 in SMALL]:
                candidates.append(pt)
            return self.get_upper_bridge(candidates, median)
        elif pmin.x > median:
            for pt in [p1 for p1, p2 in EQUAL] + [p1 for p1, p2 in LARGE] + [p1 for p1, p2 in SMALL]:
                candidates.append(pt)
            return self.get_upper_bridge(candidates, median)

    def get_lower_bridge(self, points, median):
        points.sort(key=lambda p: p.x)

        candidates = []
        pairs = []
        slopes = []

        if len(points) % 2 == 0:
            for i in range(0, len(points), 2):
                curr_pair = (points[i], points[i + 1])
                pairs.append(curr_pair)
        else:
            candidates.append(points[0])
            for i in range(1, len(points), 2):
                curr_pair = (points[i], points[i + 1])
                pairs.append(curr_pair)

        for p1, p2 in pairs:
            x1, x2 = p1.x, p2.x
            y1, y2 = p1.y, p2.y

            if x1 == x2:
                if y1 > y2:
                    candidates.append(p2)
                else:
                    candidates.append(p1)
                slopes.append(math.inf)
            else:
                slope = (y2 - y1) / (x2 - x1)
                slopes.append(slope)

        slopes = [s for s in slopes if s != math.inf]
        median_slope = self.kth_smallest(slopes, 0, len(slopes) - 1, (len(slopes) + 1) // 2)

        SMALL, EQUAL, LARGE = [], [], []

        for i, (p1, p2) in enumerate(pairs):
            x1, x2 = p1.x, p2.x
            y1, y2 = p1.y, p2.y

            if x1 != x2:
                slope = (y2 - y1) / (x2 - x1)
                if abs(slope - median_slope) < 0.001:
                    EQUAL.append((p1, p2))
                elif slope < median_slope:
                    SMALL.append((p1, p2))
                elif slope > median_slope:
                    LARGE.append((p1, p2))

        min_c = math.inf

        for x, y in [(p.x, p.y) for p in points]:
            curr_c = (y - median_slope * x)
            if curr_c < min_c:
                min_c = curr_c

        pmin, pmax = Point(math.inf, math.inf), Point(-math.inf, -math.inf)

        for x, y in [(p.x, p.y) for p in points]:
            curr_c = y - median_slope * x
            if abs(curr_c - min_c) < 0.001:
                if x < pmin.x:
                    pmin.x, pmin.y = x, y
                if x > pmax.x:
                    pmax.x, pmax.y = x, y

        if pmin.x <= median < pmax.x:
            lower_bridge = (pmin, pmax)
            return lower_bridge
        elif pmax.x <= median:
            for pt in [p2 for p1, p2 in EQUAL] + [p2 for p1, p2 in LARGE] + [p1 for p1, p2 in SMALL]:
                candidates.append(pt)
            return self.get_lower_bridge(candidates, median)
        elif pmin.x > median:
            for pt in [p1 for p1, p2 in EQUAL] + [p1 for p1, p2 in LARGE] + [p1 for p1, p2 in SMALL]:
                candidates.append(pt)
            return self.get_lower_bridge(candidates, median)

    def get_upper_hull(self, pmin, pmax, points):
        upper_hull = []
        n = len(points)
        arr = [p.x for p in points]
        median = arr[0] if n == 1 else self.kth_smallest(arr, 0, n - 1, (n + 1) // 2)
        upper_bridge = self.get_upper_bridge(points, median)

        pl, pr = upper_bridge

        if pl.x > pr.x:
            pl, pr = pr, pl

        upper_hull.append(pl)
        upper_hull.append(pr)

        if pmin != pl:
            upper_T_left = self.get_T(pmin, pl, points, False)
            left = self.get_upper_hull(pmin, pl, upper_T_left)
            upper_hull.extend(left)

        if pmax != pr:
            upper_T_right = self.get_T(pr, pmax, points, False)
            right = self.get_upper_hull(pr, pmax, upper_T_right)
            upper_hull.extend(right)

        return upper_hull

    def get_lower_hull(self, pmin, pmax, points):
        lower_hull = []
        n = len(points)
        arr = [p.x for p in points]
        median = arr[0] if n == 1 else self.kth_smallest(arr, 0, n - 1, (n + 1) // 2)
        lower_bridge = self.get_lower_bridge(points, median)

        pl, pr = lower_bridge

        if pl.x > pr.x:
            pl, pr = pr, pl

        lower_hull.append(pl)
        lower_hull.append(pr)

        if pmin != pl:
            lower_T_left = self.get_T(pmin, pl, points, True)
            left = self.get_lower_hull(pmin, pl, lower_T_left)
            lower_hull.extend(left)
        if pmax != pr:
            lower_T_right = self.get_T(pr, pmax, points, True)
            right = self.get_lower_hull(pr, pmax, lower_T_right)
            lower_hull.extend(right)

        return lower_hull

    def fit_set(self, points):
        self.points = points

    def add_point(self, pt):
        self.points.append(pt)

    def compute_hull(self):
        if len(self.points) < 3:
            print("Hull doesn't exist!!")
            return []

        pmin_u = min(self.points)
        pmax_u = max(self.points)
        pmin_l = min(self.points, key=lambda p: (p.x, -p.y))
        pmax_l = max(self.points, key=lambda p: (p.x, -p.y))

        upper_T = self.get_T(pmin_u, pmax_u, self.points, False)
        upper_hull = self.get_upper_hull(pmin_u, pmax_u, upper_T)

        lower_T = self.get_T(pmin_l, pmax_l, self.points, True)
        lower_hull = self.get_lower_hull(pmin_l, pmax_l, lower_T)

        hull_edges = upper_hull + lower_hull
        if pmin_u != pmin_l:
            hull_edges.append(pmin_l)
            hull_edges.append(pmin_u)
        if pmax_u != pmax_l:
            hull_edges.append(pmax_l)
            hull_edges.append(pmax_u)

        hull_edges.sort()
        hull = [hull_edges[0]]
        for i in range(1, len(hull_edges)):
            if hull_edges[i] != hull_edges[i - 1]:
                hull.append(hull_edges[i])

        return hull

# Usage example:
kps = KPS()
points = [Point(0, 3), Point(2, 2), Point(1, 1), Point(2, 1),
          Point(3, 0), Point(0, 0), Point(3, 3)]
kps.fit_set(points)
convex_hull = kps.compute_hull()
print(convex_hull)
