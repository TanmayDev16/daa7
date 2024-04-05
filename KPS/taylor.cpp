#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <bits/stdc++.h>
class Point
{
public:
    double x;
    double y;

    Point(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    bool operator==(const Point &other) const
    {
        return this->x == other.x && this->y == other.y;
    }

    bool operator<(const Point &other) const
    {
        return (this->x < other.x) || (this->x == other.x && this->y < other.y);
    }
};

class KPS
{
private:
    std::vector<Point> points;

    double find_median(std::vector<double> &arr)
    {
        std::sort(arr.begin(), arr.end());
        return arr[arr.size() / 2];
    }

    void swap(std::vector<Point> &arr, int a, int b)
    {
        std::swap(arr[a], arr[b]);
    }

    int partition(std::vector<Point> &arr, int l, int r, double x)
    {
        int i = l;
        for (int j = l; j < r; j++)
        {
            if (arr[j].x == x)
            {
                swap(arr, j, r);
                break;
            }
        }

        i = l;
        for (int j = l; j < r; j++)
        {
            if (arr[j].x <= x)
            {
                swap(arr, i, j);
                i++;
            }
        }
        swap(arr, i, r);
        return i;
    }

    double kth_smallest(std::vector<Point> &arr, int l, int r, int k)
    {
        if (k > 0 && k <= r - l + 1)
        {
            int n = r - l + 1;
            int i = 0;
            std::vector<double> median;
            while (i < n / 5)
            {
                median.push_back(find_median(std::vector<double>(arr.begin() + l + i * 5, arr.begin() + l + i * 5 + 5)));
                i++;
            }
            if (i * 5 < n)
            {
                median.push_back(find_median(std::vector<double>(arr.begin() + l + i * 5, arr.begin() + l + i * 5 + n % 5)));
                i++;
            }

            double med_of_med = (i == 1) ? median[i - 1] : kth_smallest(median, 0, i - 1, i / 2);

            int pos = partition(arr, l, r, med_of_med);

            if (pos - l == k - 1)
            {
                return arr[pos].x;
            }
            if (pos - l > k - 1)
            {
                return kth_smallest(arr, l, pos - 1, k);
            }
            return kth_smallest(arr, pos + 1, r, k - pos + l - 1);
        }

        return INFINITY;
    }

    double abss(double a)
    {
        return std::abs(a);
    }

    std::vector<Point> get_T(Point p1, Point p2, std::vector<Point> &points, bool flag)
    {
        std::vector<Point> upper_T;
        double slope = (p1.y - p2.y) / (p1.x - p2.x);
        for (const auto &curr_point : points)
        {
            if (p1.x < curr_point.x && curr_point.x < p2.x)
            {
                double curr_slope = (p1.y - curr_point.y) / (p1.x - curr_point.x);
                if ((!flag && curr_slope > slope) || (flag && curr_slope < slope))
                {
                    upper_T.push_back(curr_point);
                }
            }
        }
        upper_T.push_back(p1);
        upper_T.push_back(p2);
        return upper_T;
    }

    std::pair<Point, Point> get_upper_bridge(std::vector<Point> &points, double median)
    {
        std::sort(points.begin(), points.end(), [](const Point &p1, const Point &p2)
                  { return p1.x < p2.x; });

        std::vector<Point> candidates;
        std::vector<std::pair<Point, Point>> pairs;
        std::vector<double> slopes;

        if (points.size() % 2 == 0)
        {
            for (int i = 0; i < points.size(); i += 2)
            {
                pairs.push_back(std::make_pair(points[i], points[i + 1]));
            }
        }
        else
        {
            candidates.push_back(points[0]);
            for (int i = 1; i < points.size(); i += 2)
            {
                pairs.push_back(std::make_pair(points[i], points[i + 1]));
            }
        }

        for (const auto &pair : pairs)
        {
            double x1 = pair.first.x;
            double x2 = pair.second.x;
            double y1 = pair.first.y;
            double y2 = pair.second.y;

            if (x1 == x2)
            {
                if (y1 > y2)
                {
                    candidates.push_back(pair.first);
                }
                else
                {
                    candidates.push_back(pair.second);
                }
                slopes.push_back(INFINITY);
            }
            else
            {
                double slope = (y2 - y1) / (x2 - x1);
                slopes.push_back(slope);
            }
        }

        slopes.erase(std::remove(slopes.begin(), slopes.end(), INFINITY), slopes.end());
        double median_slope = kth_smallest(slopes, 0, slopes.size() - 1, (slopes.size() + 1) / 2);

        std::vector<std::pair<Point, Point>> SMALL, EQUAL, LARGE;

        for (int i = 0; i < pairs.size(); i++)
        {
            double x1 = pairs[i].first.x;
            double x2 = pairs[i].second.x;
            double y1 = pairs[i].first.y;
            double y2 = pairs[i].second.y;

            if (x1 != x2)
            {
                double slope = (y2 - y1) / (x2 - x1);
                if (std::abs(slope - median_slope) < 0.001)
                {
                    EQUAL.push_back(pairs[i]);
                }
                else if (slope < median_slope)
                {
                    SMALL.push_back(pairs[i]);
                }
                else if (slope > median_slope)
                {
                    LARGE.push_back(pairs[i]);
                }
            }
        }

        double max_c = -INFINITY;
        for (const auto &point : points)
        {
            double curr_c = (point.y - median_slope * point.x);
            if (curr_c > max_c)
            {
                max_c = curr_c;
            }
        }

        Point pmin(INFINITY, INFINITY);
        Point pmax(-INFINITY, -INFINITY);

        for (const auto &point : points)
        {
            double curr_c = point.y - median_slope * point.x;
            if (std::abs(curr_c - max_c) < 0.001)
            {
                if (point.x < pmin.x)
                {
                    pmin.x = point.x;
                    pmin.y = point.y;
                }
                if (point.x > pmax.x)
                {
                    pmax.x = point.x;
                    pmax.y = point.y;
                }
            }
        }

        if (pmin.x <= median && median < pmax.x)
        {
            return std::make_pair(pmin, pmax);
        }
        else if (pmax.x <= median)
        {
            for (const auto &pt : EQUAL)
            {
                candidates.push_back(pt.second);
            }
            for (const auto &pt : LARGE)
            {
                candidates.push_back(pt.second);
            }
            for (const auto &pt : SMALL)
            {
                candidates.push_back(pt.first);
            }
            return get_upper_bridge(candidates, median);
        }
        else if (pmin.x > median)
        {
            for (const auto &pt : EQUAL)
            {
                candidates.push_back(pt.first);
            }
            for (const auto &pt : LARGE)
            {
                candidates.push_back(pt.first);
            }
            for (const auto &pt : SMALL)
            {
                candidates.push_back(pt.first);
            }
            return get_upper_bridge(candidates, median);
        }
    }

    std::pair<Point, Point> get_lower_bridge(std::vector<Point> &points, double median)
    {
        std::sort(points.begin(), points.end(), [](const Point &p1, const Point &p2)
                  { return p1.x < p2.x; });

        std::vector<Point> candidates;
        std::vector<std::pair<Point, Point>> pairs;
        std::vector<double> slopes;

        if (points.size() % 2 == 0)
        {
            for (int i = 0; i < points.size(); i += 2)
            {
                pairs.push_back(std::make_pair(points[i], points[i + 1]));
            }
        }
        else
        {
            candidates.push_back(points[0]);
            for (int i = 1; i < points.size(); i += 2)
            {
                pairs.push_back(std::make_pair(points[i], points[i + 1]));
            }
        }

        for (const auto &pair : pairs)
        {
            double x1 = pair.first.x;
            double x2 = pair.second.x;
            double y1 = pair.first.y;
            double y2 = pair.second.y;

            if (x1 == x2)
            {
                if (y1 > y2)
                {
                    candidates.push_back(pair.second);
                }
                else
                {
                    candidates.push_back(pair.first);
                }
                slopes.push_back(INFINITY);
            }
            else
            {
                double slope = (y2 - y1) / (x2 - x1);
                slopes.push_back(slope);
            }
        }

        slopes.erase(std::remove(slopes.begin(), slopes.end(), INFINITY), slopes.end());
        double median_slope = kth_smallest(slopes, 0, slopes.size() - 1, (slopes.size() + 1) / 2);

        std::vector<std::pair<Point, Point>> SMALL, EQUAL, LARGE;

        for (int i = 0; i < pairs.size(); i++)
        {
            double x1 = pairs[i].first.x;
            double x2 = pairs[i].second.x;
            double y1 = pairs[i].first.y;
            double y2 = pairs[i].second.y;

            if (x1 != x2)
            {
                double slope = (y2 - y1) / (x2 - x1);
                if (std::abs(slope - median_slope) < 0.001)
                {
                    EQUAL.push_back(pairs[i]);
                }
                else if (slope < median_slope)
                {
                    SMALL.push_back(pairs[i]);
                }
                else if (slope > median_slope)
                {
                    LARGE.push_back(pairs[i]);
                }
            }
        }

        double min_c = INFINITY;

        for (const auto &point : points)
        {
            double curr_c = (point.y - median_slope * point.x);
            if (curr_c < min_c)
            {
                min_c = curr_c;
            }
        }

        Point pmin(INFINITY, INFINITY);
        Point pmax(-INFINITY, -INFINITY);

        for (const auto &point : points)
        {
            double curr_c = point.y - median_slope * point.x;
            if (std::abs(curr_c - min_c) < 0.001)
            {
                if (point.x < pmin.x)
                {
                    pmin.x = point.x;
                    pmin.y = point.y;
                }
                if (point.x > pmax.x)
                {
                    pmax.x = point.x;
                    pmax.y = point.y;
                }
            }
        }

        if (pmin.x <= median && median < pmax.x)
        {
            return std::make_pair(pmin, pmax);
        }
        else if (pmax.x <= median)
        {
            for (const auto &pt : EQUAL)
            {
                candidates.push_back(pt.second);
            }
            for (const auto &pt : LARGE)
            {
                candidates.push_back(pt.second);
            }
            for (const auto &pt : SMALL)
            {
                candidates.push_back(pt.first);
            }
            return get_lower_bridge(candidates, median);
        }
        else if (pmin.x > median)
        {
            for (const auto &pt : EQUAL)
            {
                candidates.push_back(pt.first);
            }
            for (const auto &pt : LARGE)
            {
                candidates.push_back(pt.first);
            }
            for (const auto &pt : SMALL)
            {
                candidates.push_back(pt.first);
            }
            return get_lower_bridge(candidates, median);
        }
    }

    std::vector<Point> get_upper_hull(Point pmin, Point pmax, std::vector<Point> &points)
    {
        std::vector<Point> upper_hull;
        int n = points.size();
        std::vector<double> arr;
        for (const auto &p : points)
        {
            arr.push_back(p.x);
        }
        double median = (n == 1) ? arr[0] : kth_smallest(arr, 0, n - 1, (n + 1) / 2);
        std::pair<Point, Point> upper_bridge = get_upper_bridge(points, median);

        Point pl = upper_bridge.first;
        Point pr = upper_bridge.second;

        if (pl.x > pr.x)
        {
            std::swap(pl, pr);
        }

        upper_hull.push_back(pl);
        upper_hull.push_back(pr);

        if (pmin != pl)
        {
            std::vector<Point> upper_T_left = get_T(pmin, pl, points, false);
            std::vector<Point> left = get_upper_hull(pmin, pl, upper_T_left);
            upper_hull.insert(upper_hull.end(), left.begin(), left.end());
        }

        if (pmax != pr)
        {
            std::vector<Point> upper_T_right = get_T(pr, pmax, points, false);
            std::vector<Point> right = get_upper_hull(pr, pmax, upper_T_right);
            upper_hull.insert(upper_hull.end(), right.begin(), right.end());
        }

        return upper_hull;
    }

    std::vector<Point> get_lower_hull(Point pmin, Point pmax, std::vector<Point> &points)
    {
        std::vector<Point> lower_hull;
        int n = points.size();
        std::vector<double> arr;
        for (const auto &p : points)
        {
            arr.push_back(p.x);
        }
        double median = (n == 1) ? arr[0] : kth_smallest(arr, 0, n - 1, (n + 1) / 2);
        std::pair<Point, Point> lower_bridge = get_lower_bridge(points, median);

        Point pl = lower_bridge.first;
        Point pr = lower_bridge.second;

        if (pl.x > pr.x)
        {
            std::swap(pl, pr);
        }

        lower_hull.push_back(pl);
        lower_hull.push_back(pr);

        if (pmin != pl)
        {
            std::vector<Point> lower_T_left = get_T(pmin, pl, points, true);
            std::vector<Point> left = get_lower_hull(pmin, pl, lower_T_left);
            lower_hull.insert(lower_hull.end(), left.begin(), left.end());
        }
        if (pmax != pr)
        {
            std::vector<Point> lower_T_right = get_T(pr, pmax, points, true);
            std::vector<Point> right = get_lower_hull(pr, pmax, lower_T_right);
            lower_hull.insert(lower_hull.end(), right.begin(), right.end());
        }

        return lower_hull;
    }

public:
    void fit_set(std::vector<Point> &points)
    {
        this->points = points;
    }

    void add_point(Point pt)
    {
        this->points.push_back(pt);
    }

    std::vector<Point> compute_hull()
    {
        if (this->points.size() < 3)
        {
            std::cout << "Hull doesn't exist!!" << std::endl;
            return std::vector<Point>();
        }

        Point pmin_u = *std::min_element(this->points.begin(), this->points.end());
        Point pmax_u = *std::max_element(this->points.begin(), this->points.end());
        Point pmin_l = *std::min_element(this->points.begin(), this->points.end(), [](const Point &p1, const Point &p2)
                                         { return std::make_pair(p1.x, -p1.y) < std::make_pair(p2.x, -p2.y); });
        Point pmax_l = *std::max_element(this->points.begin(), this->points.end(), [](const Point &p1, const Point &p2)
                                         { return std::make_pair(p1.x, -p1.y) < std::make_pair(p2.x, -p2.y); });

        std::vector<Point> upper_T = get_T(pmin_u, pmax_u, this->points, false);
        std::vector<Point> upper_hull = get_upper_hull(pmin_u, pmax_u, upper_T);

        std::vector<Point> lower_T = get_T(pmin_l, pmax_l, this->points, true);
        std::vector<Point> lower_hull = get_lower_hull(pmin_l, pmax_l, lower_T);

        std::vector<Point> hull_edges = upper_hull;
        hull_edges.insert(hull_edges.end(), lower_hull.begin(), lower_hull.end());
        if (pmin_u != pmin_l)
        {
            hull_edges.push_back(pmin_l);
            hull_edges.push_back(pmin_u);
        }
        if (pmax_u != pmax_l)
        {
            hull_edges.push_back(pmax_l);
            hull_edges.push_back(pmax_u);
        }

        std::sort(hull_edges.begin(), hull_edges.end());
        std::vector<Point> hull;
        hull.push_back(hull_edges[0]);
        for (int i = 1; i < hull_edges.size(); i++)
        {
            if (hull_edges[i] != hull_edges[i - 1])
            {
                hull.push_back(hull_edges[i]);
            }
        }

        return hull;
    }
};

int main()
{
    KPS kps;
    std::vector<Point> points = {Point(0, 3), Point(2, 2), Point(1, 1), Point(2, 1),
                                 Point(3, 0), Point(0, 0), Point(3, 3)};
    kps.fit_set(points);
    std::vector<Point> convex_hull = kps.compute_hull();
    for (const auto &point : convex_hull)
    {
        std::cout << "(" << point.x << ", " << point.y << ") ";
    }
    std::cout << std::endl;

    return 0;
}
