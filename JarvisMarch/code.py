# To find orientation of ordered triplet (p, q, r).
# The function returns following values
# 0 --> p, q and r are collinear
# 1 --> Clockwise
# 2 --> Counterclockwise
import matplotlib.pyplot as plt
def angle(p, q, r):
    val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])

    if val == 0:
        return 0  # collinear
    return 1 if val > 0 else 2  # clock or counterclock wise

# Prints convex hull of a set of n points.
def convexHull(points):
    n = len(points)
    # There must be at least 3 points
    if n < 3:
        return

    # Initialize Result
    hull = []

    # Find the leftmost point
    l = 0
    for i in range(1, n):
        if points[i][0] < points[l][0]:
            l = i

    # Start from leftmost point, keep moving counterclockwise
    # until reach the start point again. This loop runs O(h)
    # times where h is number of points in result or output.
    p = l
    while True:
        # Add current point to result
        hull.append(points[p])

        # Search for a point 'q' such that orientation(p, q, x)
        # is counterclockwise for all points 'x'. The idea
        # is to keep track of last visited most counterclockwise
        # point in q. If any point 'i' is more counterclockwise
        # than q, then update q.
        q = (p + 1) % n
        for i in range(n):
            # If i is more counterclockwise than current q, then
            # update q
            if angle(points[p], points[i], points[q]) == 2:
                q = i

        # Now q is the most counterclockwise with respect to p
        # Set p as q for next iteration, so that q is added to
        # result 'hull'
        p = q

        # While we don't come to first point
        if p == l:
            break

    # Print Result
    return hull

# Driver program to test above functions
if __name__ == "__main__":
    points = [(0, 0), (1, -4), (-1, -5), (-5, -3), (-3, -1), (-1, -3), (-2, -2), (-1, -1), (-2,-1), (-1, 1)]
    x_1 = [point[0] for point in points]
    y_1 = [point[1] for point in points]
    plt.scatter(x_1, y_1, color='blue', label='Polygon Points')
    hull = convexHull(points)
    x = [point[0] for point in hull]
    y = [point[1] for point in hull]
    
    plt.plot(x, y, color='red', linestyle='-', linewidth=2, label='Convex Hull')
    plt.show()