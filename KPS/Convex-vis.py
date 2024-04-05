import math
import matplotlib.pyplot as plt

class Coordinate:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __str__(self):
        return f"({self.x}, {self.y})"

class Polygon:
    def __init__(self, polygon):
        self.polygon = polygon
        self.convex_hull = []

    def cross_product(self, point_a, point_b, point_c):
        return (point_c.x - point_a.x) * (point_b.y - point_a.y) - (point_c.y - point_a.y) * (point_b.x - point_a.x) > 0

    def sort_points(self, point_a, point_b):
        if point_a.x != point_b.x:
            return point_a.x < point_b.x
        return point_a.y < point_b.y

    def compute_convex_hull(self):
        if len(self.polygon) < 3:
            return []

        self.polygon.sort(key=lambda point: (point.x, point.y))
        upper_hull = [self.polygon[0], self.polygon[1]]

        for i in range(2, len(self.polygon)):
            while len(upper_hull) > 1 and not self.cross_product(upper_hull[-2], upper_hull[-1], self.polygon[i]):
                upper_hull.pop()
            upper_hull.append(self.polygon[i])

        lower_hull = [self.polygon[-1], self.polygon[-2]]

        for i in range(len(self.polygon) - 3, -1, -1):
            while len(lower_hull) > 1 and not self.cross_product(lower_hull[-2], lower_hull[-1], self.polygon[i]):
                lower_hull.pop()
            lower_hull.append(self.polygon[i])

        self.convex_hull = upper_hull[:-1] + lower_hull[:-1]
        
        return self.convex_hull

    def plot_convex_hull(self):
        plt.figure(figsize=(8, 6))

        
        polygon_x = [point.x for point in self.polygon]
        polygon_y = [point.y for point in self.polygon]
        plt.scatter(polygon_x, polygon_y, color='blue', label='Polygon Points')

        
        hull_x = [point.x for point in self.convex_hull]
        hull_y = [point.y for point in self.convex_hull]
        hull_x.append(self.convex_hull[0].x)  # Connect back to the first point
        hull_y.append(self.convex_hull[0].y)
        plt.plot(hull_x, hull_y, color='red', linestyle='-', linewidth=2, label='Convex Hull')

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('Convex Hull Visualization')
        plt.legend()
        plt.grid(True)
        plt.show()


polygon_points = [
    Coordinate(3, 4),
    Coordinate(2.5, 3),
    Coordinate(2, 2),
    Coordinate(1.5, 1.5),
    Coordinate(1, 1),
    Coordinate(0.5, 0.5),
    Coordinate(0, 0),
    Coordinate(-0.5, -0.5),
    Coordinate(-1, -1),
    Coordinate(-1.5, -1.5)
]
polygon = Polygon(polygon_points)
convex_hull = polygon.compute_convex_hull()
polygon.plot_convex_hull()
