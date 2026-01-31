#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

using ld = long double;
const ld EPS = 1e-9;
const ld PI = acos(-1.0);

struct Point {
    ld x, y;
    Point operator+(const Point& other) const { return {x + other.x, y + other.y}; }
    Point operator-(const Point& other) const { return {x - other.x, y - other.y}; }
    Point operator*(ld scalar) const { return {x * scalar, y * scalar}; }
    Point operator/(ld scalar) const { return {x / scalar, y / scalar}; }
};

struct Line {
    ld a, b, c; // ax + by + c = 0
};

// 1. Get equation of line from two points
Line lineFromPoints(Point p1, Point p2) {
    ld a = p1.y - p2.y;
    ld b = p2.x - p1.x;
    ld c = -a * p1.x - b * p1.y;
    return {a, b, c};
}

// 2. Intersection of two lines
// Returns {true, point} if they intersect, {false, {0,0}} if parallel
pair<bool, Point> intersect(Line l1, Line l2) {
    ld det = l1.a * l2.b - l2.a * l1.b;
    if (abs(det) < EPS) return {false, {0, 0}}; // Parallel
    ld x = (l1.b * l2.c - l2.b * l1.c) / det;
    ld y = (l1.c * l2.a - l2.c * l1.a) / det;
    return {true, {x, y}};
}

// 3. Area of a triangle from three points
ld triangleArea(Point a, Point b, Point c) {
    return abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) / 2.0;
}

// 4. Perpendicular line to another line going through a given point
Line perpendicularLine(Line l, Point p) {
    // Original: ax + by + c = 0 -> Perpendicular: -bx + ay + d = 0
    ld newA = -l.b;
    ld newB = l.a;
    ld newC = -newA * p.x - newB * p.y;
    return {newA, newB, newC};
}

// 5. Reflection of a point over a line
Point reflectPoint(Point p, Line l) {
    ld dist = (l.a * p.x + l.b * p.y + l.c) / (l.a * l.a + l.b * l.b);
    return {p.x - 2 * l.a * dist, p.y - 2 * l.b * dist};
}

// 6. Distance between two points (Euclidean)
ld dist(Point p1, Point p2) {
    // hypot(dx, dy) is equivalent to sqrt(dx*dx + dy*dy) but safer against overflow
    return hypot(p1.x - p2.x, p1.y - p2.y);
}

struct Point3D {
    ld x, y, z;
    
    Point3D operator+(const Point3D& other) const { return {x + other.x, y + other.y, z + other.z}; }
    Point3D operator-(const Point3D& other) const { return {x - other.x, y - other.y, z - other.z}; }
    Point3D operator*(ld scalar) const { return {x * scalar, y * scalar, z * scalar}; }
    Point3D operator/(ld scalar) const { return {x / scalar, y / scalar, z / scalar}; }

    // Dot Product: u . v
    ld dot(const Point3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Cross Product: u x v (Result is a vector perpendicular to both)
    Point3D cross(const Point3D& other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }

    // Magnitude squared
    ld norm_sq() const {
        return x*x + y*y + z*z;
    }

    // Magnitude (Length)
    ld norm() const {
        return sqrt(norm_sq());
    }
};

struct Line3D {
    Point3D p; // A point on the line
    Point3D v; // Direction vector
};

// 1. Equation of line from two points
// In 3D, we store this as a Point and a Direction Vector
Line3D lineFromPoints(Point3D p1, Point3D p2) {
    return {p1, p2 - p1};
}

// 2. Intersection of two lines
// In 3D, lines are often "Skew" (they don't touch and aren't parallel).
// This function checks if they actually intersect.
pair<bool, Point3D> intersect(Line3D l1, Line3D l2) {
    Point3D p1 = l1.p, v1 = l1.v;
    Point3D p2 = l2.p, v2 = l2.v;

    Point3D p1p2 = p2 - p1;
    
    // Check if lines are coplanar: (p2-p1) . (v1 x v2) == 0
    if (abs(p1p2.dot(v1.cross(v2))) > EPS) {
        return {false, {0,0,0}}; // Skew lines (do not intersect)
    }

    // Check if parallel (cross product is 0)
    Point3D crossProd = v1.cross(v2);
    if (crossProd.norm() < EPS) {
        return {false, {0,0,0}}; // Parallel
    }

    // If coplanar and not parallel, they intersect. solve for t1:
    // P1 + t1*V1 = P2 + t2*V2  =>  t1(V1 x V2) = (P2 - P1) x V2
    // We compare magnitudes or components to find t1.
    // Using vector triple product projection logic:
    ld t1 = p1p2.cross(v2).dot(crossProd) / crossProd.norm_sq();
    
    return {true, p1 + v1 * t1};
}

// 3. Area of a triangle from three points
// Area = 0.5 * |AB x AC|
ld triangleArea(Point3D a, Point3D b, Point3D c) {
    Point3D ab = b - a;
    Point3D ac = c - a;
    return ab.cross(ac).norm() / 2.0;
}

// 4. Perpendicular line to another line going through a given point
Line3D perpendicularLine(Line3D l, Point3D p) {
    // First, find the projection of point p onto line l
    Point3D ap = p - l.p;
    ld t = ap.dot(l.v) / l.v.norm_sq();
    Point3D projection = l.p + l.v * t;

    // The perpendicular line passes through p and the projection point
    return lineFromPoints(p, projection);
}

// 5. Reflection of a point over a line
Point3D reflectPoint(Point3D p, Line3D l) {
    Point3D ap = p - l.p;
    ld t = ap.dot(l.v) / l.v.norm_sq();
    Point3D projection = l.p + l.v * t;
    
    // Reflection formula: R = 2*Projection - Original
    return projection * 2.0 - p;
}

// 6. Distance between two points (Euclidean)
ld dist(Point3D p1, Point3D p2) {
    return (p1 - p2).norm();
}

struct Plane {
    ld a, b, c, d; // Equation: ax + by + cz + d = 0
};

Plane planeFromPoints(Point3D p1, Point3D p2, Point3D p3) {
    // 1. Create two vectors on the plane
    Point3D u = p2 - p1;
    Point3D v = p3 - p1;

    // 2. Cross product gives the Normal Vector (A, B, C)
    Point3D normal = u.cross(v);
    
    // Check if points are collinear (Normal is 0,0,0)
    if (normal.norm() < EPS) {
        // Handle error: Points are collinear, they define a line, not a plane.
        return {0, 0, 0, 0}; 
    }

    // 3. Solve for D
    // ax + by + cz + d = 0  =>  d = -(ax + by + cz)
    // This is equivalent to -(normal . p1)
    ld d = -normal.dot(p1);

    return {normal.x, normal.y, normal.z, d};
}

// Returns {true, Line} if they intersect, {false, ...} if parallel
pair<bool, Line3D> intersectPlanes(Plane p1, Plane p2) {
    Point3D n1 = {p1.a, p1.b, p1.c};
    Point3D n2 = {p2.a, p2.b, p2.c};
    
    // 1. Direction of the intersection line
    Point3D v = n1.cross(n2);

    // If the cross product is nearly zero, the planes are parallel
    if (v.norm() < EPS) return {false, {}};

    // 2. Find a point on the line
    // Solve the system with one variable set to zero. 
    // We choose the most stable variable to set to zero based on the cross product.
    Point3D p;
    if (abs(v.z) > EPS) {
        // Set z = 0, solve:
        // a1*x + b1*y = -d1
        // a2*x + b2*y = -d2
        ld det = p1.a * p2.b - p2.a * p1.b;
        p.x = (p1.b * p2.d - p2.b * p1.d) / det;
        p.y = (p1.d * p2.a - p2.d * p1.a) / det;
        p.z = 0;
    } else if (abs(v.y) > EPS) {
        // Set y = 0
        ld det = p1.a * p2.c - p2.a * p1.c;
        p.x = (p1.c * p2.d - p2.c * p1.d) / det;
        p.z = (p1.d * p2.a - p2.d * p1.a) / det;
        p.y = 0;
    } else {
        // Set x = 0
        ld det = p1.b * p2.c - p2.b * p1.c;
        p.y = (p1.c * p2.d - p2.c * p1.d) / det;
        p.z = (p1.d * p2.b - p2.d * p1.b) / det;
        p.x = 0;
    }

    return {true, {p, v}};
}

// Returns {true, point} if they intersect, {false, {0,0,0}} if parallel
pair<bool, Point3D> intersectLinePlane(Line3D l, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    ld dot_prod = n.dot(l.v);

    // If the direction of the line is perpendicular to the normal, 
    // the line is parallel to the plane.
    if (abs(dot_prod) < EPS) return {false, {0, 0, 0}};

    // Solve for parameter t: a(p.x + t*v.x) + b(p.y + t*v.y) + c(p.z + t*v.z) + d = 0
    ld t = -(n.dot(l.p) + pl.d) / dot_prod;
    
    return {true, l.p + l.v * t};
}

Point3D reflectPointPlane(Point3D p, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    
    // Distance (with sign) from point to plane divided by |n|^2
    ld dist_factor = (n.dot(p) + pl.d) / n.norm_sq();
    
    // Reflection formula: P' = P - 2 * dist_factor * n
    return p - n * (2.0 * dist_factor);
}

// Shortest distance from Point to Plane (Bonus)
ld pointPlaneDist(Point3D p, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    return abs(n.dot(p) + pl.d) / n.norm();
}