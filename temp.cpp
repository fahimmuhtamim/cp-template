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

// Get equation of line from two points
Line lineFromPoints(Point p1, Point p2) {
    ld a = p1.y - p2.y;
    ld b = p2.x - p1.x;
    ld c = -a * p1.x - b * p1.y;
    return {a, b, c};
}

// Returns {true, point} if they intersect, {false, {0,0}} if parallel
pair<bool, Point> intersect(Line l1, Line l2) {
    ld det = l1.a * l2.b - l2.a * l1.b;
    if (abs(det) < EPS) return {false, {0, 0}}; // Parallel
    ld x = (l1.b * l2.c - l2.b * l1.c) / det;
    ld y = (l1.c * l2.a - l2.c * l1.a) / det;
    return {true, {x, y}};
}

ld triangleArea(Point a, Point b, Point c) {
    return abs(a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)) / 2.0;
}

// Perpendicular line to another line going through a given point
Line perpendicularLine(Line l, Point p) {
    // Original: ax + by + c = 0 -> Perpendicular: -bx + ay + d = 0
    ld newA = -l.b;
    ld newB = l.a;
    ld newC = -newA * p.x - newB * p.y;
    return {newA, newB, newC};
}

// Reflection of a point over a line
Point reflectPoint(Point p, Line l) {
    ld dist = (l.a * p.x + l.b * p.y + l.c) / (l.a * l.a + l.b * l.b);
    return {p.x - 2 * l.a * dist, p.y - 2 * l.b * dist};
}

ld dist(Point p1, Point p2) {
    // hypot(dx, dy) is equivalent to sqrt(dx*dx + dy*dy) but safer against overflow
    return hypot(p1.x - p2.x, p1.y - p2.y);
}

ld dist(Point p1, Point p2) {
    return hypot(p1.x - p2.x, p1.y - p2.y);
}

vector<Line> angleBisectors(Line l1, Line l2) {
    ld d1 = hypot(l1.a, l1.b);
    ld d2 = hypot(l2.a, l2.b);
    
    // Equations: (a1/d1 ± a2/d2)x + (b1/d1 ± b2/d2)y + (c1/d1 ± c2/d2) = 0
    Line b1 = {l1.a/d1 + l2.a/d2, l1.b/d1 + l2.b/d2, l1.c/d1 + l2.c/d2};
    Line b2 = {l1.a/d1 - l2.a/d2, l1.b/d1 - l2.b/d2, l1.c/d1 - l2.c/d2};
    
    return {b1, b2};
}

Point getIncenter(Point A, Point B, Point C) {
    ld a = dist(B, C);
    ld b = dist(A, C);
    ld c = dist(A, B);
    return (A * a + B * b + C * c) / (a + b + c);
}

Point getCircumcenter(Point A, Point B, Point C) {
    ld D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    // Check if D is nearly 0 for collinear points
    ld ux = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) + (C.x * C.x + C.y * C.y) * (A.y - B.y)) / D;
    ld uy = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) + (C.x * C.x + C.y * C.y) * (B.x - A.x)) / D;
    return {ux, uy};
}

Point getOrthocenter(Point A, Point B, Point C) {
    Point G = (A + B + C) / 3.0; // Centroid
    Point O = getCircumcenter(A, B, C);
    return (G * 3.0) - (O * 2.0);
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

// Equation of line from two points
Line3D lineFromPoints(Point3D p1, Point3D p2) {
    return {p1, p2 - p1};
}

// Intersection of two lines
pair<bool, Point3D> intersect(Line3D l1, Line3D l2) {
    Point3D p1 = l1.p, v1 = l1.v;
    Point3D p2 = l2.p, v2 = l2.v;

    Point3D p1p2 = p2 - p1;
    
    if (abs(p1p2.dot(v1.cross(v2))) > EPS) {
        return {false, {0,0,0}}; 
    }

    Point3D crossProd = v1.cross(v2);
    if (crossProd.norm() < EPS) {
        return {false, {0,0,0}}; 
    }
    ld t1 = p1p2.cross(v2).dot(crossProd) / crossProd.norm_sq();
    
    return {true, p1 + v1 * t1};
}

ld triangleArea(Point3D a, Point3D b, Point3D c) {
    Point3D ab = b - a;
    Point3D ac = c - a;
    return ab.cross(ac).norm() / 2.0;
}

// Perpendicular line to another line going through a given point
Line3D perpendicularLine(Line3D l, Point3D p) {
    Point3D ap = p - l.p;
    ld t = ap.dot(l.v) / l.v.norm_sq();
    Point3D projection = l.p + l.v * t;

    return lineFromPoints(p, projection);
}

Point3D reflectPoint(Point3D p, Line3D l) {
    Point3D ap = p - l.p;
    ld t = ap.dot(l.v) / l.v.norm_sq();
    Point3D projection = l.p + l.v * t;
    
    return projection * 2.0 - p;
}

ld dist(Point3D p1, Point3D p2) {
    return (p1 - p2).norm();
}

struct Plane {
    ld a, b, c, d; // ax + by + cz + d = 0
};

Plane planeFromPoints(Point3D p1, Point3D p2, Point3D p3) {
    Point3D u = p2 - p1;
    Point3D v = p3 - p1;

    Point3D normal = u.cross(v);
    if (normal.norm() < EPS) {
        return {0, 0, 0, 0}; 
    }
    ld d = -normal.dot(p1);

    return {normal.x, normal.y, normal.z, d};
}

pair<bool, Line3D> intersectPlanes(Plane p1, Plane p2) {
    Point3D n1 = {p1.a, p1.b, p1.c};
    Point3D n2 = {p2.a, p2.b, p2.c};
    
    Point3D v = n1.cross(n2);
    if (v.norm() < EPS) return {false, {}};
    Point3D p;
    if (abs(v.z) > EPS) {
        ld det = p1.a * p2.b - p2.a * p1.b;
        p.x = (p1.b * p2.d - p2.b * p1.d) / det;
        p.y = (p1.d * p2.a - p2.d * p1.a) / det;
        p.z = 0;
    } else if (abs(v.y) > EPS) {
        ld det = p1.a * p2.c - p2.a * p1.c;
        p.x = (p1.c * p2.d - p2.c * p1.d) / det;
        p.z = (p1.d * p2.a - p2.d * p1.a) / det;
        p.y = 0;
    } else {
        ld det = p1.b * p2.c - p2.b * p1.c;
        p.y = (p1.c * p2.d - p2.c * p1.d) / det;
        p.z = (p1.d * p2.b - p2.d * p1.b) / det;
        p.x = 0;
    }

    return {true, {p, v}};
}

pair<bool, Point3D> intersectLinePlane(Line3D l, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    ld dot_prod = n.dot(l.v);
    if (abs(dot_prod) < EPS) return {false, {0, 0, 0}};
    ld t = -(n.dot(l.p) + pl.d) / dot_prod;
    
    return {true, l.p + l.v * t};
}

Point3D reflectPointPlane(Point3D p, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    ld dist_factor = (n.dot(p) + pl.d) / n.norm_sq();

    return p - n * (2.0 * dist_factor);
}

// Shortest distance from Point to Plane (Bonus)
ld pointPlaneDist(Point3D p, Plane pl) {
    Point3D n = {pl.a, pl.b, pl.c};
    return abs(n.dot(p) + pl.d) / n.norm();
}
