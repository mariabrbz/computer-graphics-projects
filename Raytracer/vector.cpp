#include <cmath>

class Vector {
    private:
    double coords[3];

    public:
    explicit Vector(double x = 0., double y = 0., double z = 0.) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    };

    Vector& operator+=(const Vector& b) {
        coords[0] += b[0];
        coords[1] += b[1];
        coords[2] += b[2];
        return *this;
    }

    const double& operator [] (int i) const {return coords[i];}
    double& operator [] (int i) {return coords[i];}
};

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
} 

Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
} 

Vector operator-(const Vector& b) {
    return Vector(-b[0], -b[1], -b[2]);
} 

Vector operator*(double t, const Vector& a) {
    return Vector(a[0]*t, a[1]*t, a[2]*t);
}

Vector operator*(const Vector& a, const Vector& b) {            //define V*U as the componentwise multiplication
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

Vector operator/(const Vector& a, double t) {
    return Vector(a[0]/t, a[1]/t, a[2]/t);
}

double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

double norm(const Vector& a) {
    return sqrt(dot(a, a));
}

double sq_norm(const Vector& a) {
    return dot(a, a);
}