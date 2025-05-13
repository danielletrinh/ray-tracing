#include"primitive.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "../helpers/solver.h"

#define ran() ( double( rand() % 32768 ) / 32768 )

const int BEZIER_MAX_DEGREE = 5;
const int Combination[BEZIER_MAX_DEGREE + 1][BEZIER_MAX_DEGREE + 1] =
{	0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0,
    1, 2, 1, 0, 0, 0,
    1, 3, 3, 1, 0, 0,
    1, 4, 6, 4, 1, 0,
    1, 5, 10,10,5, 1
};

const int MAX_COLLIDE_TIMES = 10;
const int MAX_COLLIDE_RANDS = 10;

std::pair<double, double> ExpBlur::GetXY()
{
    double x,y;
    x = ran();
    // x in [0, 1), but with a higher prob to be a small value
    x = pow(2, x)-1;
    y = rand();
    return std::pair<double, double>(x*cos(y),x*sin(y));
}

// ====================================================

Material::Material() {
    color = absor = Color();
    refl = refr = 0;
    diff = spec = 0;
    rindex = 0;
    drefl = 0;
    texture = NULL;
    blur = new ExpBlur();
}

void Material::Input( std::string var , std::stringstream& fin ) {
    if ( var == "color=" ) color.Input( fin );
    if ( var == "absor=" ) absor.Input( fin );
    if ( var == "refl=" ) fin >> refl;
    if ( var == "refr=" ) fin >> refr;
    if ( var == "diff=" ) fin >> diff;
    if ( var == "spec=" ) fin >> spec;
    if ( var == "drefl=" ) fin >> drefl;
    if ( var == "rindex=" ) fin >> rindex;
    if ( var == "texture=" ) {
        std::string file; fin >> file;
        texture = new Bmp;
        texture->Input( file );
    }
    if ( var == "blur=" ) {
        std::string blurname; fin >> blurname;
        if(blurname == "exp")
            blur = new ExpBlur();
    }
}

// ====================================================

Primitive::Primitive() {
    sample = rand();
    material = new Material;
    next = NULL;
}

Primitive::Primitive( const Primitive& primitive ) {
    *this = primitive;
    material = new Material;
    *material = *primitive.material;
}

Primitive::~Primitive() {
    delete material;
}

void Primitive::Input( std::string var , std::stringstream& fin ) {
    material->Input( var , fin );
}

// -----------------------------------------------

Sphere::Sphere() : Primitive() {
    De = Vector3( 0 , 0 , 1 );
    Dc = Vector3( 0 , 1 , 0 );
}

void Sphere::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O=" ) O.Input( fin );
    if ( var == "R=" ) fin >> R;
    if ( var == "De=" ) De.Input( fin );
    if ( var == "Dc=" ) Dc.Input( fin );
    Primitive::Input( var , fin );
}

CollidePrimitive Sphere::Collide( Vector3 ray_O , Vector3 ray_V ) {
    ray_V = ray_V.GetUnitVector();
    Vector3 P = ray_O - O;
    double b = -P.Dot( ray_V );
    double det = b * b - P.Module2() + R * R;
    CollidePrimitive ret;

    if ( det > EPS ) {
        det = sqrt( det );
        double x1 = b - det  , x2 = b + det;

        if ( x2 < EPS ) return ret;
        if ( x1 > EPS ) {
            ret.dist = x1;
            ret.front = true;
        } else {
            ret.dist = x2;
            ret.front = false;
        }
    } else {
        return ret;
    }

    ret.C = ray_O + ray_V * ret.dist;
    ret.N = ( ret.C - O ).GetUnitVector();
    if ( ret.front == false ) ret.N = -ret.N;
    ret.isCollide = true;
    ret.collide_primitive = this;
    return ret;
}

Color Sphere::GetTexture(Vector3 crash_C) {
    Vector3 I = ( crash_C - O ).GetUnitVector();
    double a = acos( -I.Dot( De ) );
    double b = acos( std::min( std::max( I.Dot( Dc ) / sin( a ) , -1.0 ) , 1.0 ) );
    double u = a / PI , v = b / 2 / PI;
    if ( I.Dot( Dc * De ) < 0 ) v = 1 - v;
    return material->texture->GetSmoothColor( u , v );
}

// -----------------------------------------------

void Plane::Input( std::string var , std::stringstream& fin ) {
    if ( var == "N=" ) N.Input( fin );
    if ( var == "R=" ) fin >> R;
    if ( var == "Dx=" ) Dx.Input( fin );
    if ( var == "Dy=" ) Dy.Input( fin );
    Primitive::Input( var , fin );
}

CollidePrimitive Plane::Collide( Vector3 ray_O , Vector3 ray_V ) {
    ray_V = ray_V.GetUnitVector();
    N = N.GetUnitVector();
    double d = N.Dot( ray_V );
    CollidePrimitive ret;
    if ( fabs( d ) < EPS ) return ret;
    double l = ( N * R - ray_O ).Dot( N ) / d;
    if ( l < EPS ) return ret;

    ret.dist = l;
    ret.front = ( d < 0 );
    ret.C = ray_O + ray_V * ret.dist;
    ret.N = ( ret.front ) ? N : -N;
    ret.isCollide = true;
    ret.collide_primitive = this;
    return ret;
}

Color Plane::GetTexture(Vector3 crash_C) {
    double u = crash_C.Dot( Dx ) / Dx.Module2();
    double v = crash_C.Dot( Dy ) / Dy.Module2();
    return material->texture->GetSmoothColor( u , v );
}

// -----------------------------------------------

void Square::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O=" ) O.Input( fin );
    if ( var == "Dx=" ) Dx.Input( fin );
    if ( var == "Dy=" ) Dy.Input( fin );
    Primitive::Input( var , fin );
}

CollidePrimitive Square::Collide( Vector3 ray_O , Vector3 ray_V ) {
    CollidePrimitive ret;

    ray_V = ray_V.GetUnitVector();
    auto N = (Dx * Dy).GetUnitVector();
    double d = N.Dot(ray_V);

    if (fabs(d) < EPS) {
        return ret;
    }

    // solve equation
    double t = (O - ray_O).Dot(N) / d;
    if (t < EPS) {
        return ret;
    }
    auto P = ray_O + ray_V * t;

    // check whether inside square
    double DxLen2 = Dx.Module2();
    double DyLen2 = Dy.Module2();

    double x2 = abs((P - O).Dot(Dx));
    double y2 = abs((P - O).Dot(Dy));
    if (x2 > DxLen2 || y2 > DyLen2) {
        return ret;
    }

    ret.dist = t;
    ret.front = (d < 0);
    ret.C = P;
    ret.N = (ret.front) ? N : -N;
    ret.isCollide = true;
    ret.collide_primitive = this;
    return ret;
}

Color Square::GetTexture(Vector3 crash_C) {
    double u = (crash_C - O).Dot( Dx ) / Dx.Module2() / 2 + 0.5;
    double v = (crash_C - O).Dot( Dy ) / Dy.Module2() / 2 + 0.5;
    return material->texture->GetSmoothColor( u , v );
}

// -----------------------------------------------

void Cylinder::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O1=" ) O1.Input( fin );
    if ( var == "O2=" ) O2.Input( fin );
    if ( var == "R=" ) fin>>R;
    Primitive::Input( var , fin );
}

CollidePrimitive Cylinder::Collide(Vector3 ray_O, Vector3 ray_V) {
    CollidePrimitive ret;

    // Get cylinder axis and height
    Vector3 axis = O2 - O1;
    double height = axis.Module();
    axis = axis.GetUnitVector();

    // Check intersection with cylinder body
    Vector3 AO = ray_O - O1;
    Vector3 AOxAxis = AO * axis;
    Vector3 ray_VxAxis = ray_V * axis;

    double a = ray_VxAxis.Module2();
    double b = 2 * AOxAxis.Dot(ray_VxAxis);
    double c = AOxAxis.Module2() - R * R;

    double delta = b * b - 4 * a * c;
    double t = -1;
    Vector3 P, N;
    bool isSide = true;

    if (delta > EPS) {
        double t1 = (-b - sqrt(delta)) / (2 * a);
        double t2 = (-b + sqrt(delta)) / (2 * a);

        if (t1 > EPS) {
            P = ray_O + ray_V * t1;
            double h = (P - O1).Dot(axis);
            if (h >= 0 && h <= height) {
                t = t1;
            }
        }

        if (t2 > EPS && (t < 0 || t2 < t)) {
            P = ray_O + ray_V * t2;
            double h = (P - O1).Dot(axis);
            if (h >= 0 && h <= height) {
                t = t2;
            }
        }
    }

    // Check intersection with end caps
    double denom1 = ray_V.Dot(axis);
    if (fabs(denom1) > EPS) {
        double t_cap = (O1 - ray_O).Dot(axis) / denom1;
        if (t_cap > EPS && (t < 0 || t_cap < t)) {
            Vector3 P_cap = ray_O + ray_V * t_cap;
            Vector3 v = P_cap - O1;
            if (v.Module2() <= R * R) {
                t = t_cap;
                P = P_cap;
                N = -axis;  // Bottom cap normal points down
                isSide = false;
                ret.front = (denom1 < 0);
            }
        }
    }

    double denom2 = ray_V.Dot(axis);
    if (fabs(denom2) > EPS) {
        double t_cap = (O2 - ray_O).Dot(axis) / denom2;
        if (t_cap > EPS && (t < 0 || t_cap < t)) {
            Vector3 P_cap = ray_O + ray_V * t_cap;
            Vector3 v = P_cap - O2;
            if (v.Module2() <= R * R) {
                t = t_cap;
                P = P_cap;
                N = axis;  // Top cap normal points up
                isSide = false;
                ret.front = (denom2 < 0);
            }
        }
    }

    if (t < 0) return ret;

    if (isSide) {
        Vector3 hitToAxis = P - O1;
        double h = hitToAxis.Dot(axis);
        Vector3 axisPoint = O1 + axis * h;
        N = (P - axisPoint).GetUnitVector();
        ret.front = (ray_V.Dot(N) < 0);
    }

    ret.dist = t;
    ret.C = P;
    ret.N = N;
    ret.isCollide = true;
    ret.collide_primitive = this;
    return ret;
}

Color Cylinder::GetTexture(Vector3 crash_C) {
    // Calculate texture coordinates
    Vector3 axis = (O2 - O1).GetUnitVector();
    double h = (crash_C - O1).Dot(axis);
    double theta = atan2((crash_C - O1 - axis * h).Dot(axis.GetAnVerticalVector()),
                        (crash_C - O1 - axis * h).Dot(axis * axis.GetAnVerticalVector()));

    double u = theta / (2 * PI);
    double v = h / (O2 - O1).Module();

    return material->texture->GetSmoothColor(u, v);
}

// -----------------------------------------------

void Bezier::Input( std::string var , std::stringstream& fin ) {
    if ( var == "O1=" ) O1.Input( fin );
    if ( var == "O2=" ) O2.Input( fin );
    if ( var == "P=" ) {
        double newR, newZ;
        fin>>newZ>>newR;
        R.push_back(newR);
        Z.push_back(newZ);
        degree++;
    }
    if ( var == "Cylinder" ) {
        double maxR = 0;
        for (int i = 0; i < R.size(); i++) {
            if (R[i] > maxR) {
                maxR = R[i];
            }
        }
        boundingCylinder = new Cylinder(O1, O2, maxR);
        N = (O2 - O1).GetUnitVector();
        Nx = N.GetAnVerticalVector().GetUnitVector();
        Ny = (N * Nx).GetUnitVector();
    }
    Primitive::Input( var , fin );
}

CollidePrimitive Bezier::Collide(Vector3 ray_O, Vector3 ray_V) {
    CollidePrimitive ret;

    // Basic validation
    if (degree < 0 || R.empty() || Z.empty() || !boundingCylinder) return ret;

    // Normalize ray direction
    ray_V = ray_V.GetUnitVector();

    // First check with bounding cylinder for early rejection
    CollidePrimitive cyl_col = boundingCylinder->Collide(ray_O, ray_V);
    if (!cyl_col.isCollide) return ret;

    // Get axis information
    Vector3 axis = O2 - O1;
    double axisLen = axis.Module();
    if (axisLen < EPS) return ret;
    Vector3 axisDir = axis.GetUnitVector();

    // Binary search for intersection
    const int MAX_STEPS = 50;
    double t_min = 0;
    double t_max = cyl_col.dist;
    double best_t = -1;
    double best_dist = 1e6;
    double best_u = 0.5;

    for (int step = 0; step < MAX_STEPS; step++) {
        double t = (t_min + t_max) * 0.5;
        Vector3 P = ray_O + ray_V * t;

        // Project point onto axis
        double s = (P - O1).Dot(axisDir);
        double u = s / axisLen;
        u = std::max(0.0, std::min(1.0, u)); // Clamp to valid range

        // Get radius at current position
        auto [z, r] = valueAt(u);

        // Calculate distance from point to surface of revolution
        Vector3 axisPoint = O1 + axisDir * s;
        Vector3 radialVec = P - axisPoint;
        double radialDist = radialVec.Module();
        double dist = fabs(radialDist - r);

        if (dist < best_dist) {
            best_dist = dist;
            best_t = t;
            best_u = u;
        }

        // Update search range
        if (radialDist > r) {
            t_min = t;  // Ray point is outside surface
        } else {
            t_max = t;  // Ray point is inside surface
        }

        // Exit early if we're close enough
        if (dist < 1e-6) break;
    }

    // If we found a good intersection
    if (best_t > 0 && best_dist < 0.001) {
        Vector3 P = ray_O + ray_V * best_t;
        double u = best_u;

        // Calculate normal at intersection point
        Vector3 axisPoint = O1 + axisDir * (u * axisLen);
        Vector3 radialVec = P - axisPoint;

        // Handle point on axis
        Vector3 radialDir;
        if (radialVec.Module() < EPS) {
            radialDir = axisDir.GetAnVerticalVector().GetUnitVector();
        } else {
            radialDir = radialVec.GetUnitVector();
        }

        // Calculate derivative of radius along curve
        const double du = 0.001;
        auto [z1, r1] = valueAt(std::max(0.0, u - du));
        auto [z2, r2] = valueAt(std::min(1.0, u + du));
        double dr_du = (r2 - r1) / (2 * du);

        // Calculate surface normal
        Vector3 profileTangent = axisDir + radialDir * dr_du;
        profileTangent = profileTangent.GetUnitVector();
        Vector3 circleTangent = (radialDir * axisDir).GetUnitVector();
        ret.N = (profileTangent * circleTangent).GetUnitVector();

        // Ensure normal points outward
        if (radialVec.Dot(ret.N) < 0) ret.N = -ret.N;

        ret.dist = best_t;
        ret.C = P;
        ret.isCollide = true;
        ret.front = (ray_V.Dot(ret.N) < 0);
        ret.collide_primitive = this;
    }

    return ret;
}

Color Bezier::GetTexture(Vector3 crash_C) {
    Vector3 axis = (O2 - O1).GetUnitVector();

    // Calculate height along axis (v coordinate)
    double height = (crash_C - O1).Dot(axis);
    double v = height / (O2 - O1).Module();
    v = std::max(0.0, std::min(1.0, v)); // Clamp to [0,1]

    // Calculate angle around axis (u coordinate)
    Vector3 radial = crash_C - (O1 + axis * height);

    // Ensure consistent coordinate system
    if (Nx.Module() < 0.001 || Ny.Module() < 0.001) {
        // Recalculate coordinate system if needed
        Nx = axis.GetAnVerticalVector().GetUnitVector();
        Ny = (axis * Nx).GetUnitVector();
    }

    // Calculate angle using consistent axes
    double cosAngle = radial.Dot(Nx) / radial.Module();
    double sinAngle = radial.Dot(Ny) / radial.Module();

    // Map angle to [0,1] range without discontinuity
    double u = atan2(sinAngle, cosAngle) / (2 * PI);
    if (u < 0) u += 1.0; // Ensure u is in [0,1]

    return material->texture->GetSmoothColor(u, v);
}

std::pair<double, double> Bezier::valueAt(double u)
{
    if (Z.empty() || R.empty()) {
        return std::make_pair(0.0, 0.0);
    }
    return valueAt(u, Z, R);
}


std::pair<double, double> Bezier::valueAt(double u, const std::vector<double>& xs, const std::vector<double>& ys)
{
    if (xs.empty() || ys.empty() || xs.size() != ys.size()) {
        return std::make_pair(0.0, 0.0);
    }

    const int deg = xs.size() - 1;
    if (deg < 0 || deg >= BEZIER_MAX_DEGREE) {
        return std::make_pair(xs[0], ys[0]); // Default to first point
    }

    // Clamp u to [0,1]
    u = std::max(0.0, std::min(1.0, u));

    double x = 0;
    double y = 0;

    for (int i = 0; i <= deg; i++) {
        double bernstein = Combination[deg][i] * pow(u, i) * pow(1 - u, deg - i);
        x += bernstein * xs[i];
        y += bernstein * ys[i];
    }

    // Ensure radius is non-negative
    y = std::max(0.0, y);

    return std::make_pair(x, y);
}